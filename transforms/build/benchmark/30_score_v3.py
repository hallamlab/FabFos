"""The reference scorer, v3 — the executable half of the contract.

Reads an answer key Y and an implementation's `observations.tsv`, emits three
outputs:

    directional specificity   SCORED   (§5.1) — target vs own off-target, signed by Y
    potency                   SCORED   (§5.2) — diagonal vs enumerated null
    essentiality diagnostic   REPORTED (§5.3) — does effect separate essential genes

WHAT CHANGED FROM v1, AND WHY THIS FILE EXISTS
----------------------------------------------
v1's scorer imported `Y_DIR` from `_bench.py`, which pinned

    BENCH = DATA / "ecspr/benchmark/v1"

so the scorer was WELDED to one benchmark version: pointing it at another key
meant editing a shared module that every builder also imported. Here the key
directory is a CLI argument, `--key-dir`, defaulting to `canon.BENCH_V3_Y`.
Everything else about the measurement is unchanged — this is a re-pointing, not
a re-definition, so v1 and v3 numbers remain comparable.

The default resolves through `canon` exactly as `_v3bench.py` does (repo `src/`
on `sys.path`, then `from fabfos import canon`). `canon.BENCH_V3_Y` raises
CanonError while the key is unbuilt, and that refusal is intended: scoring
against an absent key must never quietly produce an empty result. It is
resolved lazily, at `main()` time and only when `--key-dir` was not given, so
merely importing this module never trips it.

RANK-ONLY BY CONSTRUCTION
-------------------------
Every number is a Mann-Whitney / van-Elteren AUC or a null percentile, so it
depends only on the ORDER of `effect`. `--selftest` proves it: it re-scores under
a family of monotone maps and asserts bit-identical output. Any scalar an
implementation emits — ΔI_eff, a p-value, a rank — is read the same way.

DIRECTION COMES FROM Y, NOT FROM THE ARM
----------------------------------------
The incumbent signed every score with `sgn = -1 if arm.startswith('lof') else 1`.
Here the sign is `dir` on the Y row, so the same code path scores a GOF gain and
a LOF loss without knowing which arm it is in. A LOF target is expected to go
DOWN, so its "correctly directed" effect is `-effect`; a GOF target UP, so it is
`+effect`. Specificity then asks: is a condition's correctly-directed target
effect ranked above its own off-target effects?

NO NULL POOL IN v3
------------------
v3 ships no enumerated null, so potency — a strictly vs-random claim — is
UNDEFINED, by design rather than by omission. `--null` stays optional and the
potency table reports `status = undefined:no_null_pool` with `auc = NaN`. It is
never a zero and never an empty frame: a silent 0.0 there would read as "no
potency", which is a scientific claim we have no evidence for.

SELF-CONTAINED TREE
-------------------
`_panelstats.py` sits beside this file as a VERBATIM copy of
`fig-validate/main/ecspr/validation/_panelstats.py` (stratified_auc,
mannwhitney_u, cluster_bootstrap, null_percentile) and is imported locally, so
the v3 tree reaches sideways into no sibling scope.

Usage:
  mamba run -n ml python 30_score_v3.py --obs <observations.tsv> [--key-dir DIR]
                                        [--null NULL.tsv] [--mode signed|unsigned]
  mamba run -n ml python 30_score_v3.py --selftest
Env: ml
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import _panelstats as ps  # noqa: E402  (verbatim copy, beside this file)

# The build tree this scorer writes into by default. Outputs go into the REPO,
# not into the data library, for the same reason `_v3bench.BUILD_OUT` does: the
# library is assembled from declarations by the staging script, and a file
# written into it directly is one no record describes and no hash covers. v1
# defaulted to `Y_DIR.parent / "scored"`, which put results next to the key —
# fine when the key lived in a scratch tree, wrong now that it lives in the
# library. Pass `--out` to override.
BUILD_OUT = HERE / "v3_build"

ELEMENTS = ["C", "N", "S", "P"]
FACETS = ["netA_iML1515", "netA_iECDH10B", "netB_iML1515", "netB_iECDH10B"]
DIR_SIGN = {"+": 1.0, "-": -1.0, "0": 0.0}

# The columns each Y table must carry. v3's `conditions.tsv` adds `target_basis`,
# `obs_id` and `host_gem_is_proxy` beyond v1's set, and future keys may add more,
# so this is a REQUIRED-SUBSET check, never an equality check and never a
# positional one: every read below is by name. An unexpected extra column is a
# key that carries more provenance than the scorer needs, which is not an error.
# A missing REQUIRED column, though, is silent nonsense downstream — an absent
# `expected_dir` would make every condition's axis 0.0 and score every arm at
# chance — so it is named and refused up front.
REQUIRED_COLS = {
    "conditions.tsv": {"condition_id", "arm", "tier", "element", "expected_dir",
                       "target_mnxm", "essential_on_glucose_minimal"},
    "panel_edges.tsv": {"edge_id", "element"},
    "expectations.tsv": {"condition_id", "edge_id", "role", "dir"},
}


# ---------------------------------------------------------------------------
# loading + the default rule
# ---------------------------------------------------------------------------
def _read_key_table(key_dir: Path, name: str) -> pd.DataFrame:
    """Read one Y table as strings, checking the required columns are present."""
    path = key_dir / name
    if not path.exists():
        raise SystemExit(f"answer key incomplete: {path} does not exist")
    df = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    missing = REQUIRED_COLS[name] - set(df.columns)
    if missing:
        raise SystemExit(
            f"{path}: missing required column(s) {sorted(missing)}. "
            f"Present: {sorted(df.columns)}")
    return df


def load_Y(key_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    cond = _read_key_table(key_dir, "conditions.tsv")
    panel = _read_key_table(key_dir, "panel_edges.tsv")
    exp = _read_key_table(key_dir, "expectations.tsv")
    return cond, panel, exp


def expectation_lookup(exp: pd.DataFrame) -> dict:
    """{(condition_id, edge_id): (role, dir)}. Absent pairs default via `resolve`."""
    return {(r.condition_id, r.edge_id): (r.role, r["dir"])
            for _, r in exp.iterrows()}


def resolve(lookup: dict, cid: str, eid: str) -> tuple[str, str]:
    """The declared default rule: absent (condition, edge) is off_target / 0."""
    return lookup.get((cid, eid), ("off_target", "0"))


# ---------------------------------------------------------------------------
# §5.1 directional specificity
# ---------------------------------------------------------------------------
def directional_specificity(obs: pd.DataFrame, cond: pd.DataFrame, panel: pd.DataFrame,
                            lookup: dict, facet: str, mode: str, element: str | None):
    """Per condition: does the target move furthest along its OWN expected axis?

    The comparison must be apples-to-apples. In signed mode BOTH the target and
    the off-target cells are projected onto the condition's expected target
    direction — `s_c * effect`, where `s_c` is +1 for a GOF condition (target
    should rise) and -1 for a LOF one (target should fall). The question is then
    "is the target's movement-in-the-expected-direction larger than the
    off-targets' movement-along-that-same-axis". Signing only the target — and
    comparing it to raw off-targets — gives any consistently-signed signal
    (e.g. a hub-toucher) a free target advantage, which is a scorer artifact, not
    specificity. In unsigned mode every cell is |effect|, so direction is untested.

    Returns (auc, ci_lo, ci_hi, n_conditions, strata_by_condition).
    """
    o = obs[obs.facet == facet]
    if element:
        edges = set(panel[panel.element == element].edge_id)
        o = o[o.edge_id.isin(edges)]
    cond_dir = dict(zip(cond.condition_id, cond.expected_dir))

    # Group by condition ONCE — scanning all cells per condition is O(N*conds).
    strata_by_cond: dict[str, list] = {}
    for cid, grp in o.groupby("condition_id"):
        s_c = DIR_SIGN.get(cond_dir.get(cid, "0"), 0.0)   # the condition's own axis
        pos, neg = [], []
        for eid, v in zip(grp.edge_id, grp.effect.astype(float)):
            role, _d = resolve(lookup, cid, eid)
            proj = abs(v) if mode == "unsigned" else s_c * v
            (pos if role == "target" else neg).append(proj)
        if pos and neg:
            strata_by_cond[cid] = [(np.array(pos), np.array(neg))]

    if not strata_by_cond:
        return float("nan"), float("nan"), float("nan"), 0, {}
    strata = [s[0] for s in strata_by_cond.values()]
    auc = ps.stratified_auc(strata)
    lo, hi = ps.cluster_bootstrap(strata_by_cond) if len(strata_by_cond) > 1 \
        else (float("nan"), float("nan"))
    return auc, lo, hi, len(strata_by_cond), strata_by_cond


# ---------------------------------------------------------------------------
# §5.2 potency — diagonal vs the enumerated null pool, per column
# ---------------------------------------------------------------------------
def potency(obs: pd.DataFrame, null: pd.DataFrame | None, facet: str,
            element: str | None, panel: pd.DataFrame):
    """Target diagonal against the null pool on the same column (anchor×product).

    Returns (auc, n_columns, status). The status string is the point of this
    signature: potency has THREE distinct undefined-ness modes and collapsing
    them all to a bare NaN loses which one you hit.

      ok                      — measured against a real pool
      undefined:no_null_pool  — no `--null` was given at all. v3's normal case:
                                the key ships no enumerated null, so a vs-random
                                claim cannot be made. NOT a score of zero, NOT a
                                score of 0.5 — no evidence either way.
      undefined:no_overlap    — a null WAS given but shares no scorable column
                                with the observations on this facet/element.
                                That is a mismatch worth seeing, and it must not
                                be confused with the case above.
    """
    if null is None:
        return float("nan"), 0, "undefined:no_null_pool"
    o = obs[obs.facet == facet]
    n = null[null.facet == facet]
    if element:
        edges = set(panel[panel.element == element].edge_id)
        o = o[o.edge_id.isin(edges)]
        n = n[n.edge_id.isin(edges)]
    # Pre-group the null ONCE. Filtering `n[n.edge_id == eid]` inside the loop is
    # an O(|null|) scan per column — catastrophic on a 20M-row null.
    npool = {eid: g.effect.astype(float).values for eid, g in n.groupby("edge_id")}
    strata = []
    for eid, grp in o.groupby("edge_id"):
        pool = npool.get(eid)
        pos = grp.effect.astype(float).values
        if pos.size and pool is not None and pool.size:
            strata.append((pos, pool))
    if not strata:
        return float("nan"), 0, "undefined:no_overlap"
    return ps.stratified_auc(strata), len(strata), "ok"


# ---------------------------------------------------------------------------
# §5.3 the essentiality diagnostic — REPORTED, never gated
# ---------------------------------------------------------------------------
def essentiality_diagnostic(obs: pd.DataFrame, cond: pd.DataFrame, facet: str):
    """Does |effect| on the auxotrophy edge separate essential from non-essential?

    Within element and within facet, never pooled. GPR-silent conditions MUST
    have been emitted as effect=0 (contract §3); if they were dropped, the
    non-essential silents vanish and this measures nothing — so we report how
    many of each class actually carried an observation.

    Only `yes`/`no` rows enter. v3's key also carries `uncertain` (and empty on
    the whole GOF arm); those are excluded rather than folded into `no`, because
    a guessed label here would move a reported diagnostic without any evidence.
    """
    lof = cond[(cond.arm == "lof") & (cond.essential_on_glucose_minimal.isin(["yes", "no"]))]
    o = obs[obs.facet == facet]
    # each LOF condition's strongest movement on ITS OWN target product
    best = {}
    tgt = {r.condition_id: r.target_mnxm for _, r in lof.iterrows()}
    for cid, grp in o[o.condition_id.isin(tgt)].groupby("condition_id"):
        best[cid] = float(grp.effect.astype(float).abs().max())

    rows = []
    for el, g in lof.groupby("element"):
        ess = [best[c] for c in g[g.essential_on_glucose_minimal == "yes"].condition_id if c in best]
        non = [best[c] for c in g[g.essential_on_glucose_minimal == "no"].condition_id if c in best]
        if ess and non:
            auc = ps.auc(np.array(ess), np.array(non))
        else:
            auc = float("nan")
        rows.append(dict(element=el, facet=facet, auc_essential_vs_not=auc,
                         n_essential=len(ess), n_nonessential=len(non),
                         n_essential_total=int((g.essential_on_glucose_minimal == "yes").sum()),
                         n_nonessential_total=int((g.essential_on_glucose_minimal == "no").sum())))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# driver
# ---------------------------------------------------------------------------
def score(obs: pd.DataFrame, null: pd.DataFrame | None, mode: str,
          key_dir: Path) -> dict:
    cond, panel, exp = load_Y(key_dir)
    lookup = expectation_lookup(exp)
    # `tier` may legitimately carry no 'reference' rows at all (v3's key is
    # entirely `scored`); the filter is written to survive that rather than to
    # assume the tier vocabulary of any one key version.
    scored_elements = [e for e in ELEMENTS
                       if e in set(cond[cond.tier != "reference"].element)]

    spec_rows, pot_rows, ess_frames = [], [], []
    for facet in FACETS:
        if facet not in set(obs.facet):
            continue
        for el in [None] + scored_elements:
            auc, lo, hi, n, _ = directional_specificity(
                obs, cond, panel, lookup, facet, mode, el)
            spec_rows.append(dict(facet=facet, element=el or "ALL", mode=mode,
                                  auc=auc, ci_lo=lo, ci_hi=hi, n_conditions=n,
                                  passes=(not np.isnan(lo)) and lo > 0.5))
            pauc, pn, pstatus = potency(obs, null, facet, el, panel)
            pot_rows.append(dict(facet=facet, element=el or "ALL",
                                 auc=pauc, n_columns=pn, status=pstatus))
        ess_frames.append(essentiality_diagnostic(obs, cond, facet))

    return dict(
        specificity=pd.DataFrame(spec_rows),
        potency=pd.DataFrame(pot_rows),
        essentiality=pd.concat(ess_frames, ignore_index=True) if ess_frames else pd.DataFrame(),
    )


def _selftest(key_dir: Path) -> int:
    """Rank-only: any strictly monotone transform leaves every number unchanged."""
    ps._selftest()
    rng = np.random.default_rng(0)
    cond, panel, exp = load_Y(key_dir)
    lookup = expectation_lookup(exp)

    # synthetic observations on one facet: a perfect oracle
    facet = "netA_iML1515"
    pe = panel[panel.element == "C"]
    scond = cond[(cond.arm == "lof") & (cond.element == "C")].head(20)
    assert len(scond) and len(pe), \
        f"{key_dir}: no LOF/C conditions or no C panel edges — cannot self-test"
    rows = []
    for _, c in scond.iterrows():
        for _, e in pe.iterrows():
            role, d = resolve(lookup, c.condition_id, e.edge_id)
            base = 5.0 if role == "target" else rng.uniform(0, 0.5)
            # LOF target should go down; oracle makes it strongly negative
            val = -base if (role == "target") else (rng.uniform(-0.2, 0.2))
            rows.append(dict(facet=facet, condition_id=c.condition_id,
                             edge_id=e.edge_id, effect=val))
    obs = pd.DataFrame(rows)

    def spec_all(o):
        a, *_ = directional_specificity(o, cond, panel, lookup, facet, "signed", "C")
        return a

    # Signed-mode invariance holds under SIGN-PRESERVING strictly-monotone maps
    # only — those with f(0)=0 and f increasing (x^3, arcsinh, scaling, x*|x|).
    # exp/log are NOT sign-preserving: they map every value positive and thereby
    # change the answer to "did the target go DOWN", which is real information,
    # not a nuisance. The contract states this precisely (§3). Unsigned mode
    # reads |effect|, so there ANY positive monotone map is invariant — tested
    # separately below.
    base = spec_all(obs)
    signed_maps = [("x^3", lambda x: x**3),
                   ("5x", lambda x: 5.0 * x),
                   ("arcsinh", lambda x: np.arcsinh(x)),
                   ("x|x|", lambda x: x * np.abs(x))]
    for name, fn in signed_maps:
        o2 = obs.copy()
        o2["effect"] = fn(obs.effect.values)
        got = spec_all(o2)
        assert abs(got - base) < 1e-12 or (np.isnan(got) and np.isnan(base)), \
            f"NOT rank-invariant under sign-preserving {name}: {base} vs {got}"
        print(f"  signed rank-invariant under {name}: {got:.6f} == {base:.6f}")

    # unsigned mode: |effect|, so exp/log (positive monotone) are invariant too
    def spec_unsigned(o):
        a, *_ = directional_specificity(o, cond, panel, lookup, facet, "unsigned", "C")
        return a
    # unsigned reads |effect|, so the invariant transforms are monotone in the
    # MAGNITUDE: exp(|x|), log(|x|+1). (abs(exp(x)) is NOT monotone in |x|, since
    # exp is not even — an implementation that emits a signed scalar and then
    # exponentiates it has changed its magnitude ordering, which is a real change.)
    ub = spec_unsigned(obs)
    for name, fn in [("exp|x|", lambda x: np.exp(np.abs(x))),
                     ("log|x|+1", lambda x: np.sign(x) * np.log(np.abs(x) + 1.0))]:
        o2 = obs.copy()
        o2["effect"] = fn(obs.effect.values)
        got = spec_unsigned(o2)
        assert abs(got - ub) < 1e-12 or (np.isnan(got) and np.isnan(ub)), \
            f"NOT rank-invariant under magnitude-monotone {name} in unsigned mode: {ub} vs {got}"
        print(f"  unsigned rank-invariant under {name}: {got:.6f} == {ub:.6f}")

    # the oracle must pass directional specificity
    assert base > 0.99, f"perfect oracle should score ~1.0 on specificity, got {base}"
    print(f"  perfect oracle specificity AUC = {base:.4f}  (expect ~1.0)")

    # The no-null path is part of the contract, not an edge case — assert it
    # reports undefined rather than crashing or returning a number.
    pauc, pn, pstatus = potency(obs, None, facet, "C", panel)
    assert np.isnan(pauc) and pn == 0 and pstatus == "undefined:no_null_pool", \
        f"no-null potency must be explicitly undefined, got ({pauc}, {pn}, {pstatus})"
    print(f"  potency with no null pool: auc=NaN status={pstatus}")

    print(f"  key: {key_dir}")
    print("SELFTEST PASSED")
    return 0


def _default_key_dir() -> Path:
    """`canon.BENCH_V3_Y`, resolved LAZILY.

    Deferred to call time on purpose. canon raises CanonError for a declared but
    absent symbol, and resolving at import (or in an argparse `default=`) would
    make that failure fire for `--help` and for an explicit `--key-dir` pointing
    somewhere perfectly valid. `_v3bench.py` puts the repo `src/` on sys.path
    the same way; this file does it locally so it does not depend on _v3bench.
    """
    repo = HERE.parents[2]          # transforms/build/benchmark -> repo root
    src = str(repo / "src")
    if src not in sys.path:
        sys.path.insert(0, src)
    from fabfos import canon        # noqa: PLC0415 (deliberately late)
    return Path(canon.BENCH_V3_Y)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--obs", type=Path, help="observations.tsv")
    ap.add_argument("--key-dir", type=Path, default=None,
                    help="answer-key directory (conditions/panel_edges/"
                         "expectations .tsv). Default: canon.BENCH_V3_Y")
    ap.add_argument("--null", type=Path, help="null pool tsv (facet,edge_id,effect)")
    ap.add_argument("--mode", choices=["signed", "unsigned"], default="signed")
    ap.add_argument("--out", type=Path, default=None,
                    help=f"output directory. Default: {BUILD_OUT / 'scored'}")
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()

    key_dir = args.key_dir if args.key_dir is not None else _default_key_dir()
    out = args.out if args.out is not None else BUILD_OUT / "scored"

    if args.selftest:
        return _selftest(key_dir)
    if not args.obs:
        ap.error("--obs is required unless --selftest")

    obs = pd.read_csv(args.obs, sep="\t", dtype={"facet": str, "condition_id": str,
                                                 "edge_id": str, "effect": float})
    null = pd.read_csv(args.null, sep="\t") if args.null else None
    if null is None:
        # Said once, loudly, at the top — a reader who skips the potency table's
        # status column should still not mistake a NaN row for a failed run.
        print("no --null given: potency is UNDEFINED (no null pool), "
              "reported as status=undefined:no_null_pool, not as a score.")
    res = score(obs, null, args.mode, key_dir)

    out.mkdir(parents=True, exist_ok=True)
    for name, df in res.items():
        df.to_csv(out / f"{name}_{args.mode}.tsv", sep="\t", index=False)
        print(f"\n=== {name} ({args.mode}) ===")
        print(df.to_string(index=False))
    print(f"\nkey   <- {key_dir}")
    print(f"wrote -> {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

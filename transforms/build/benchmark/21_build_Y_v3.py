"""Freeze v3's Y -- the answer key: conditions and expectations.

    conditions.tsv     one row per (perturbation, element): what is done, to
                       which host, which product should move, and the stratifiers
    expectations.tsv   sparse (condition, edge) -> role x dir

THE GOVERNING RULE: Y DECLARES BIOLOGY, NEVER LIVENESS
------------------------------------------------------
Every row is derived from a measured observation. Nothing in this builder's
input path reads ECSPr output -- no `out/`, no scored matrix, no
`analysis_report.json`. `50_audit_v3.py` proves that mechanically. A Y derived
from the implementation cannot fail, and an expectation set to "whatever the
incumbent returns" is not an expectation.

WHERE v3's DIRECTION COMES FROM
-------------------------------
`measured`, NOT `expected`. The `expected` column is populated on only 10 of the
382 GOF rows; `measured` is populated on 325. up -> `+`, down -> `-`, anything
else -> NO DIRECTION AND NO TARGET CELL, reported as coverage and never guessed.

The resulting distribution is lopsided and it belongs here rather than being
discovered at score time: 322 up / 3 down / 57 unknown. With the 166 uniformly
negative LOF conditions that is 491 of 548 carrying a usable direction. Nearly
all sign contrast comes from the LOF arm plus 3 GOF rows, so the GOF arm is
close to a pure "does this raise conductance toward its own target more than
elsewhere" test. That is a property of the observation set, not a defect, and
the result has to say so.

HOW A CONDITION GETS ITS TARGETS
--------------------------------
Two sources, in order of authority:

  1. MECHANICAL -- the direct non-currency products of the reactions the
     condition INSERTS, read from X. This is v1's own rule and it is what
     carries the arm: 380 of 382 GOF rows have a perturbed reaction in X.
  2. RESOLVED -- the free-text target string mapped to an id in
     target_resolution.tsv. Enrichment; partial by nature, reported as coverage.

v1 took its GOF targets from `_gof_specs.GOF_SPECS`, a module in the
metabolic-modelling scope. v3 does not, and cannot: reaching outside the tree
would break the self-containment that is the point of this version. v3's own
observation tables carry `add_mnxr` / `del_mnxr` per row, so the derivation is
mechanical and local.

ONE ROW PER (OBSERVATION, ELEMENT)
----------------------------------
The panel is element-resolved, so an expectation only means something on one
element. A perturbation whose targets carry both C and N is genuinely two
claims and is emitted as two conditions, `<obs_id>__C` and `<obs_id>__N`. v1
could collapse this because each curated spec declared a single element; v3's
observations do not, and inventing one would silently drop the other claim.

THE LOF ARM IS UNIFORM
----------------------
All 166 Keio conditions carry `dir=-`, isozyme controls included. Y states the
CONDUCTANCE claim -- a knockout removes a route, so conductance to the affected
product falls. It states no growth claim.

`essential_on_glucose_minimal` and `is_neg` ride along as STRATIFIER columns.
They select rows for the essentiality diagnostic; they set no direction and gate
nothing. Whether conductance tracks essentiality is a measured question, not an
assumption baked into the answer key.

Usage:  mamba run -n ml python 21_build_Y_v3.py
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _v3bench import (  # noqa: E402
    CURRENCY_NAMES, ELEMENTS, GT_DIR, OBS_DIR, X_DIR, Y_DIR, canon, write_tsv,
)

TARGETS = Path(canon.BENCH_V3_DECISIONS) / "target_resolution.tsv"
CITE_KEIO = "Baba et al. Mol Syst Biol 2006; doi:10.1038/msb4100050 (PMC1681482)"

# `measured` -> direction. Anything not listed carries NO direction and NO
# target cell: an unknown outcome is not a null result, and scoring it as one
# would manufacture 57 false negatives.
DIR_OF = {"up": "+", "increase": "+", "increased": "+", "higher": "+",
          "down": "-", "decrease": "-", "decreased": "-", "lower": "-"}


def _split(cell: str) -> list[str]:
    return [x.strip() for x in re.split(r"[;,]", str(cell)) if x.strip()]


def x_tables():
    met = pd.read_csv(X_DIR / "metabolites.tsv", sep="\t", dtype=str).fillna("")
    in_x = set(met["mnxm"])
    names = dict(zip(met["mnxm"], met["name"]))
    currency = {m for m, n in zip(met["mnxm"], met["name"])
                if n.strip().lower() in CURRENCY_NAMES}
    pats = {el: re.compile(el + r"(?![a-z])") for el in ELEMENTS}
    elmap = {m: {el for el, p in pats.items() if p.search(f)}
             for m, f in zip(met["mnxm"], met["formula"]) if f}
    rp = pd.read_csv(X_DIR / "reaction_participants.tsv", sep="\t", dtype=str).fillna("")
    prod_of = rp[rp.role == "product"].groupby("mnxr")["mnxm"].apply(list).to_dict()
    return in_x, names, currency, elmap, prod_of


def conditions() -> pd.DataFrame:
    in_x, names, currency, elmap, prod_of = x_tables()

    tr = pd.read_csv(TARGETS, sep="\t", dtype=str).fillna("")
    resolved: dict[tuple[str, str], set[str]] = {}
    for _, r in tr[tr.mnxm != ""].iterrows():
        resolved.setdefault((r.arm, r.raw_string), set()).add(r.mnxm)

    rows = []

    # ---- GOF -----------------------------------------------------------------
    gof = pd.read_csv(OBS_DIR / "gof_observations.tsv", sep="\t", dtype=str).fillna("")
    n_nodir = n_notarget = 0
    for _, o in gof.iterrows():
        direction = DIR_OF.get(str(o.measured).strip().lower(), "")
        if not direction:
            n_nodir += 1
            continue

        targets: dict[str, str] = {}          # mnxm -> how it was derived
        for mnxr in _split(o.add_mnxr):
            for m in prod_of.get(mnxr, []):
                if m in in_x and m not in currency:
                    targets.setdefault(m, "mechanical")
        for m in resolved.get(("gof", str(o.target).strip()), ()):
            if m in in_x:
                targets.setdefault(m, "resolved")
        if not targets:
            n_notarget += 1
            continue

        per_el: dict[str, list[str]] = {}
        for m in targets:
            for el in elmap.get(m, ()):
                per_el.setdefault(el, []).append(m)

        for el, ms in per_el.items():
            rows.append(dict(
                condition_id=f"{o.obs_id}__{el}", arm="gof", tier="scored",
                ptype=o.perturbation, host=o.host, gem=o.host_gem,
                n_units=o.n_add, is_control="",
                citation=f"{o.dataset_source} doi:{o.doi}" if o.doi else o.dataset_source,
                note=o.note, gene=o.gene_set, element=el,
                target_mnxm=",".join(sorted(ms)),
                target_name=";".join(names.get(m, "") for m in sorted(ms)),
                target_basis=",".join(sorted({targets[m] for m in ms})),
                expected_dir=direction,
                essential_on_glucose_minimal="", is_neg="",
                obs_id=o.obs_id, host_gem_is_proxy=o.host_gem_is_proxy,
            ))

    print(f"  GOF: {len(gof)} observations -> "
          f"{sum(1 for r in rows if r['arm'] == 'gof')} conditions")
    print(f"       {n_nodir} carry no usable direction in `measured` "
          f"(reported as coverage, never scored as null)")
    print(f"       {n_notarget} have a direction but no target that lands in X")

    # ---- LOF -----------------------------------------------------------------
    lof = pd.read_csv(OBS_DIR / "lof_observations.tsv", sep="\t", dtype=str).fillna("")
    pheno = pd.read_csv(GT_DIR / "pheno_edges.tsv", sep="\t", dtype=str).fillna("")
    lof_gt = pd.read_csv(GT_DIR / "lof_ground_truth.tsv", sep="\t", dtype=str).fillna("")
    ess = dict(zip(lof_gt.gene_name, lof_gt.essential_on_glucose_minimal))
    prod = {r.gene: (r.element, r.dst_mnxm, r["product"], r.is_neg)
            for _, r in pheno.iterrows()}

    n_lof_nt = 0
    for _, o in lof.iterrows():
        # v3 writes the perturbation INTO the gene token -- `argA:del` -- while
        # pheno_edges.tsv and lof_ground_truth.tsv key on the bare gene name.
        # Joining on the raw token matches 0 of 166 and every row then falls
        # through to the free-text path, where roughly half resolve by luck.
        # That failure reads as ordinary coverage loss, which is why it is
        # stripped explicitly here rather than left to a fillna.
        gene = str(o.gene_set).strip().split(":", 1)[0]
        el, dst, pname, is_neg = prod.get(gene, ("", "", "", ""))
        if not (el and dst):
            # fall back to the resolved raw source-metabolite string
            ms = [m for m in resolved.get(("lof", str(o.source_metabolite_raw).strip()), ())
                  if m in in_x]
            # Pick the element the resolved targets actually CARRY. Defaulting
            # to carbon would put the condition on a panel element its product
            # has no atoms of, which scores as a failure to respond rather than
            # as a cell that was never meaningful. A product carrying none of
            # C/N/S/P is not expressible on this panel at all.
            els = sorted({e for m in ms for e in elmap.get(m, ())})
            if not (ms and els):
                n_lof_nt += 1
                continue
            el = els[0]
            ms = [m for m in ms if el in elmap.get(m, ())]
            dst, pname = ",".join(sorted(ms)), ""
        rows.append(dict(
            condition_id=f"lof_{gene}", arm="lof", tier="scored",
            ptype="deletion", host=o.host, gem=o.host_gem,
            n_units=o.n_dead, is_control="", citation=CITE_KEIO, note=o.note,
            gene=gene, element=el, target_mnxm=dst, target_name=pname,
            target_basis="pheno_edges",
            expected_dir="-",                     # UNIFORM -- see module docstring
            essential_on_glucose_minimal=ess.get(gene, ""), is_neg=is_neg,
            obs_id=o.obs_id, host_gem_is_proxy="",
        ))
    print(f"  LOF: {len(lof)} observations -> "
          f"{sum(1 for r in rows if r['arm'] == 'lof')} conditions "
          f"({n_lof_nt} with no target product)")

    return pd.DataFrame(rows)


def expectations(cond: pd.DataFrame, panel: pd.DataFrame) -> pd.DataFrame:
    """Sparse (condition, edge) -> role x dir. Only non-default rows are stored.

    THE DECLARED DEFAULT: any (condition, edge) pair absent from this table is
    role=off_target, dir=0. The scorer materialises the full cross so the
    two-field semantics holds exactly rather than by convention.

    A condition's TARGET cells are every (anchor, its_product) edge on its own
    element -- the panel is anchor x product, so a condition targets a product
    across the whole anchor gradient. That gradient IS the reading: response is
    inversely ordered with anchor degree, so one condition yields one curve
    rather than one number.
    """
    by_el = {el: panel[panel.element == el] for el in ELEMENTS}
    rows = []
    for _, c in cond.iterrows():
        if not c.element or not c.target_mnxm or not c.expected_dir:
            continue
        sub = by_el.get(c.element)
        if sub is None:
            continue
        wanted = [t for t in str(c.target_mnxm).split(",") if t]
        for _, e in sub[sub.product_mnxm.isin(wanted)].iterrows():
            rows.append(dict(
                condition_id=c.condition_id, edge_id=e.edge_id,
                role="target", dir=c.expected_dir, tier=c.tier,
                element=c.element, anchor=e.anchor, anchor_tier=e.anchor_tier,
                basis=c.citation,
            ))
    return (pd.DataFrame(rows)
            .drop_duplicates(subset=["condition_id", "edge_id"]))


def main() -> int:
    Y_DIR.mkdir(parents=True, exist_ok=True)
    panel = pd.read_csv(Y_DIR / "panel_edges.tsv", sep="\t")

    print(f"CONDITIONS\n{'=' * 78}")
    cond = conditions()
    write_tsv(cond, Y_DIR / "conditions.tsv")

    print(f"\n  by arm      : {dict(cond.arm.value_counts())}")
    print(f"  by direction: {dict(cond.expected_dir.value_counts())}")
    print(f"  by element  : {dict(cond.element.value_counts())}")

    print(f"\nEXPECTATIONS\n{'=' * 78}")
    exp = expectations(cond, panel)
    write_tsv(exp, Y_DIR / "expectations.tsv")

    n_cross = len(cond) * len(panel)
    print(f"  {len(exp):,} stored of {n_cross:,} full-cross cells "
          f"({100 * len(exp) / n_cross:.3f}% -- the rest are the declared default)")
    covered = exp.condition_id.nunique()
    print(f"  {covered} of {len(cond)} conditions have at least one target cell")
    if covered < len(cond):
        miss = sorted(set(cond.condition_id) - set(exp.condition_id))
        print(f"  {len(miss)} condition(s) have a target the PANEL cannot express "
              f"-- reported, never dropped: {', '.join(miss[:6])}"
              f"{' ...' if len(miss) > 6 else ''}")

    (Y_DIR / "_y_provenance.json").write_text(json.dumps(dict(
        version="v3",
        n_conditions=int(len(cond)),
        n_expectations=int(len(exp)),
        default_cell="role=off_target, dir=0",
        direction_source=("the observations' `measured` column, NOT `expected` "
                          "-- `expected` is populated on only 10 of 382 GOF rows"),
        unknown_direction=("carries no direction and no target cell; reported as "
                           "coverage, never scored as a null result"),
        target_derivation=[
            "mechanical: non-currency products of the reactions the condition inserts, from X",
            "resolved: the free-text target string mapped via target_resolution.tsv",
        ],
        lof_arm="uniform dir=-, isozyme controls included; a conductance claim, not a growth claim",
        stratifiers=["essential_on_glucose_minimal", "is_neg", "host_gem_is_proxy"],
        independence=("no ECSPr output on any input path; enforced by 50_audit_v3.py"),
    ), indent=2) + "\n")

    print(f"\nY frozen -- {len(cond):,} conditions, {len(exp):,} expectation cells")
    return 0


if __name__ == "__main__":
    sys.exit(main())

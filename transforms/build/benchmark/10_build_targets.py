"""Resolve every free-text target string in v3 to a metabolite id.

WHY THIS IS A TABLE AND NOT A DICT IN A BUILDER
-----------------------------------------------
v3's own decision tables (`host_lineage.tsv`, `gene_reaction_overrides.tsv`)
established the idiom: a judgement call is committed as DATA with the evidence
that decided it, not buried as a literal in code. This is the third such table.
Each row carries `method` (how it resolved) and `reason` (why that is the right
answer), so the table is reviewable without reading this script.

WHAT IT RESOLVES
----------------
267 distinct strings: 184 gain-of-function `target` values plus 90
loss-of-function `source_metabolite_raw` values, less 7 that appear in both.

The LOF arm is read from the RAW column, not the normalized one (90 distinct
against 42). 46 of those cells name more than one metabolite -- `L-isoleucine +
L-valine`, `FMN/FAD` -- and the normalized column collapses exactly the
multi-target information the sparse expectations format is built to express as
several target rows.

RESOLUTION ORDER, strongest first. Every hit must land in `X/metabolites.tsv`:
the benchmark's own sufficiency gate requires every target id to exist in X, so
a resolution to a metabolite outside X is not a resolution, it is a deferred
failure.

  1. exact name in X                 -- unambiguous
  2. normalized name in X            -- punctuation/case/greek-letter fold only
  3. normalized name in MetaNetX chem_prop, then required to be in X
  4. a MetaNetX cross-reference synonym, then required to be in X
  5. curated alias (below), for names no database spells our way
  6. UNRESOLVED -- reported as a coverage miss, never silently dropped

Ambiguity is a REFUSAL, not a coin-flip: if a normalized name matches more than
one metabolite in X, the row resolves to nothing and says so. Picking the first
would be a wrong target that scores as though it were right.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _v3bench import (  # noqa: E402
    DEC_DIR, OBS_DIR, ROOT, X_DIR, canon, norm_name, split_components, write_tsv,
)

# Written into the REPO beside this builder, not into the library. The library
# is built from declarations by the staging script; writing into it directly
# would put a file in the tree that no record describes and no hash covers.
# It is declared as a library item, so it lands in decisions/ on the next build.
OUT = Path(__file__).resolve().parent / "target_resolution.tsv"

# ---------------------------------------------------------------------------
# curated aliases -- names that no MetaNetX synonym spells our way.
# Each entry is a claim that needs the reason beside it, so they live here
# rather than being quietly folded into norm_name().
# ---------------------------------------------------------------------------
CURATED: dict[str, tuple[str, str]] = {
    # normalized string            (mnxm,        reason)
    "biomass": ("", "growth itself, not a metabolite; carries no target cell"),
    "growth": ("", "phenotype, not a metabolite; carries no target cell"),
    "growthrate": ("", "phenotype, not a metabolite; carries no target cell"),
    "none": ("", "explicitly no target -- the non-required isozyme controls"),
    "na": ("", "explicitly no target"),
}


def load_x_index() -> tuple[dict, dict, set]:
    met = pd.read_csv(X_DIR / "metabolites.tsv", sep="\t", dtype=str)
    in_x = set(met["mnxm"])
    exact, normed = {}, {}
    for mnxm, name in zip(met["mnxm"], met["name"].fillna("")):
        if not name:
            continue
        exact.setdefault(name, set()).add(mnxm)
        normed.setdefault(norm_name(name), set()).add(mnxm)
    return exact, normed, in_x


def load_mnx_names(in_x: set) -> dict:
    """normalized name -> {mnxm}, restricted to metabolites present in X."""
    out: dict[str, set] = {}
    cp = pd.read_csv(canon.CHEM_PROP, sep="\t", comment="#", header=None,
                     usecols=[0, 1], names=["mnxm", "name"], dtype=str)
    for mnxm, name in zip(cp["mnxm"], cp["name"].fillna("")):
        if name and mnxm in in_x:
            out.setdefault(norm_name(name), set()).add(mnxm)
    return out


def load_xref_names(in_x: set) -> dict:
    out: dict[str, set] = {}
    cx = pd.read_csv(canon.CHEM_XREF, sep="\t", comment="#", header=None,
                     usecols=[1, 2], names=["mnxm", "desc"], dtype=str)
    for mnxm, desc in zip(cx["mnxm"], cx["desc"].fillna("")):
        if desc and mnxm in in_x:
            out.setdefault(norm_name(desc), set()).add(mnxm)
    return out


def resolve_one(raw: str, exact, normed, mnx, xref) -> tuple[str, str, str]:
    """-> (mnxm, method, reason). Empty mnxm means unresolved or no-target."""
    n = norm_name(raw)
    if n in CURATED:
        mnxm, reason = CURATED[n]
        return mnxm, "curated", reason
    if raw in exact and len(exact[raw]) == 1:
        return next(iter(exact[raw])), "exact_name_in_X", "verbatim match on the X metabolite name"
    for table, method, what in (
        (normed, "normalized_name_in_X", "case/punctuation fold of the X metabolite name"),
        (mnx, "metanetx_name", "MetaNetX chem_prop name, and the id is present in X"),
        (xref, "metanetx_xref", "MetaNetX cross-reference synonym, and the id is present in X"),
    ):
        hits = table.get(n)
        if not hits:
            continue
        if len(hits) > 1:
            return "", "ambiguous", (
                f"{what} matched {len(hits)} distinct metabolites "
                f"({', '.join(sorted(hits)[:4])}...); refused rather than guessed"
            )
        return next(iter(hits)), method, what
    return "", "unresolved", "no exact, folded, MetaNetX or cross-reference name matched an X metabolite"


def main() -> int:
    exact, normed, in_x = load_x_index()
    print(f"X carries {len(in_x)} metabolites, {len(normed)} distinct folded names")
    mnx = load_mnx_names(in_x)
    xref = load_xref_names(in_x)
    print(f"MetaNetX contributes {len(mnx)} folded names and {len(xref)} synonyms that land in X")

    gof = pd.read_csv(OBS_DIR / "gof_observations.tsv", sep="\t", dtype=str)
    lof = pd.read_csv(OBS_DIR / "lof_observations.tsv", sep="\t", dtype=str)

    strings: dict[tuple[str, str], set] = {}
    for arm, col, df in (("gof", "target", gof), ("lof", "source_metabolite_raw", lof)):
        for v in df[col].dropna():
            v = v.strip()
            if v:
                strings.setdefault((arm, v), set())

    rows = []
    for (arm, raw) in sorted(strings):
        comps = split_components(raw)
        if not comps:
            rows.append(dict(arm=arm, raw_string=raw, component=raw, mnxm="",
                             mnxm_name="", method="no_target",
                             reason="cell names no metabolite"))
            continue
        for c in comps:
            mnxm, method, reason = resolve_one(c, exact, normed, mnx, xref)
            name = ""
            if mnxm:
                hit = [k for k, v in exact.items() if mnxm in v]
                name = hit[0] if hit else ""
            rows.append(dict(arm=arm, raw_string=raw, component=c, mnxm=mnxm,
                             mnxm_name=name, method=method, reason=reason,
                             is_compound_cell=len(comps) > 1))

    out = pd.DataFrame(rows)
    write_tsv(out, OUT)

    n_raw = len(strings)
    resolved = out[out.mnxm != ""]
    by_method = out.method.value_counts().to_dict()
    covered = out[out.mnxm != ""].groupby(["arm", "raw_string"]).size().reset_index()
    miss = out[(out.mnxm == "") & (~out.method.isin({"curated", "no_target"}))]

    print(f"\n{n_raw} distinct raw strings -> {len(out)} components")
    print(f"resolved to a metabolite: {len(resolved)} components "
          f"across {len(covered)} of {n_raw} raw strings")
    print("\nby method:")
    for m, c in sorted(by_method.items(), key=lambda kv: -kv[1]):
        print(f"  {m:24s} {c:4d}")
    if len(miss):
        print(f"\n{len(miss)} component(s) UNRESOLVED -- reported, never dropped:")
        for _, r in miss.head(40).iterrows():
            print(f"  [{r.arm}] {r.component!r}  ({r.method})")
        if len(miss) > 40:
            print(f"  ... and {len(miss)-40} more; see {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

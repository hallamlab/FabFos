"""Assemble tier 4 = tier 3 + the MetaCyc curated-AAM increment, STRICTLY ADDITIVE.

Tier 3 (`atom_pairs_working_tail.parquet`) is the base and is not modified. Only
MetaCyc pairs for reactions tier 3 does not already carry are added.

The 15,539 reactions where MetaCyc could REPLACE a predicted map with a curated one
are deliberately excluded. They agree with tier 3 at 90.8% exact atom-to-atom and
98.6% molecule-level, so the disagreement is small -- but replacing is not adding,
it changes reactions that current results depend on, and it would need its own
before/after benchmark rather than riding along inside a coverage change. That path
is recorded in the scope journal for a decision.

Gate before writing: an added reaction must not already exist in tier 3 for that
element, and every added atom pair must address an atom via the same canonical-rank
system tier 3 uses -- otherwise the two halves of the graph address different atoms.
The rank system is enforced upstream in `metacyc_pairs.canonical_ranks`, which is a
verbatim copy of the producer's; here we re-assert the additive property and refuse
on any overlap rather than silently letting one side win.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

REF = Path("/home/tony/agentic_workspace/data/scadc/ecspr_reference/mnxref-4_5")
TMP = Path(__file__).parent
OUT = TMP / "atom_pairs_tier4.parquet"

ELEMENTS = ["C", "N", "S", "P"]


def main() -> int:
    t3 = pd.read_parquet(REF / "recovery/atom_pairs_working_tail.parquet")
    mc = pd.read_parquet(TMP / "metacyc_atom_pairs.parquet")
    print(f"tier 3          : {len(t3):,} pairs, {t3.mnxr.nunique():,} reactions")
    print(f"metacyc curated : {len(mc):,} pairs, {mc.mnxr.nunique():,} reactions")

    # additive only: per (mnxr, element), tier 3 wins if it has the reaction at all
    have = set(zip(t3.mnxr, t3.element))
    add = mc[~pd.Series(list(zip(mc.mnxr, mc.element)), index=mc.index).isin(have)]
    print(f"\nadditive slice  : {len(add):,} pairs, {add.mnxr.nunique():,} reactions")
    dropped = mc.mnxr.nunique() - add.mnxr.nunique()
    print(f"withheld (tier 3 already has them; upgrade path deferred): {dropped:,} reactions")

    t4 = pd.concat([t3, add], ignore_index=True)

    # --- gates -----------------------------------------------------------------
    ok = True
    # Duplicate pairs are part of the EXISTING convention, not a defect: a metabolite
    # with stoichiometry >1 contributes the same (metabolite, rank) twice, and tier 3
    # carries 3,283 such rows (0.13%) on its own. So the gate is not "no duplicates" --
    # it is "the merge introduces no CROSS-SOURCE collision", i.e. the additive slice
    # must not address an atom pair tier 3 already addresses. That is the thing that
    # would double-count weight; a source's own internal repeats already balance.
    key = ["mnxr", "element", "substrate", "product", "sub_idx", "prod_idx"]
    collide = len(set(map(tuple, add[key].values)) & set(map(tuple, t3[key].values)))
    if collide:
        print(f"GATE FAIL: {collide:,} atom pairs collide across sources")
        ok = False
    print(f"  gate: cross-source collisions {collide}  "
          f"(tier3 internal dups {t3.duplicated(subset=key).sum():,}, "
          f"added internal dups {add.duplicated(subset=key).sum():,})")
    for e in ELEMENTS:
        a, b = t3[t3.element == e].mnxr.nunique(), t4[t4.element == e].mnxr.nunique()
        if b < a:
            print(f"GATE FAIL: {e} lost reactions ({a:,} -> {b:,})")
            ok = False
    if (t4.sub_idx < 0).any() or (t4.prod_idx < 0).any():
        print("GATE FAIL: negative atom index")
        ok = False
    if not ok:
        print("\nREFUSED -- not writing tier 4")
        return 1

    t4.to_parquet(OUT, index=False)
    print(f"\n=== tier 4 written: {OUT} ===")
    print(f"{'el':<4}{'tier3':>10}{'tier4':>10}{'delta':>9}{'pairs t4':>12}")
    for e in ELEMENTS:
        a = t3[t3.element == e].mnxr.nunique()
        b = t4[t4.element == e].mnxr.nunique()
        n = (t4.element == e).sum()
        print(f"{e:<4}{a:>10,}{b:>10,}{b-a:>+9,}{n:>12,}")
    print(f"\ntotal pairs {len(t3):,} -> {len(t4):,}   "
          f"reactions {t3.mnxr.nunique():,} -> {t4.mnxr.nunique():,}")
    print(f"curated-source pairs in tier 4: {(t4.source=='metacyc_aam').sum():,}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""Structure audit of target_resolution.tsv -- confidence, not a re-resolution.

10_build_targets.py resolves a target string to a metabolite by NAME, through
four tiers (exact in X, folded in X, chem_prop name, chem_xref synonym). It
returns on the FIRST tier that hits, and reports ambiguity only WITHIN that
tier. So a string that matches exactly one metabolite by X-name and three more
by xref synonym resolves silently, with nothing recording that the other three
existed. That is the error mode this audit is for: not the 5 refusals, which
already abstain, but the rows that look clean because nobody looked further
down.

The independent axis is structure. For every row we collect the candidate
metabolites across ALL four tiers, then collapse them by InChIKey connectivity
layer (the first block, before the first hyphen). Two metabolites sharing a
skeleton are the same molecule differing only in protonation, charge or
stereochemistry -- for an atom-transfer network those are the same node for our
purposes, so a "multi-match" that collapses to one skeleton was never really
ambiguous. Two metabolites with different skeletons are different molecules,
and picking one by tier precedence is a coin-flip the resolver did not disclose.

Verdicts:
  UNIQUE     one candidate, or all candidates share one InChIKey skeleton.
             The name cannot be pointing somewhere else.
  COLLISION  candidates span >1 skeleton -- genuinely different molecules.
             Whatever the resolver picked, it picked by tier order, not by
             evidence. These are the review queue.
  NOSTRUCT   at least one candidate carries no InChIKey, so the skeleton test
             cannot run and the row is unverifiable by this method.

What this does NOT prove: that the resolved metabolite is the one the benchmark
author meant. Every tier here is still MetaNetX name/synonym matching, and
MetaNetX is the only structure source on disk -- so this bounds the internal
ambiguity of the mapping, it does not validate it against an outside authority.
A UNIQUE verdict means "no other metabolite answers to this name", not "this is
the right metabolite".

Table-only. Writes structure_audit.tsv beside target_resolution.tsv; changes
nothing the scorer reads.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from _v3bench import X_DIR, canon, norm_name, write_tsv  # noqa: E402

HERE = Path(__file__).parent
SKEL_UNKNOWN = ""


def inchikey_layers() -> tuple[dict[str, str], dict[str, str]]:
    """mnxm -> (connectivity layer, full InChIKey); '' when absent.

    The two layers fail in OPPOSITE directions, which is why both are kept:

      * The first block is STEREO-BLIND -- stereo lives in the second block.
        (S)-lactate and (R)-lactate share a skeleton. So a skeleton match does
        NOT rule out an enantiomer swap, which for a metabolic network is a
        real error (different enzymes act on D- and L-lactate).
      * The first block is PROTONATION-SENSITIVE -- H counts are part of the
        connectivity layer. FMN at charge 0 and FMN at charge -3 get different
        skeletons despite being the same compound. So a skeleton mismatch does
        NOT prove two different molecules.

    Neither layer alone is the right test. We report both and let the verdict
    name which risk applies.
    """
    cp = pd.read_csv(canon.CHEM_PROP, sep="\t", comment="#", header=None,
                     usecols=[0, 7], names=["mnxm", "ik"], dtype=str)
    skel, full = {}, {}
    for mnxm, ik in zip(cp["mnxm"], cp["ik"].fillna("")):
        skel[mnxm] = ik.split("-")[0] if ik else SKEL_UNKNOWN
        full[mnxm] = ik
    return skel, full


def all_tier_candidates() -> dict[str, set[str]]:
    """normalized name -> every X metabolite any tier would offer for it.

    Deliberately a UNION, where the resolver takes the first hit. The union is
    what tells us whether the resolver had a choice it did not report.
    """
    met = pd.read_csv(X_DIR / "metabolites.tsv", sep="\t", dtype=str)
    in_x = set(met["mnxm"])
    cand: dict[str, set[str]] = {}

    def add(name: str, mnxm: str) -> None:
        if name and mnxm in in_x:
            cand.setdefault(norm_name(name), set()).add(mnxm)

    for mnxm, name in zip(met["mnxm"], met["name"].fillna("")):
        add(name, mnxm)
    cp = pd.read_csv(canon.CHEM_PROP, sep="\t", comment="#", header=None,
                     usecols=[0, 1], names=["mnxm", "name"], dtype=str)
    for mnxm, name in zip(cp["mnxm"], cp["name"].fillna("")):
        add(name, mnxm)
    cx = pd.read_csv(canon.CHEM_XREF, sep="\t", comment="#", header=None,
                     usecols=[1, 2], names=["mnxm", "desc"], dtype=str)
    for mnxm, desc in zip(cx["mnxm"], cx["desc"].fillna("")):
        add(desc, mnxm)
    return cand


def main() -> int:
    tr = pd.read_csv(HERE / "target_resolution.tsv", sep="\t", dtype=str).fillna("")
    skel, full_ik = inchikey_layers()
    cand = all_tier_candidates()

    rows = []
    for r in tr.itertuples(index=False):
        comp = r.component
        hits = sorted(cand.get(norm_name(comp), set()))
        skels = {skel.get(m, SKEL_UNKNOWN) for m in hits}
        has_unknown = SKEL_UNKNOWN in skels
        real = skels - {SKEL_UNKNOWN}

        if not hits:
            verdict, note = "NOCAND", "no tier offers a candidate for this string"
        elif has_unknown:
            verdict = "NOSTRUCT"
            note = (f"{len(hits)} candidate(s), at least one with no InChIKey -- "
                    "skeleton test cannot run")
        elif len(real) == 1:
            # One skeleton. Now ask the stereo question the skeleton cannot:
            # do the candidates differ in the stereo layer of the InChIKey?
            stereo = {full_ik.get(m, "") for m in hits}
            if len(stereo) > 1:
                verdict = "STEREO_SPLIT"
                note = (f"{len(hits)} candidate(s) share skeleton {sorted(real)[0]} "
                        f"but differ in the InChIKey stereo layer -- same "
                        "connectivity, possibly different stereoisomer")
            else:
                verdict = "UNIQUE"
                note = (f"{len(hits)} candidate(s), one InChIKey "
                        f"({sorted(real)[0]}) -- no other metabolite answers "
                        "to this name")
        else:
            verdict = "COLLISION"
            note = (f"{len(hits)} candidate(s) span {len(real)} distinct InChIKey "
                    f"skeletons -- different connectivity OR only different "
                    "protonation; formula comparison decides which")

        # Did the resolver's pick even survive the union? A pick outside its own
        # candidate set would mean the two paths disagree about the name.
        picked = r.mnxm
        in_cand = "" if not picked else ("yes" if picked in hits else "NO")

        rows.append({
            "arm": r.arm,
            "component": comp,
            "method": r.method,
            "picked_mnxm": picked,
            "picked_in_candidates": in_cand,
            "n_candidates": len(hits),
            "n_skeletons": len(real),
            "verdict": verdict,
            "candidates": ";".join(hits[:8]),
            "note": note,
        })

    out = pd.DataFrame(rows)
    write_tsv(out, HERE / "structure_audit.tsv")

    resolved = out[out["picked_mnxm"] != ""]
    print(f"{len(out)} rows audited, {len(resolved)} carry a picked metabolite\n")
    print("-- verdict over ALL rows")
    print(out["verdict"].value_counts().to_string())
    print("\n-- verdict over rows the resolver RESOLVED (the ones that move AUCs)")
    print(resolved["verdict"].value_counts().to_string())
    print("\n-- resolved rows by method x verdict")
    print(pd.crosstab(resolved["method"], resolved["verdict"]).to_string())

    bad = resolved[resolved["picked_in_candidates"] == "NO"]
    print(f"\n-- picks outside the union of candidates: {len(bad)}")
    if len(bad):
        print(bad[["component", "picked_mnxm", "method"]].to_string(index=False))

    refused = out[out["method"] == "ambiguous"]
    print(f"\n-- the {len(refused)} refusals, re-examined by structure")
    for r in refused.itertuples(index=False):
        print(f"  {r.component:<24} {r.verdict:<10} {r.note}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

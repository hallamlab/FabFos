"""Turn MetaCyc's CURATED atom maps into ECSPr atom pairs.

`atom-mappings-smiles.dat` holds 16,890 expert-assigned reaction atom maps -- every
atom carries a map number, and per-element mass balance holds at 99.8-100%. It has
never been read by anything in this tree. The AAM ensemble's "MetaCyc member" was
structures, not these maps.

The hard part is not the chemistry, it is the ADDRESSING. Tier 3's atom pairs are
keyed by (mnxr, element, substrate MNXM, product MNXM, sub_idx, prod_idx), where the
index is the atom's position in *that MNXM's own SMILES parse*. MetaCyc gives mapped
SMILES fragments with no MNXM labels and its own atom ordering. So each fragment has
to be identified as an MNXM and then re-indexed into that MNXM's frame. A silent
error here does not raise -- it emits a plausible pair between the wrong atoms.

Hence the validation this script exists to support: 15,539 of these reactions are
ALREADY in tier 3 with consensus pairs derived independently by MCS+RXNMapper. If the
addressing is right, the two must largely agree on those. Disagreement is the signal
that the alignment is wrong, and it is checked before anything is banked.

Substitutions applied to make RDKit accept the source, and their justification:
  [R:n] / [R1:n]      -> [*:n]   R-group: a real atom of unspecified identity.
  [a protein:n] etc.  -> [*:n]   free-text carrier names. These are exactly the
                                 placeholder class: a dummy contributes no C/N/S/P,
                                 which is the same semantics `ecspr_aam_rescue` gives
                                 a placeholder, and the balance gate still judges the
                                 concrete atoms.
  [L-cysteine:n]      -> LEFT ALONE, deliberately. Cysteine carries real C, N and S;
                                 dummying it would silently delete sulfur atoms. The
                                 954 records that still fail are refused, not coerced.
"""
from __future__ import annotations

import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

D = Path("/home/tony/agentic_workspace/data/scadc/references")
MX = D / "metanetx"
REF = Path("/home/tony/agentic_workspace/data/scadc/ecspr_reference/mnxref-4_5")
OUT = Path(__file__).parent

ELEMENTS = ["C", "N", "S", "P"]
TERM = re.compile(r"(\d+)\s+(MNXM\d+|WATER|BIOMASS)")
RGRP = re.compile(r"\[R(\d*):(\d+)\]")
TEXT = re.compile(r"\[(?:a|an|the)\s+[^\]:]+:(\d+)\]")


def fix(s: str) -> str:
    return TEXT.sub(r"[*:\1]", RGRP.sub(r"[*:\2]", s))


def canon(smiles: str) -> str | None:
    """Canonical SMILES with map numbers and charges stripped, for identity matching."""
    m = Chem.MolFromSmiles(smiles, sanitize=False)
    if m is None:
        return None
    for a in m.GetAtoms():
        a.SetAtomMapNum(0)
    try:
        return Chem.MolToSmiles(m)
    except Exception:
        return None


def load_metanetx():
    smi, eq = {}, {}
    for line in open(MX / "chem_prop.tsv"):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) > 8 and f[8]:
            smi[f[0]] = f[8]
    for line in open(MX / "reac_prop.tsv"):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) > 1 and f[1]:
            eq[f[0]] = f[1]
    mc2mnx = {}
    for line in open(MX / "reac_xref.tsv"):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 2 or ":" not in f[0]:
            continue
        db, val = f[0].split(":", 1)
        if db.startswith("metacyc"):
            mc2mnx.setdefault(val, f[1])
    return smi, eq, mc2mnx


def side_participants(side: str) -> list[str]:
    out = []
    for n, m in TERM.findall(side):
        out.extend([m] * int(n))
    return out


def frag_to_mnxm(frag_mol, cand_canon: dict[str, list[str]]):
    """Identify a mapped fragment as one of the reaction's declared participants."""
    c = canon(Chem.MolToSmiles(frag_mol))
    if c is None:
        return None
    hits = cand_canon.get(c)
    return hits[0] if hits else None


def canonical_ranks(mol):
    """Atom -> canonical rank. COPIED VERBATIM from `ecspr_atom_pairs.canonical_ranks`.

    Not reimplemented, and not replaced with a substructure match against the MNXM's
    own parse -- that was the first attempt here and it agreed with tier 3 on only
    22.6% of atom pairs while agreeing on 98.6% of molecule-level transfers, which is
    the exact signature the original docstring describes: the molecules were right and
    the atom addressing was wrong.

    `(metabolite, GetIdx())` is not an atom -- the same metabolite is written
    differently in different reactions, so a graph keyed on it welds unrelated atoms
    together. `CanonicalRankAtoms` is invariant to input ordering, so `(met, rank)` IS
    an atom. Sanitize unconditionally: ranks depend on perceived aromaticity, so
    ranking some molecules sanitized and others raw gives one metabolite two rank
    systems and reintroduces the bug.

    This must stay byte-identical in behaviour to the producer, or MetaCyc-derived
    pairs and tier-3 pairs address different atoms and cannot share a graph.
    """
    m2 = Chem.Mol(mol)
    for a in m2.GetAtoms():
        a.SetAtomMapNum(0)
    try:
        Chem.SanitizeMol(m2)
        return list(Chem.CanonicalRankAtoms(m2, breakTies=True))
    except Exception:
        return None


def build_side(side_smiles: str, participants: list[str], smi: dict):
    """-> list of (mnxm, {map_num: mnxm_atom_idx}, Counter(element counts))"""
    mol = Chem.MolFromSmiles(side_smiles, sanitize=False)
    if mol is None:
        return None
    cand_canon = defaultdict(list)
    for p in participants:
        s = smi.get(p)
        if s:
            c = canon(s)
            if c:
                cand_canon[c].append(p)
    used = Counter()
    out = []
    for frag_idx in Chem.GetMolFrags(mol):
        sub = Chem.PathToSubmol(mol, []) if False else None
        frag = Chem.RWMol(mol)
        keep = set(frag_idx)
        for i in sorted(set(range(mol.GetNumAtoms())) - keep, reverse=True):
            frag.RemoveAtom(i)
        frag = frag.GetMol()
        order = sorted(frag_idx)          # RemoveAtom preserves relative order
        c = canon(Chem.MolToSmiles(frag))
        if c is None:
            continue
        pool = [p for p in cand_canon.get(c, []) if used[p] < participants.count(p)]
        if not pool:
            continue
        mnxm = pool[0]
        used[mnxm] += 1
        ranks = canonical_ranks(frag)
        if ranks is None:                 # genuinely unrankable -> refuse, never approximate
            continue
        m2i = {}
        for local, global_i in enumerate(order):
            a = mol.GetAtomWithIdx(global_i)
            if a.GetAtomMapNum() > 0:
                m2i[a.GetAtomMapNum()] = (ranks[local], a.GetSymbol())
        out.append((mnxm, m2i))
    return out


def main() -> int:
    limit = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    smi, eq, mc2mnx = load_metanetx()

    maps = {}
    for line in open(D / "metacyc26_flatfiles/atom-mappings-smiles.dat", encoding="latin-1"):
        if "\t" not in line:
            continue
        k, v = line.rstrip("\n").split("\t", 1)
        if ">>" in v:
            maps[k] = v
    print(f"[mc] curated maps: {len(maps):,}", flush=True)

    rows, stats = [], Counter()
    items = [(k, v) for k, v in maps.items() if k in mc2mnx]
    if limit:
        items = items[:limit]
    print(f"[mc] with an MNXR: {len(items):,}", flush=True)

    for i, (mcid, rxn) in enumerate(items, 1):
        if i % 1000 == 0:
            print(f"[mc]   {i:,}/{len(items):,} banked={len(rows):,}", flush=True)
        mnxr = mc2mnx[mcid]
        e = eq.get(mnxr)
        if not e or "=" not in e:
            stats["skip:no_equation"] += 1
            continue
        lhs_p, rhs_p = (side_participants(x) for x in e.split("=", 1))
        f = fix(rxn)
        if ">>" not in f:
            stats["skip:no_arrow"] += 1
            continue
        ls, rs = f.split(">>", 1)
        L = build_side(ls, lhs_p, smi)
        R = build_side(rs, rhs_p, smi)
        if not L or not R:
            stats["skip:unaligned"] += 1
            continue
        # map number -> (mnxm, idx, symbol) on each side
        lmap, rmap = {}, {}
        for mnxm, m2i in L:
            for mn, (idx, sym) in m2i.items():
                lmap.setdefault(mn, (mnxm, idx, sym))
        for mnxm, m2i in R:
            for mn, (idx, sym) in m2i.items():
                rmap.setdefault(mn, (mnxm, idx, sym))
        n_before = len(rows)
        for mn, (sm, si, sym) in lmap.items():
            if sym not in ELEMENTS or mn not in rmap:
                continue
            pm, pi, psym = rmap[mn]
            if psym != sym:
                stats["drop:element_changed"] += 1
                continue
            rows.append((mnxr, sym, sm, pm, si, pi, 1.0, "curated", "metacyc_aam", 1.0))
        if len(rows) > n_before:
            stats["banked"] += 1
        else:
            stats["skip:no_pairs"] += 1

    out = pd.DataFrame(rows, columns=["mnxr", "element", "substrate", "product",
                                      "sub_idx", "prod_idx", "pair_w", "method",
                                      "source", "confidence"])
    out.to_parquet(OUT / "metacyc_atom_pairs.parquet", index=False)
    print(f"\n=== {len(out):,} atom pairs from {stats['banked']:,} reactions ===")
    for k, v in sorted(stats.items()):
        print(f"  {k:<28} {v:>7,}")
    if len(out):
        print("\n  per element:")
        for e in ELEMENTS:
            s = out[out.element == e]
            print(f"    {e}  pairs={len(s):>8,}  reactions={s.mnxr.nunique():>6,}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

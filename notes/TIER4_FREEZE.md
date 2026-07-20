# Tier 4 atom-pair universe — frozen 2026-07-20

    atom_pairs_tier4.parquet
      sha256 7b3b217f91373f3141404767f9be3ad20a87be10d756ccef6a7ecdb1311abe93
      2,455,235 pairs / 63,621 reactions / 15,963,618 bytes

    recovery/metacyc_atom_pairs.parquet          (the increment, unmerged form)
      sha256 97ce397535434851054d92acf1533e5edaf127929e2636d41f28046d03670cbf
      105,891 pairs / 5,246 reactions

## What it is

Tier 3 (`recovery/atom_pairs_working_tail.parquet`), unmodified, plus a **strictly
additive** increment of curated atom maps from MetaCyc.

| element | tier 3 | tier 4 | delta |
|---|---|---|---|
| C | 61,766 | 61,989 | +223 |
| N | 47,400 | 47,532 | +132 |
| S | 15,570 | 15,624 | +54 |
| P | 34,846 | 34,922 | +76 |

Added rows carry `method=curated`, `source=metacyc_aam`, `confidence=1.0`, and are
the only rows in the file not derived from prediction.

## Where the increment came from

`references/metacyc26_flatfiles/atom-mappings-smiles.dat` — 16,890 expert-assigned
reaction atom maps that **had never been read by anything in this tree** (grep for
the filename across the workspace returns no consumer). The AAM ensemble's "MetaCyc
member" was structures, not these maps.

Quality of the source, measured: every atom carries a map number (100%), and
per-element mass balance holds at 99.8% C / 99.9% N / 99.9% S / 100.0% P. For
comparison, the RXNMapper-derived T1 increment attempted the same day balanced
carbon in 12.5% of candidates.

Reach: 16,845 of 16,890 dereference to an MNXR via `reac_xref`.

## The validation that licenses this

15,539 of these reactions already exist in tier 3 with pairs derived independently
by MCS+RXNMapper. Comparing the two on the shared set:

- **molecule-level agreement 98.6%** — which molecule donates atoms to which
- **exact atom-to-atom agreement 90.8%**

Two independent methods converging at that rate is what justifies trusting the
1,306 reactions where only MetaCyc has an answer.

**The 90.8% was 22.6% on the first attempt.** The cause was atom addressing, not
chemistry: the first pass aligned fragments to MNXM structures with
`GetSubstructMatch`, yielding `GetIdx()`-based indices. `ecspr_atom_pairs.py`
documents at length why that is wrong — the same metabolite is written differently in
different reactions, so `(metabolite, GetIdx())` is a different atom depending on
which reaction you read it from. Switching to `CanonicalRankAtoms(breakTies=True)` on
a sanitized molecule — a verbatim copy of the producer's `canonical_ranks` — moved
agreement to 90.8%. **Any future producer feeding this file must use that same rank
system, or the two halves of the graph address different atoms.**

## What was deliberately NOT done

- **The 15,539-reaction upgrade path is excluded.** MetaCyc could replace predicted
  maps with curated ones for 25% of the carbon universe. That is a replacement, not
  an addition; it changes reactions current results depend on and needs its own
  before/after benchmark rather than riding inside a coverage change. Recorded in the
  scope journal for a decision.
- **954 of 16,890 source records are refused, not coerced.** After substituting
  `[R:n]` and free-text carrier names (`[a protein:n]`) with dummy atoms, these still
  fail to parse; they contain `[L-cysteine:n]`-style named residues and Fe-S centres.
  Cysteine carries real C, N and S, so dummying it would silently delete sulfur.
- **11,389 records extracted no pairs** because their fragments could not be
  identified as declared MetaNetX participants — MetaCyc and MetaNetX differ on
  protonation and tautomer. Recoverable in principle by matching on InChIKey
  connectivity layer rather than canonical SMILES; not attempted.

## Gates applied

- additive only: an added `(mnxr, element)` must be absent from tier 3
- **cross-source atom-pair collisions: 0** — the added pairs address no atom tier 3
  already addresses. Note the gate is *not* "no duplicate pairs": tier 3 carries
  3,283 internal duplicates (0.13%) from stoichiometry >1, which is the existing
  convention, and an earlier version of this gate wrongly refused the build over them.
- no negative atom indices; no element may lose reactions

## Not yet done — this freeze does not include

Graph rebuild and canon repoint. `ecspr_atom_graph.py` is a library consumed by the
transform harness, not a CLI, so building `mnx_bipartite_*.pkl` from this parquet and
re-pinning `canon.py` requires that harness. **Until that happens this file is frozen
but unread** — which is precisely the condition tier 3 has been in since 2026-07-17,
and the reason this whole exercise started. It should not be left that way.

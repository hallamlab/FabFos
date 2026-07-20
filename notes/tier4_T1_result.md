# T1 executed: the free-bucket premise does not survive a real balance gate

Measured 2026-07-20. Raw output: `tier4_t1_free_increment.tsv` (1,458 rows).

## What was claimed

The plan rated T1 the highest-confidence lever: **868 reactions, "high" confidence,
needs nothing**. Three groups — 691 blank-reason closures, 117 already-balanced
`unmapped_computable`, ~60 `generic_rgroup` — that supposedly only needed a mapping
run, because the closure ledger had refused them without attempting them.

The supporting evidence was a 24/24 sample of blank-reason reactions that "map and
balance".

## What actually happened

1,497 candidates were built from MetaNetX equations and mapped with RXNMapper.
1,458 mapped. Then the gate was applied correctly, and the number collapsed.

**The first run's balance verdict was vacuous.** It tested `net[element] == 0`,
which is trivially true when the reaction contains no atoms of that element at all.
That reported 1,393/1,458 reactions "balancing" sulfur — for reactions that are
mostly sulfur-free. Requiring the element to be *present as well as conserved*:

| element | present | present **and** balanced |
|---|---|---|
| C | 1,458 | **183** |
| N | 979 | 696 |
| S | 228 | 182 |
| P | 532 | 388 |

Carbon is the element that carries the universe (tier 3 = 61,766 carbon reactions).
**183 net-new carbon reactions, not 868.** Confirmed net-new: zero overlap with
tier 3.

By bucket, against what each was predicted to yield:

| bucket | banked | carbon-balanced | predicted |
|---|---|---|---|
| blank_reason | 687 | **20** | ~691 |
| unmapped_computable | 703 | **118** | 117 |
| generic_rgroup | 68 | **45** | ~60 |

## The two conclusions this forces

**1. The "bookkeeping defect" was mostly my own measurement error.** The plan called
the 691 blank-reason closures "the clearest instance of missed vs reverted" and
claimed 24/24 balanced on sample. Under a presence-requiring gate, 20 of 687 balance
carbon — 2.9%. The 24/24 was the vacuous gate counting sulfur-free reactions as
sulfur-balanced. The prior curation's closure of this bucket was substantially
correct. What remains true is narrower: those rows should carry a reason rather than
a blank, which is a bookkeeping fix worth roughly 20 carbon reactions, not 691.

The one bucket that held up is `unmapped_computable` at 118 ≈ the 117 independently
identified as already balancing. That agreement is the reassuring signal here — two
routes to the same number.

**2. Banking these would reproduce the defect the plan set out to fix.** The
carbon-balanced 183 have RXNMapper confidence **mean 0.281, median 0.156, 79% below
0.5**. The plan's central charge against tier 1 is that its strip lane carries mean
confidence 0.536 with 45% below 0.5 and gates on none of it. This increment is
*worse on both figures*. Applying a 0.5 floor leaves **39 reactions**.

Against the standing constraint — "make sure they are genuine improvements" — 183
reactions at median confidence 0.156 is not one.

## What this implies for T2–T6

T1's 868 was derived by classifying names into levers, not by walking them. It came
in 4.7× high. T2 (487, "high") and T3 (497, "medium-high") were estimated the same
way and carry the same unquantified error. The plan's headline — ~2,150 reactions
over tier 3, landing near tier 1's count — should not be treated as load-bearing
until each remaining lever is executed rather than estimated.

This does not invalidate the audit's core finding, which was independently measured
and stands: three universes exist, the curation landed in one that nothing reads,
and tier 1's surplus is 84% an unbalanced-by-construction strip lane whose
confidence gate was specified and never implemented. The repoint from tier 1 to a
gated universe remains the real work. What has changed is the expectation that the
gap can be closed on reaction count — on present evidence it largely cannot, and
tier 4 would win on edge quality while losing on the headline number.

## Also fixed

`ecspr_atom_pairs.py:495` read the structure crosswalk with `comment="#"`, silently
truncating any curated SMILES containing a triple bond (`C#N`, `C#C`) at the bond.
`ecspr_aam_rescue.py` reads the same file and was already fixed; the two modules
could disagree about what a curated row says. Ported (submodule `508cb77`).

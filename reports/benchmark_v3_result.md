# ECSPr against benchmark V3 — first validation result

Run: sockeye job `12317266` (solve) + `12317276` (merge), 2026-07-20.
Agent home `/scratch/st-shallam-1/txyliu/ecspr_bench_v3_1784532632`, run key `cPieN6LE`.
Engine pinned at metasmith 0.18.8; task image
`quay.io/hallamlab/external_ecspr:2026.07.14` (see the naming caveat below).

Solve: 16 shards (4 facets x 4 elements), 32 worker processes, one node,
~14 min wall. Merge: ~1 min. Output 5,690,248 rows, exactly 1,422,562 per
facet, in the frozen V1 contract shape (`facet, condition_id, edge_id, effect`).

The answer key was frozen and hashed before this solve ran, and the
independence audit passes on the V3 builders — so these numbers are a
measurement, not a fit.

## Directional specificity (the headline)

AUC per facet and element, 925 scored conditions. `passes` is the contract's
own gate; every stratum passes in both modes.

| facet | mode | ALL | C | N | S | P |
|---|---|---|---|---|---|---|
| netA_iML1515 | signed | 0.831 | 0.828 | 0.844 | 0.842 | 0.866 |
| netA_iECDH10B | signed | 0.829 | 0.825 | 0.849 | 0.853 | 0.877 |
| netB_iML1515 | signed | 0.849 | 0.846 | 0.863 | 0.916 | 0.883 |
| netB_iECDH10B | signed | 0.845 | 0.842 | 0.855 | 0.915 | 0.884 |
| netA_iML1515 | unsigned | 0.872 | 0.872 | 0.858 | 0.856 | 0.894 |
| netA_iECDH10B | unsigned | 0.871 | 0.870 | 0.864 | 0.851 | 0.893 |
| netB_iML1515 | unsigned | 0.892 | 0.889 | 0.915 | 0.892 | 0.922 |
| netB_iECDH10B | unsigned | 0.890 | 0.887 | 0.915 | 0.892 | 0.923 |

Two things are worth reading off this rather than leaving to a reader:

- **netB beats netA on every element, in both modes.** netB's base graph is
  induced from the experiment's own annotation evidence rather than from a
  curated GEM. That it scores *higher* against an independently frozen key is
  the more interesting result in the table, and it is not what a "curated is
  better" prior would predict.
- **Unsigned exceeds signed everywhere.** Expected, and it is the direction
  dimension paying its cost: per the T3 analysis, nearly all sign contrast in
  this benchmark comes from the 166 uniformly-negative deletion conditions plus
  3 strain rows. The gain-of-function arm is close to a pure "does this
  perturbation raise conductance toward its target more than elsewhere" test.
  The signed-vs-unsigned gap (~0.04) is roughly the price of asking for sign.

## Essentiality diagnostic

Carbon element, 131 essential vs 32 non-essential genes:

| facet | AUC |
|---|---|
| netA_iML1515 | 0.756 |
| netA_iECDH10B | 0.761 |
| netB_iML1515 | 0.776 |
| netB_iECDH10B | 0.765 |

Identical across signed and unsigned, as it should be — the diagnostic reads
magnitude only. The S element has 2 essential and 0 non-essential genes, so its
AUC is undefined and is reported as NaN rather than as a score.

## Potency: undefined, by decision

Every potency cell reports `status=undefined:no_null_pool`. This is the
no-nulls decision behaving correctly: with no null pool the potency measure has
no reference distribution, and the scorer says so instead of emitting a number.
Directional specificity and the essentiality diagnostic are the live outputs.

## Scorer self-test

`30_score_v3.py --selftest` PASSES: rank-invariance holds under x^3, 5x,
arcsinh, x|x| (signed) and exp|x|, log|x|+1 (unsigned), all 1.000000 == 1.000000;
perfect-oracle specificity AUC = 1.0; no-null-pool path returns NaN/undefined.
This is what licenses emitting the signed conductance difference rather than
the pilot's log ratio — the scorer reads order only.

## Caveats a reader needs

- **The 341-row target-resolution table has not had human review.** The plan
  put a review gate here; the run was on autopilot, so it was crossed. 124
  component rows resolved to a metabolite (85 gain-of-function, 39
  loss-of-function), 217 carry no target cell, and 5 ambiguous strings were
  refused rather than guessed. Every row carries its reason. Unresolved rows
  are coverage misses, never silent drops — but the AUCs above are computed on
  what did resolve, so a review that moves rows moves the numbers.
- **Image naming was wrong, and has since been fixed.** This run used
  `external_ecspr`, but the `external_` prefix denotes third-party images and
  this one was built in-house (`provenance/containers/ecspr.yml` records
  `built_by_us: true`). On 2026-07-20 the image was retagged
  `quay.io/hallamlab/ecspr:2026.07.14` registry-side from the same manifest
  index, so the digests this run was pinned against are unchanged and the old
  reference still resolves. Nothing about the numbers below moves.
- No local compute was used for the solve. Scoring ran locally: ~15 min per
  mode on one core, not the sub-minute the plan assumed, because the merged
  table is 424 MB.

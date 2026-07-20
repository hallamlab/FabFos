# provenance/ — declarations

This tree holds **declarations**: what each data dependency is, where it came
from, and what its bytes should hash to. It is in git, reviewable in a diff.

The counterpart is `.awm/data/ref/provenance/`, which holds **receipts**: what
actually happened on a specific machine during a specific build (host, wall
clock, realized RNG state). Git holds the claim; the library holds the event.
**A disagreement between them is the alarm** — that is the whole reason both
exist.

## Layout

```
provenance/
├── data/        one record per library item
│   └── _declared.yml    the single source-of-truth table (hand-edited)
├── containers/  one record per image (generated — see below)
└── env/         conda specs (placeholder; see the provision.py callout)
```

`provenance/build/` deliberately does **not** exist. A `kind: derived` record
already carries its producing transform, that transform's hash, its inputs with
their hashes, its parameters and its determinism class — that *is* build
provenance, attached to the artifact rather than floating in a parallel tree.
The build DAG itself is `transforms/build/`, versioned in git. A second,
hand-maintained statement of the same fact is precisely the pathology this
refactor removes.

## The four kinds

Every record carries a common header (`id`, `path`, `type`, `tier`, `bytes`,
`sha256`) plus exactly one discriminated block:

| kind | means | must carry |
|---|---|---|
| `external` | someone else made it; we fetched it | `url`, `version`, `retrieved`, `license` |
| `licensed` | fetched, but **not redistributable** | acquisition steps, `redistributable: false`, `degraded_mode` |
| `derived` | one of our build transforms made it | `producer`, `inputs`, `parameters`, `determinism` |
| `container` | an image | `reference` + `digest`, or `recipe` + recipe hash |

`licensed` records must state, in `degraded_mode`, exactly what changes when the
item is absent — a fresh machine without a BioCyc licence has to know what it is
losing, not merely that something failed.

`derived` records for Monte Carlo artifacts declare `determinism: seeded`, not
`reproducible`. They are bit-identical only under the same BitGenerator and
thread count. **Verify them by sha256; never by regeneration.**

## Tiers

Tier answers *who is allowed to produce this*, not how big it is.

- `external/` — fetched from upstream
- `external/licensed/` — fetched, never redistributed
- `derived/` — a `transforms/build/` transform is its only legitimate producer
- `validation/` — frozen expected outputs; gates read them, nothing writes them
- `sif/` — images

## Regenerating the container records

```
python transforms/build/_stage/resolve_container_digests.py \
    src/metasmith_libraries provenance/containers
```

Records under `provenance/containers/` are generated; edit the `.oci` resources
in the library, not the records.

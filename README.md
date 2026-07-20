# FabFos

An automated pipeline for resolving inserts from pooled fosmid DNA, rebuilt on
the [metasmith](https://github.com/hallamlab/Metasmith) workflow framework.

FabFos is now a thin front end: it describes the fosmid pipeline as a metasmith
**library** (typed transforms over containers/conda envs) and lets metasmith
plan and execute the workflow.

```bash
fabfos \
    --reads /.../interleaved.fastq.gz --interleaved \
    --background /.../host_background_genome.fasta \
    --vector /.../plasmid_backbone.fasta \
    --endf /.../forward_ends.fasta --endr /.../reverse_ends.fasta \
    --output ./example_out
```

The pipeline resolves: reads → QC/trim (bbduk) → optional host filtering
(minimap2/samtools, forced when `--background` is given) → assembly (megahit) →
non-redundant contigs (length filter + dedup, blast). When `--endf/--endr` are
given, an orthogonal blast step tags those contigs with the fosmid end
mappings; with `--vector`, a pool-size estimate (minimap2/samtools/vsearch) is
also produced. ORF annotation, when requested, runs on the non-redundant
contigs rather than the raw assembly.

## Repository layout

```
setup.py  dev.sh  envs/  conda_recipe/   # project files (src-layout: the package is in src/)
src/
  fabfos/               # THE package: CLI front end + canon.py
  metasmith/            # submodule → Metasmith @ release (the framework)
  metasmith_libraries/  # submodule → MetasmithLibraries @ feat/fabfos (the fosmid transforms)
transforms/{build,run}/ # build = compile data dependencies; run = consume them
containers/<name>/      # only images WE build; public ones are pinned, not vendored
examples/               # worked templates for calling metasmith directly
reports/dag/            # rendered workflow DAGs (from the planner, not hand-drawn)
provenance/             # declarations: where every dependency came from + its sha256
_old/                   # the previous Snakemake implementation, kept for reference
```

`src/metasmith_libraries` is pinned to `feat/fabfos`, which is the branch
carrying the `fosmids/` transform domain and the ECSPr prerequisite chain —
`release` has neither. `src/metasmith` is pinned to `release`, the only branch
its remote publishes. Both also carry a secondary `awm` git remote (the sibling
local bare repos) for push-free local sync during co-development; the commits
they are pinned to are on the GitHub remotes, so a plain
`git clone --recurse-submodules` works.

## Data dependencies

Everything the pipeline reads — MetaNetX, UniRef50, KOfam, ESM-C weights, the
frozen MetaNetX-4.5 pre-bake, the direction ensemble, the frozen nulls, the
validation set — lives in one metasmith `DataInstanceLibrary` at
`.awm/data/ref/` (911 items, ~46 GB), tiered `external/` (786),
`external/licensed/` (5), `derived/` (46) and `validation/` (74). A `sif/` tier
is described in the tier tables but has never been created — see Known gaps.

`provenance/data/_declared.yml` is the hand-edited source of truth; the records
beside it are **generated** and carry each item's size, sha256 and per-file
digests. Build and check the library with:

```bash
python transforms/build/_stage/build_ref_library.py --stage --place --index --verify
python transforms/build/_stage/check_provenance.py    # AC3/AC4 + hash cross-checks
python transforms/build/_stage/check_canon.py         # canon resolves through the manifest
```

On a machine where the sources already exist, items are **hardlinked** into the
library, so it costs zero bytes and the originals keep working unchanged.

`src/fabfos/canon.py` is the one code↔data interface: scalars are plain
constants, and every data path resolves through the library manifest lazily via
a module-level `__getattr__`. Importing it for a scalar needs neither the
library nor metasmith. A missing key raises `CanonError` naming the symbol —
there is deliberately no fallback to an absolute path.

## Develop

```bash
git clone --recurse-submodules <repo>
cd <repo>
./dev.sh --ibase                 # create the `fabfos` conda env (python + metasmith)
./dev.sh -r --plan-only \        # resolve the DAG without executing
    -r reads.fq.gz -i -o out --vector backbone.fna
```

`./dev.sh -r` runs the CLI from source against the sibling metasmith +
library submodules. `./dev.sh -b` bundles the library into the package
(`src/fabfos/_library`) for shipping; `-bp`/`-bc` build the wheel / conda
package.

## Method version

The package version (`src/fabfos/version.txt`) versions the CLI. The **method**
version (`src/fabfos/method_version.txt`, currently `0.3.0`) versions the
composition that decides what a number out of this pipeline *means*: canon's
content and status, the transform library's commit and dirty flag, metasmith's
own version, every container digest, the sha256 of the data library's index,
the type contract, and the planner's domain list. A CLI bugfix is not a new
method; repinning the engine library is one even if no fabfos source changed.

```bash
fabfos --method-version      # 0.3.0+434c471
fabfos --describe-method     # the full hashed document -- diff two to see WHICH part moved
fabfos --require-method 0.3.0+434c471   # fail unless the live method matches
```

Stamping refuses while a container **on the method path** is unresolved. It is
scoped to the path rather than all 44 records, because refusing over `stringtie`
or `sra-tools` — which this method never invokes — would be a false claim in
the other direction. Every digest is still hashed; only the refusal is scoped.

## Status

The core pipeline plans and runs: reads → QC → host filter → assembly →
non-redundant contigs → pool-size estimate. `fabfos --plan-only --dag <path>`
renders it (5 steps).

The **full ECSPr chain now resolves end to end**, reference inserts →
significance, 12 steps:

```
build_ec_bridge · orfcall_inserts · build_uniprot_bridge · clean_lane ·
uniref_lane · kofam_lane · compile_evidence · evidence_weights ·
addition_weights · base_graphs · solve_directed · significance
```

`examples/ecspr_full_dag.py` builds and renders it, and is the worked template
for calling metasmith directly. The rendered DAG is committed at
`reports/dag/ecspr_full.svg`. Planning does not need the staged files to exist,
so it renders on a machine holding none of the 46 GB.

Lane selection is deliberate, not structural: `ecsprNetB`/`ecsprNetA` both
produce `ecspr::base_graphs`, and `ecsprDirected`/`ecsprUndirected` both produce
the axes reports, so loading both members of a pair would let the planner choose
by tiebreak. `library.domains_for()` drops the unselected lane.

### Known gaps

- **12 staged inputs still resolve to absolute paths in the incumbent tree**
  rather than through the data library — `biomass_axes`, `direction_ratios`,
  `metanetx_chem_prop`, `ko_to_mnxr`, the two `rhea2uniprot` tables, the KOfam
  pair, `uniref50_dmnd`, `metanetx_reac_xref`, `evidence_source`, and the
  reference inserts. `examples/ecspr_full_dag.py` prints this list every run.
  Each one is a reason the method is not yet portable, and closing them is the
  next revision's main work.
- **`ecspr::compute_profile` is undeclared**, so `ecspr::ablation_importance`
  cannot be targeted alongside significance.
- **No transform produces `ecspr::frozen_null`.** It is staged from canon's
  explicit file list — never a glob, because the scorer discovers its draw
  sizes by listing that directory, so a stray retired size would silently widen
  the null basis.
- **The `sif/` tier still does not exist**, despite being named in the tier
  tables. Both images are on quay now, so it is a convenience rather than the
  fallback it was written as.
- **The licensed BioCyc PGDBs are absent** from the data tree — only the
  MetaCyc flatfiles survive, so the direction ensemble's curated member cannot
  be rebuilt as-is. `provenance/data/biocyc.pgdbs.yml` states what degrades.
- **`references/kegg/` is empty**, so KEGG-keyed rollups are not reproducible
  from what is on disk.
- **Three on-path transforms use `envs::*.condaenv`** (`orfcall_inserts`,
  `kofam_lane`, `uniref_lane`), so the canonical chain runs under
  `Runtime.MAMBA`. Under a container runtime the file's *contents* are handed
  to the runtime as an image URI. Three further off-path transforms use it too.
- **An ambient `PYTHONPATH` entry can shadow the installed metasmith** with an
  older incompatible copy; run the examples with `PYTHONPATH=src`.

See the legacy pipeline in `_old/` for the original Snakemake implementation
and full argument/output reference.

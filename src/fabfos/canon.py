"""The canonical ECSPr basis, as data.

This module is the only place these values exist. Prose may point at a name
defined here; prose may not restate its value. That rule is the whole point of
the file, and it is checked mechanically by `check_no_transcribed_numbers.py`.

The failure this corrects was not an absence of documentation. There were four
live documents each declaring itself the canonical method, and they contradicted
each other, because nothing checked them against anything. A fifth well-written
document would reproduce the disease exactly. So the countermeasure is not a
better document: it is that there is nothing left to transcribe. Import the name.

Two consumers, two directions:

  * Experiment drivers (`run_experiment.py`) read these values and stage them as
    inputs. Values flow OUT of this module into the engine.
  * Figure generators call the assertion helpers to refuse a table that is not
    the canonical one.

The engine library (`metasmith-libraries/fabfos`) MUST NOT import this module.
It is a consumer of the basis, never a reader of it — the reverse edge, library
depending on a specific experiment, is the disease we are removing. The engine
derives what it needs from its staged inputs.

Status is `STATUS` below, and it is load-bearing: while PROVISIONAL this file
merely *describes what the incumbent claims*. It becomes the definition only when
the parity gate proves the FabFos path reproduces the incumbent. See MIGRATION.md.
"""
from __future__ import annotations

from pathlib import Path

# =====================================================================
# Status
# =====================================================================
# CANONICAL since 2026-07-14, when the parity gate went green on both lanes:
# the engine reproduces the incumbent tables to ~1e-16 (machine epsilon) on every
# joined row, with exact at_floor agreement and the split contigs present.
# Re-provable at any time, in seconds, on CPU:
#
#     python parity/run_parity.py --all
#
# This is a handoff, not just a test. From this line the incumbent tables are a
# FROZEN REFERENT -- read, never authoritative -- and this file is the definition.
# It proves reproduction, not truth: it shows the engine reproduces the incumbent,
# and does not revisit whether the incumbent is right.
STATUS = "CANONICAL"
STATUS_SINCE = "2026-07-14"
STATUS_NOTE = (
    "Parity gate green on both lanes at ~1e-16 vs the incumbent; "
    "the incumbent is now a frozen referent. See MIGRATION.md."
)

# =====================================================================
# Roots -- resolved through the data-dependency library, not hardcoded
# =====================================================================
# Every reference path below is a KEY into `.awm/data/ref/`, the metasmith
# DataInstanceLibrary built by `transforms/build/_stage/build_ref_library.py`.
# Nothing here is an absolute path to this machine.
#
# Resolution is lazy (PEP 562 module-level __getattr__) and memoized, so
# importing canon for a scalar -- `canon.KMAX`, `canon.DRAW_SIZES` -- costs
# nothing and works with no library present and no metasmith installed. That
# matters: the paper worktrees import this module for its constants and do not
# have metasmith.
#
# Failure is LOUD. A missing key raises CanonError naming both the canon symbol
# and the key it wanted. There is deliberately NO fallback to an absolute path:
# a fallback would make the whole migration untestable, because everything would
# keep working whether or not the library was correct.
import os as _os

_LIB_ENV = "FABFOS_REF"
_DEFAULT_LIB = Path(__file__).resolve().parents[2] / ".awm" / "data" / "ref"


class CanonError(RuntimeError):
    """A canon symbol could not be resolved through the data library."""


def library_root() -> Path:
    root = Path(_os.environ.get(_LIB_ENV, _DEFAULT_LIB))
    if not (root / "_metadata").is_dir():
        raise CanonError(
            f"no data library at [{root}]. Set ${_LIB_ENV}, or build it with\n"
            f"    python transforms/build/_stage/build_ref_library.py --stage --place --index"
        )
    return root


# -----------------------------------------------------------------
# Experiment roots -- deliberately NOT library items
# -----------------------------------------------------------------
# Run results on actual use cases are out of scope for the data-dependency
# library: they are outputs of a particular experiment, not dependencies of the
# method. But canon still has to name them, so they get env-overridable roots
# rather than the hardcoded absolutes they used to be. Point $FABFOS_DATA at
# your own run tree and the ORFS_FAA / EVIDENCE_TABLE / INSERTS_FNA / AXES_TSV
# constants follow.
DATA = Path(_os.environ.get("FABFOS_DATA", "/home/tony/agentic_workspace/data/scadc"))
INCUMBENT_ROOT = Path(_os.environ.get(
    "FABFOS_INCUMBENT",
    "/home/tony/agentic_workspace/projects/scadc/metabolic-modelling"
    "/main/metabolic-modelling/04_reaction_network",
))
INCUMBENT_CACHE = INCUMBENT_ROOT / "cache"
# The hand-made byte-copy backup. The live cache carries retired draw sizes and
# an un-suffixed draws file left by the overwrite incident, so the live cache is
# never a safe source; this directory is. It is still not safe to GLOB -- it also
# holds retired sizes. Curate by explicit list (see FROZEN_NULL_FILES).
INCUMBENT_K1000 = INCUMBENT_CACHE / "K1000"
ENGINE_LIB = Path(_os.environ.get(
    "FABFOS_ENGINE_LIB",
    "/home/tony/agentic_workspace/projects/metasmith-libraries/fabfos",
))

_manifest_cache: dict | None = None


def _manifest() -> dict:
    """path -> type name, read once.

    Prefers metasmith's loader; falls back to reading _metadata/index.yml
    directly so a worktree without metasmith can still resolve paths.
    """
    global _manifest_cache
    if _manifest_cache is not None:
        return _manifest_cache
    root = library_root()
    try:
        from metasmith.models.libraries import DataInstanceLibrary
        lib = DataInstanceLibrary.Load(root)
        _manifest_cache = {str(k): v for k, v in lib.manifest.items()}
    except Exception:
        import yaml
        index = root / "_metadata" / "index.yml"
        if not index.exists():
            raise CanonError(f"library at [{root}] has no _metadata/index.yml")
        raw = yaml.safe_load(index.open()) or {}
        # index.yml IS the manifest -- a flat `path: type` mapping written in
        # YAML's explicit-key form. It has no `manifest:` wrapper, so asking
        # for one yields {} and every symbol then fails as "not in the
        # manifest", pointing at the declaration instead of at this reader.
        # Accept the wrapper if a future writer adds one; otherwise take the
        # document itself.
        man = raw.get("manifest", raw) if isinstance(raw, dict) else {}
        _manifest_cache = {
            str(k): (v["type"] if isinstance(v, dict) else v)
            for k, v in man.items()
        }
    # An EMPTY manifest is a broken read, never a legitimately empty library:
    # every caller is asking for a path that must exist. Failing here names the
    # real fault; failing later names an innocent symbol.
    if not _manifest_cache:
        raise CanonError(
            f"library at [{root}] resolved an EMPTY manifest. The library is "
            f"not built, or _metadata/index.yml is not in the expected "
            f"`path: type` form. This is a reader/library fault, not a bad "
            f"symbol -- do not chase the declaration."
        )
    return _manifest_cache


def _resolve(symbol: str, key: str) -> Path:
    man = _manifest()
    root = library_root()
    if key in man:
        return root / key
    # a directory item is not itself a manifest key -- its members are
    prefix = key.rstrip("/") + "/"
    if any(k.startswith(prefix) for k in man):
        return root / key
    raise CanonError(
        f"canon.{symbol} wants library key [{key}], which is not in the manifest "
        f"at [{root}]. Either the library is stale (rebuild it) or the "
        f"declaration in provenance/data/_declared.yml is wrong."
    )


# canon symbol -> library key. The ONE table mapping code to data.
_PATHS: dict[str, str] = {
    "REFERENCE_ROOT":              "derived/mnxref-4_5",
    "REFERENCE_GRAPH_DIR":         "derived/mnxref-4_5/graph",
    "REFERENCE_SOLVE_DIR":         "derived/mnxref-4_5/solve",
    "REFERENCE_SOLVE_DIRECTED_DIR":"derived/mnxref-4_5/solve_directed",
    "REFERENCE_DIRECTION":         "derived/mnxref-4_5/direction.parquet",
    "REFERENCE_MANIFEST":          "derived/mnxref-4_5/MANIFEST.json",
    "REFERENCE_ATOM_PAIRS":        "derived/mnxref-4_5/atom_pairs.parquet",
    "REFERENCE_LEDGER":            "derived/mnxref-4_5/closure_ledger.parquet",
    "REFERENCE_NULL_UNDIRECTED_DIR": "derived/mnxref-4_5/null",
    "REFERENCE_NULL_DIRECTED_DIR": "derived/mnxref-4_5/null_directed",
    "AXES_TESTABLE_JSON":          "derived/mnxref-4_5/solve/axes_testable.json",
    "SIGNIFICANCE_DIR":            "validation/significance",
    "VALIDATION_DIR":              "validation/dual_network",
    "REAC_PROP":                   "external/metanetx/4.5/reac_prop.tsv",
    "CHEM_PROP":                   "external/metanetx/4.5/chem_prop.tsv",
    "CHEM_XREF":                   "external/metanetx/4.5/chem_xref.tsv",
    "REAC_XREF":                   "external/metanetx/4.5/reac_xref.tsv",
    "DIR_REAC_PROP":               "external/metanetx/4.5/reac_prop.tsv",
    "DIR_CHEM_PROP":               "external/metanetx/4.5/chem_prop.tsv",
    "DIR_CHEM_XREF":               "external/metanetx/4.5/chem_xref.tsv",
    "NETA_GEM":                    "external/gem/iECDH10B_1368.json",
    "DIR_DATA":                    "derived/direction",
    "DIR_TABLE":                   "derived/direction/direction_annotation.parquet",
    "DIR_CALIBRATION":             "derived/direction/calibration.parquet",
    "DIR_CURATED":                 "derived/direction/curated_per_mnxr.parquet",
    "METACYC_FLATFILES":           "external/licensed/metacyc26_flatfiles",
    "UNIREF50_DMND":               "derived/uniref50/uniref50.dmnd",
    "ESMC_WEIGHTS":                "external/esmc/esmc_600m.tgz",
    "KOFAM_KO_LIST":               "external/kofam/ko_list",
    "KOFAM_PROFILES":              "external/kofam/profiles.tar.gz",
    # Aliases the original file defined by plain assignment. They must route
    # through __getattr__ too: a module-level binding always wins over
    # __getattr__, so leaving the assignment in would silently restore the old
    # absolute path while every test still passed.
    "BIPARTITE_DIR":               "derived/mnxref-4_5/graph",
    "SOLVE_BASE_DIR":              "derived/mnxref-4_5/solve",
    # ---- the X/Y benchmark, v3 ----
    # One SELF-CONTAINED tree: X, the contract shape, the ground truth the key
    # is derived from, the baseline and v1's provenance all live inside v3, so
    # deleting a sibling version cannot break this one. BENCH_V3_Y is the
    # answer key and is declared MISSING until the v3 key is built and frozen;
    # touching it before then raises CanonError naming the symbol, which is the
    # intended refusal -- scoring against an absent key must never quietly
    # produce an empty result.
    "BENCH_V3_ROOT":               "validation/benchmark/v3",
    "BENCH_V3_OBSERVATIONS":       "validation/benchmark/v3/observations",
    "BENCH_V3_DECISIONS":          "validation/benchmark/v3/decisions",
    "BENCH_V3_X":                  "validation/benchmark/v3/X",
    "BENCH_V3_CONTRACT":           "validation/benchmark/v3/contract",
    "BENCH_V3_GROUND_TRUTH":       "validation/benchmark/v3/ground_truth",
    "BENCH_V3_BASELINE":           "validation/benchmark/v3/baseline",
    "BENCH_V3_V1_PROVENANCE":      "validation/benchmark/v3/v1",
    "BENCH_V3_BASE_GRAPHS":        "validation/benchmark/v3/base_graphs",
    "BENCH_V3_UNIVERSE":           "validation/benchmark/v3/universe",
    "BENCH_V3_Y":                  "validation/benchmark/v3/Y",
}

# Declared, but absent from every machine we have looked at. Named so the
# failure says what is missing and how to get it, instead of FileNotFoundError
# on a path nobody recognises.
_UNAVAILABLE: dict[str, str] = {
    "DIR_METACYC_PGDB": (
        "metacyc26.pgdb is not on this machine and was not found anywhere in the "
        "data tree. It is a licensed BioCyc artifact (subscription required, not "
        "redistributable). Only metacyc26_flatfiles survives -- see "
        "canon.METACYC_FLATFILES and provenance/data/biocyc.pgdbs.yml."
    ),
    "DIR_ECOCYC_PGDB": (
        "ecocyc26.pgdb is not on this machine and was not found anywhere in the "
        "data tree. See provenance/data/biocyc.pgdbs.yml for acquisition and for "
        "what the direction ensemble loses without it."
    ),
}


def __getattr__(name: str) -> Path:
    """PEP 562: resolve data paths on first touch, never at import."""
    if name in _UNAVAILABLE:
        raise CanonError(f"canon.{name}: {_UNAVAILABLE[name]}")
    if name in _PATHS:
        return _resolve(name, _PATHS[name])
    if name == "REFERENCE_NULL_DIR":
        return __getattr__(
            "REFERENCE_NULL_DIRECTED_DIR" if _DIRECTED else "REFERENCE_NULL_UNDIRECTED_DIR"
        )
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__() -> list:
    return sorted(list(globals()) + list(_PATHS) + list(_UNAVAILABLE) + ["REFERENCE_NULL_DIR"])

# =====================================================================
# Solver orientation
# =====================================================================
# The canonical solver ORIENTATION. Undirected until 2026-07-18; flipped to "directed"
# after the fixed diode solver (softplus-smoothed rectification, engine 933d93e -> a
# strictly convex energy with a unique, start-independent minimiser) regenerated the
# directed null on the honest reference graph and it cleared the replacement gates: the
# symmetric-limit parity (directed force-noop reproduces the undirected reference to ~1e-6,
# directed/directed_parity.py), the engine self-tests (worst |d|~1e-15), and observed<->null
# lockstep (identical testable-axis sets per element). The undirected solve + null are
# RETAINED as the frozen symmetric-limit referent the parity gates join against -- they are
# no longer the graph+solver the canonical solve runs on. Flip this one constant to revert;
# every resolver below follows it.
CANONICAL_ORIENTATION = "directed"
_ORIENTATIONS = ("undirected", "directed")
assert CANONICAL_ORIENTATION in _ORIENTATIONS
_DIRECTED = CANONICAL_ORIENTATION == "directed"

# The per-reaction direction ratios (exp(dG'/RT), MNXR-keyed) the directed solve rectifies
# each edge with. Frozen beside the reference; ratio 1.0 == no evidence == reversible ==
# undirected parity. Baked into the frozen directed reports/null below -- pinned here as the
# provenance of which direction table produced the canonical directed solve, and hash-checked
# in assert_canonical_reference() so a silent regenerate cannot change the canonical answer.
# [resolved through the library] REFERENCE_DIRECTION = REFERENCE_ROOT / "direction.parquet"
REFERENCE_DIRECTION_SHA256 = (
    "c80009055601e64342ba56aaf48547f609f2621c47fa4f495abff8c6fd6a8fc5"
)

# The canonical null: the reference-GRAPH null (reference/build_reference_null.py for the
# undirected referent; directed/build_directed_null.py for the directed canonical), the
# matched pair to reference_axes_report -- same graph AND same orientation, so delta_obs and
# the null it is scored against are never on different solvers. Same draws/ORF weights across
# both orientations (graph-independent, reused verbatim). Curated by the explicit
# FROZEN_NULL_FILES list below -- never a glob.
# [resolved through the library] REFERENCE_NULL_UNDIRECTED_DIR = REFERENCE_ROOT / "null"
# [resolved through the library] REFERENCE_NULL_DIRECTED_DIR = REFERENCE_ROOT / "null_directed"
# [resolved through the library] REFERENCE_NULL_DIR = REFERENCE_NULL_DIRECTED_DIR if _DIRECTED else REFERENCE_NULL_UNDIRECTED_DIR

# =====================================================================
# Scorer
# =====================================================================
SCORER = "sig_mix"
SCORER_DESC = "full-mixture survival function"
# Retired scorers. Named here so tooling can refuse them by name rather than by
# a human remembering which of the sibling tables is current.
RETIRED_SCORERS = ("sig_negbin", "sig_emp")

# =====================================================================
# Fosmid basis
# =====================================================================
FOSMID_BASIS = 199
ORFS_FAA = DATA / "fabfos_2026" / f"orfs_{FOSMID_BASIS}.faa"
# Split contigs. Present in the basis and dropped silently by a `\w`-based ORF
# regex, because `\w` excludes the dot. Their presence in a scored table is a
# positive check that the ORF counter is the rsplit form.
SPLIT_CONTIGS = ("C00310.A", "C00310.B", "C00708.A", "C00708.B")

# The figure basis is a filter applied at figure time, not a second table:
# score once on FOSMID_BASIS and subset. BH-q is not load-bearing (ranking is by
# effect size, which is per-row and basis-independent), so subsetting is safe.
FIGURE_BASIS_OPEN = 132

# =====================================================================
# Null draws
# =====================================================================
# The sizes the SCORER interpolates across. NOT the sizes the incumbent null
# GENERATOR declares -- those two forked inside one directory, which is one of
# the failures this work exists to remove.
DRAW_SIZES = (14, 25, 30, 35, 43)
# The grid is the ORF-count distribution of the 199 canonical inserts sampled at the
# p10/p30/p50/p70/p90 percentiles (orfs_199.faa ORF calls; median 30 ORFs, 0.96 ORF/kb).
# Superseded the earlier (14,28,34,42,51) when the directed model was frozen: those knots
# ran to n_orf~p95 (irregular spacing p10/p39/p62/p88/p95); this grid is evenly spaced in
# ORF-count rank so the scorer interpolates on a uniform lattice over the real insert range.
# What the incumbent null GENERATOR declares, against which the scorer above had
# already moved on. Recorded so the fork is visible as data rather than as a
# discovery someone has to make twice.
RETIRED_DRAW_SIZES = (21, 28, 34, 42, 51, 56)
DRAW_K = 1000
DRAW_SEED = 42
DRAW_N_ORFS = 30
MIN_CONTIG_ORFS = 30

STYLES = ("A", "B", "D", "E")
PRIMARY_STYLE = "D"
STYLE_DESC = {
    "A": "uniform over all metagenome ORFs",
    "B": "uniform over ORFs on contigs above the contig-ORF minimum",
    "D": "contiguous adjacent ORFs on one contig (operon-like)",
    "E": "uniform over evidence-bearing ORFs only",
}

LANES = ("reff", "ieff")
CANONICAL_LANE = "ieff"

# The explicit curated list. Never glob the cache or the backup: both hold
# retired sizes beside the canonical ones. The filename stem follows the canonical
# orientation: the directed null is "{lane}_null_directed_N{n}", the undirected referent
# "{lane}_null_canonical_N{n}" -- so a directed run can never silently score against the
# undirected null (different files, not a flag on the same file).
_NULL_STEM = "directed" if _DIRECTED else "canonical"
FROZEN_NULL_FILES = tuple(
    f"{lane}_null_{_NULL_STEM}_N{n}.tsv" for lane in LANES for n in DRAW_SIZES
)
FROZEN_DRAWS_FILES = tuple(
    f"null_canonical_N{n}_draws.parquet" for n in DRAW_SIZES
)

# =====================================================================
# Fitter / scorer knobs
# =====================================================================
ELEMENTS = ("C", "N", "S", "P")
KMAX = 4          # pinned; an implicit cap makes the incumbent comparison dishonest
FIT_SEED = 0
EPS = 1e-12
Q_THRESHOLD = 0.05

# =====================================================================
# Axes
# =====================================================================
AXIS_SET = "set4"
RETIRED_AXIS_SETS = ("set2", "set2cat")
# Lives in the PUBLISH tree while being read as a pipeline INPUT. That inversion
# is logged in MIGRATION.md; the path is recorded here so nothing has to guess it.
AXES_TSV = DATA / "figures" / "publish" / "03_model" / "ecspr_scadc" / f"biomass_edges_{AXIS_SET}.tsv"
# The axis DEFINITIONS (source/sink species per axis). Graph-INDEPENDENT: it is set4
# itself, so it is exactly AXES_N regardless of which graph the solve runs on. Assert
# canonical axes against THIS.
AXES_JSON = INCUMBENT_CACHE / f"biomass_dag_axes_{AXIS_SET}.json"
# The TESTABLE subset: which of the AXES_JSON axes have both endpoints in the solve
# graph's LCC. Graph-DEPENDENT, so it moves with the canonical graph -- on the honest
# reference graph one carbon axis (a phospholipid endpoint the star reached only via a
# fabricated transit) is no longer testable, so this is a strict subset of set4 and is
# NOT asserted equal to AXES_PER_ELEMENT. It lives beside the reference solve it
# describes.
# [resolved through the library] AXES_TESTABLE_JSON = REFERENCE_SOLVE_DIR / "axes_testable.json"
AXES_N = 79
AXES_PER_ELEMENT = {"C": 40, "N": 23, "S": 7, "P": 9}

# =====================================================================
# Evidence
# =====================================================================
CLEAN_FLOOR = 0.01   # F1-optimal; CLEAN does not abstain, so an unguarded lane over-nominates

# The compiled per-(source, orf, mnxr) evidence behind the canonical basis.
# Two byte-identical copies exist -- one in the code tree, one here in the data
# tree. Prefer this one: data belongs in the data tree, and a pipeline input read
# out of a code checkout is the same inversion as the axis table below.
EVIDENCE_TABLE = (DATA / "fabfos_2026_199" / "ecspr_clean" / "evidence_network"
                  / "evidence_table_clean.parquet")

# The nucleotide basis. Note the filename does NOT record the count -- the sibling
# that does is the retired subset. Assert the count; do not read the name.
INSERTS_FNA = DATA / "fabfos_2026" / "putative_inserts.fna"

# The reactions each fosmid injects onto the host base. Derived from the evidence
# table; staged rather than recomputed, so that a Network A run and a Network B run
# differ ONLY in the base and the comparison between them is about the host.
ADDITION_WEIGHTS = (DATA / "fabfos_2026_199" / "ecspr_clean" / "evidence_network"
                    / "fosmid_addition_weights.pkl")

# =====================================================================
# Network A -- the curated-GEM host
# =====================================================================
# The second, independent reconstruction of the host. Network B induces the host
# from the 4-lane annotation of its ORFs; Network A crosswalks a curated GEM's
# reactome to current MNXR and induces it with UNIFORM conductance. They share no
# derivation path, which is the entire point: a finding that holds on both rests on
# neither.
#
# The sibling curated model is K-12 (iML1515) and is a DIFFERENT STRAIN with a
# different genotype. It is the validation lane's referent, not the SCADC host.
# Experiments assert the model id from inside the JSON rather than trusting the
# filename -- the filename is what drifted, everywhere else in this file.
# [resolved through the library] NETA_GEM = DATA / "metabolic_modelling" / "networks" / "iECDH10B_1368.json"
NETA_GEM_STRAIN = "DH10B"
# [resolved through the library] REAC_XREF = DATA / "references" / "metanetx" / "reac_xref.tsv"

# =====================================================================
# The frozen AAM + direction reference
# =====================================================================
# ECSPr's atom-atom mapping and directionality are NOT dynamic: every reaction in
# any network is a MetaNetX id, so AAM(mnxr) and direction(mnxr) are static
# functions of that id. They are pre-baked ONCE over the whole MetaNetX universe as
# a frozen, version-pinned, MNXR-keyed asset (main/fabfos/reference/reference.py),
# driven to 100% adjudication by the closure ledger. The method module owns the
# schema, verdict vocabulary, closure, and the reac_prop content pin; it is the sole
# authority on those, exactly as this file is on the basis.
#
# What canon adds is the EXPERIMENT->UNIVERSE binding. reference.assert_reference
# proves the manifest is SELF-consistent (its recorded reac_prop hash equals the
# current file's), but a fresh rebuild against a DIFFERENT MetaNetX release is
# self-consistent too. This basis was validated against exactly one frozen universe;
# canon records that universe's reac_prop hash so a swap underneath fails at the
# experiment boundary, not silently. The value lives in the reference MANIFEST.json
# (written by build_ledger.py); it is data, not prose, and is pinned here so a
# consumer proves the binding with assert_canonical_reference().
REFERENCE_REAC_PROP_SHA256 = (
    "8582cc187d03ce127f8e914f8f298282ac9036f1448918117af044ac55b980db"
)

# =====================================================================
# Reference graphs
# =====================================================================
# Derived from ELEMENTS, deliberately. The dir also holds a bipartite map for an
# element outside the basis, so a glob of `mnx_bipartite_*.pkl` silently widens the
# run. Same failure as globbing the null cache, different directory.
BIPARTITE_FILES = tuple(f"mnx_bipartite_{e}.pkl" for e in ELEMENTS)
# The per-element universe bipartite the fosmid-addition map projects onto. Now the
# honest reference graph (see REFERENCE_GRAPH_DIR); the incumbent cache is a frozen
# referent, no longer the graph the solve runs on.
# [resolved through the library] BIPARTITE_DIR = REFERENCE_GRAPH_DIR

# The per-element host base graphs the solve runs on. `ADDITION_WEIGHTS` above is
# the other half; there is deliberately no second name for it here, because a
# duplicate name for one value is the disease this file exists to cure. These are the
# reference-fed bases beside the frozen reference solve.
# [resolved through the library] SOLVE_BASE_DIR = REFERENCE_SOLVE_DIR
BASE_GRAPH_FILES = tuple(f"base_{e}.pkl" for e in ELEMENTS)

# =====================================================================
# The solver gate
# =====================================================================
# WHY THIS EXISTS. Neither existing gate runs the solver. `run_parity.py` feeds the
# scorer `delta_obs` from `{lane}_axes_report.tsv` via --report-dir, and every
# experiment spec stages `ecspr::{reff,ieff}_axes_report` as an INPUT -- scadc_main
# says so outright ("the staged, canonical solve"). Both cover (delta_obs, nulls) ->
# p/q. The solve itself was gated only by a prose claim in transforms/ecspr/solve.py
# ("verified ... max|abs| 1.4e-8") -- untested, and LOOSER than PARITY_TOL.
#
# WHY IT GATES AGAINST A DENSE REBUILD RATHER THAN THE INCUMBENT TABLE. Measured:
# the engine solver agrees with a dense ground-truth rebuild to ~1e-14, while the
# incumbent `axes_report` agrees only to ~1e-8. `g_base` matches at ~1e-16, so the
# graph is identical and the divergence is entirely in the Woodbury update -- and it
# is the INCUMBENT that carries it. Gating the solver against the incumbent would
# therefore pin the accurate implementation to the inaccurate one. The referent is
# the rebuild: no Woodbury, no clique elimination, no shared factorisation.
SOLVER_GATE_TOP_CELLS = 3      # per (element, axis): the largest observed deltas
SOLVER_GATE_RANDOM_CELLS = 2   # per (element, axis): typical rows, not just extremes
SOLVER_GATE_AXES = 2           # axes sampled per element
SOLVER_GATE_SEED = 7

# =====================================================================
# Compute
# =====================================================================
# Staged as a content-hashed input so it enters the task hash. See MIGRATION.md
# on why context.params cannot carry this.
COMPUTE_CPU = "device: cpu\ndtype: float64\n"
COMPUTE_GPU = "device: cuda\ndtype: float64\n"

# =====================================================================
# Significance table schema
# =====================================================================
# The full-mixture SF interpolates across flanking anchors, so there is no single
# matched size. `matched_N` belonged to the retired nearest-size scorer and its
# presence in a table is positive evidence that the table is stale.
SIG_COLUMNS = (
    "fosmid", "element", "axis_id", "n_orfs", "N_lo", "N_hi", "w", "null",
    "delta_obs", "p", "p_emp", "k_lo", "k_hi", "n_bg", "p_floor", "at_floor",
    "q", "survives",
)
SIG_STALE_COLUMNS = ("matched_N",)

# =====================================================================
# Tolerances
# =====================================================================
# Committed BEFORE the runs they gate. A threshold chosen after seeing the number
# is a rationalization, not a gate.
PARITY_TOL = 1e-9        # library scorer vs the incumbent table, CPU vs CPU
GPU_CPU_TOL = 1e-6       # GPU float64 vs CPU float64; agrees closely, not bit-exactly

# =====================================================================
# Direction annotator
# =====================================================================
# Per-reaction conductance directionality: every base-graph reaction gets a ratio
# g_reverse/g_forward = exp(dG'/RT), fused from three members spanning method
# families -- eQuilibrator (measured + group contribution), dGbyG (learned GNN),
# and BioCyc REACTION-DIRECTION (curated physiology, orientation-aligned to MNXR).
# No evidence -> dG'=0 -> ratio 1.0 -> reversible -> parity with the undirected
# model. The ratio is the median transform exp(mu_eff/RT); the sign is expressed
# in MNXR equation orientation. All knobs live here, committed before the run.
import math as _math

DIR_R = 8.314e-3                        # kJ/mol/K
DIR_T = 298.15                          # K
DIR_RT = DIR_R * DIR_T                  # the ratio's natural scale, ~2.48 kJ/mol
DIR_DECADE = DIR_RT * _math.log(10.0)   # one decade of conductance, ~5.71 kJ/mol

# Floors, committed BEFORE the run. TAU_SHARED is the TECRDB common-mode error the
# correlated eQ/dGbyG pair's spread cannot see; it floors a PREDICTED thermo vote
# (group-contribution arm / dGbyG) and the fused pair, so two correlated predictors
# never vote as two independent. TAU_CUR_FLOOR: a curated category alone resolves no
# better than one decade. Both are physical, not tuned.
DIR_TAU_SHARED = DIR_DECADE
DIR_TAU_CUR_FLOOR = DIR_DECADE
DIR_S_MEAS_FLOOR = 0.1                  # kJ/mol; numerical only -- a real measurement
                                        # is trusted at its own sigma
DIR_SIGMA_CEILING = 100.0               # kJ/mol; a wider eQ uncertainty is no
                                        # information -> the reaction is eQ-silent

# The reversible-default prior width = robust marginal spread of measured dG' on the
# eQuilibrator reactant-contribution arm (1.4826*MAD). ESTIMATOR committed here; the
# VALUE is frozen from the calibration run that produced it (464 measured reactions,
# marginal median 0.000 -> no orientation offset). It must fall in the plausibility
# band or it is a finding, not a constant. The robust spread runs BELOW the
# outlier-inflated std, i.e. toward more shrinkage / more reversible -- the safe side.
DIR_SIGMA_0 = 9.505                     # kJ/mol
DIR_SIGMA_0_BAND = (5.0, 40.0)          # outside => stop, it is a finding

# The ratio must stay a FINITE two-way conductance ratio -- never a hard one-way
# gate (the standing ruling). A handful of macromolecular/polymer reactions carry a
# genuine |dG'| of thousands of kJ/mol, whose exp() underflows to 0.0 (an infinite
# gate). |dG'| is clamped to this bound: beyond ~the steepest realistic single-
# reaction drive in metabolism, the flux-force is saturated, and clamping keeps the
# ratio finite and > 0 (~3e-18 .. 3e17). Physical bound, committed independent of the
# data; clamped rows are flagged, not hidden.
DIR_DG_CLAMP = 100.0                     # kJ/mol

# The five curated REACTION-DIRECTION values, in MNXR orientation. A sixth token
# would be a KeyError at the aligner, not a silent default (which is how the ~7%
# right-to-left corpus would otherwise invert).
DIR_CATEGORIES = ("PHYSIOL-LEFT-TO-RIGHT", "LEFT-TO-RIGHT", "REVERSIBLE",
                  "PHYSIOL-RIGHT-TO-LEFT", "RIGHT-TO-LEFT")

# Emitted schema. dir_tier > 0 means "carries directional information", NOT "usable":
# every row is usable and ratio==1.0 (reversible) is a real physical statement, so a
# consumer must read ALL rows. A `df[df.dir_tier>0]` filter would silently drop every
# reversible reaction -- turning "default reversible" into "default absent".
DIR_COLUMNS = ("mnxr", "dG_prime", "sigma", "ratio",
               "dir_tier", "dir_method", "dir_confidence")

# Inputs and artifacts. The curated member reads the BioCyc pgdbs read-only; the
# MetaNetX crosswalks are shared with the annotation lanes. REAC_XREF is above.
# [resolved through the library] DIR_METACYC_PGDB = Path(
#     "/home/tony/agentic_workspace/projects/self-improvement/main/staging/metacyc26.pgdb")
# [resolved through the library] DIR_ECOCYC_PGDB = Path(
#     "/home/tony/agentic_workspace/projects/self-improvement/main/staging/ecocyc26.pgdb")
# [resolved through the library] DIR_REAC_PROP = DATA / "references" / "metanetx" / "reac_prop.tsv"
# [resolved through the library] DIR_CHEM_PROP = DATA / "references" / "metanetx" / "chem_prop.tsv"
# [resolved through the library] DIR_CHEM_XREF = DATA / "references" / "metanetx" / "chem_xref.tsv"
# [resolved through the library] DIR_DATA = DATA / "direction"                    # produced artifacts (shared data tree)
# [resolved through the library] DIR_TABLE = DIR_DATA / "direction_annotation.parquet"   # the per-reaction annotator
# [resolved through the library] DIR_CALIBRATION = DIR_DATA / "calibration.parquet"
# [resolved through the library] DIR_CURATED = DIR_DATA / "curated_per_mnxr.parquet"

# =====================================================================
# The canonical solve reports (the reference graph)
# =====================================================================
def reference_axes_report(lane: str = CANONICAL_LANE) -> Path:
    """The canonical solve the experiment stages as delta_obs.

    The frozen reference solve on the honest reference graph, in the canonical
    ORIENTATION (CANONICAL_ORIENTATION): the DIRECTED solve
    (REFERENCE_SOLVE_DIRECTED_DIR/{lane}_axes_report.tsv) since 2026-07-18, the undirected
    solve (REFERENCE_SOLVE_DIR/...) before. Built at 0 pending and hash-pinned via
    assert_canonical_reference(); read-only frozen data, never rebuilt by a consumer. This
    replaced incumbent_axes_report() as the staged solve when the reference graph became
    canonical (2026-07-17); the orientation flipped to directed once the fixed diode
    solver's null cleared its gates. For the undirected symmetric-limit referent the parity
    gates join against, use undirected_axes_report().
    """
    _require_lane(lane)
    d = REFERENCE_SOLVE_DIRECTED_DIR if _DIRECTED else REFERENCE_SOLVE_DIR
    return d / f"{lane}_axes_report.tsv"


def undirected_axes_report(lane: str = CANONICAL_LANE) -> Path:
    """The frozen UNDIRECTED reference solve -- the symmetric-limit referent the directed
    parity gates join against (directed force-noop must reproduce THIS, cell for cell, to
    ~1e-6). Once CANONICAL_ORIENTATION is 'directed' this is NO LONGER the staged canonical
    solve (that is reference_axes_report()); it is the fixed referent that proves the
    directed model degrades to the undirected one wherever direction is unknown."""
    _require_lane(lane)
    return REFERENCE_SOLVE_DIR / f"{lane}_axes_report.tsv"


# =====================================================================
# The incumbent referent (frozen, historical -- NOT the canonical graph)
# =====================================================================
# Kept so the SCORER-parity gate can still join the frozen incumbent inputs to the
# frozen incumbent outputs (a pure regression of the mixture-SF port, independent of
# which graph is canonical). These are no longer the graph the solve runs on.
def incumbent_sig_table(lane: str = CANONICAL_LANE) -> Path:
    """The incumbent significance table the SCORER-parity gate joins against."""
    _require_lane(lane)
    return INCUMBENT_ROOT / "reff" / f"{SCORER}_{lane}.tsv"


def incumbent_axes_report(lane: str = CANONICAL_LANE) -> Path:
    """The incumbent star-graph solve. Frozen historical referent; use
    reference_axes_report() for the canonical (reference-graph) solve."""
    _require_lane(lane)
    return INCUMBENT_CACHE / f"{lane}_axes_report.tsv"


def _require_lane(lane: str) -> None:
    if lane not in LANES:
        raise ValueError(f"unknown lane {lane!r}; expected one of {LANES}")


# =====================================================================
# Assertion helpers
# =====================================================================
# The two checks a consumer runs before trusting a table. They exist so that a
# stale input fails loudly at the point of use, instead of silently producing a
# figure that asserts the opposite of its own annotations -- which has already
# shipped once.

class CanonError(AssertionError):
    """A table is not the canonical one. Always loud, never a warning."""


def _read(table):
    import pandas as pd
    if hasattr(table, "columns"):
        return table
    p = Path(table)
    if not p.exists():
        raise CanonError(f"table does not exist: {p}")
    return pd.read_csv(p, sep="\t")


def assert_canonical_significance(table, *, basis: int = FOSMID_BASIS,
                                  lane: str | None = None):
    """Refuse a significance table that is not the canonical scorer's output.

    Pass a path or a DataFrame. Returns the DataFrame so it can wrap a read:

        sig = canon.assert_canonical_significance(path)

    Checks, in the order that gives the most useful error first:
      * the filename does not name a retired scorer;
      * the schema is the mixture scorer's, not the retired nearest-size one;
      * the split contigs survived the ORF counter;
      * the fosmid count is the expected basis.
    """
    if isinstance(table, (str, Path)):
        name = Path(table).name
        for retired in RETIRED_SCORERS:
            if name.startswith(retired):
                raise CanonError(
                    f"{name} is the retired {retired} scorer's table. "
                    f"The canonical scorer is {SCORER} ({SCORER_DESC}); "
                    f"see canon.incumbent_sig_table()."
                )
    df = _read(table)

    stale = [c for c in SIG_STALE_COLUMNS if c in df.columns]
    if stale:
        raise CanonError(
            f"table carries retired column(s) {stale} -- that schema belongs to "
            f"the nearest-size scorer, which the {SCORER_DESC} replaced. "
            f"Regenerate against the canonical scorer."
        )
    missing = [c for c in SIG_COLUMNS if c not in df.columns]
    if missing:
        raise CanonError(
            f"table is missing canonical column(s) {missing}. "
            f"Expected canon.SIG_COLUMNS."
        )

    if lane is not None:
        _require_lane(lane)

    fosmids = set(df["fosmid"].astype(str))
    missing_splits = [c for c in SPLIT_CONTIGS if c not in fosmids]
    if missing_splits:
        raise CanonError(
            f"split contigs absent from the table: {missing_splits}. These are "
            f"dropped silently by a `\\w`-based ORF-id regex (the dot is not a "
            f"word character); their absence means the ORF counter is the wrong "
            f"one, not that the data lacks them."
        )

    if len(fosmids) != basis:
        raise CanonError(
            f"table covers {len(fosmids)} fosmids; expected {basis}. "
            f"A different count means a different basis -- subset at figure time "
            f"rather than scoring a second table."
        )
    return df


def assert_canonical_axes(table):
    """Refuse an axis table that is not the canonical axis set.

    Accepts the axes TSV, the axes JSON, or the testable JSON (path), or a
    DataFrame of the TSV. Asserts the count and the per-element split rather
    than trusting the filename, because the filename is what drifted.
    """
    if isinstance(table, (str, Path)):
        p = Path(table)
        for retired in RETIRED_AXIS_SETS:
            if f"_{retired}." in p.name or p.name.endswith(f"_{retired}.json"):
                raise CanonError(
                    f"{p.name} is the retired {retired} axis set; the canonical "
                    f"set is canon.AXIS_SET."
                )
        if p.suffix == ".json":
            import json
            data = json.loads(p.read_text())
            # testable json: {element: [axis_id, ...]}; axes json: {axis_id: {...}}
            if data and all(isinstance(v, list) for v in data.values()):
                per_el = {k: len(v) for k, v in data.items()}
                if per_el != AXES_PER_ELEMENT:
                    raise CanonError(
                        f"axis per-element split {per_el} != canon.AXES_PER_ELEMENT "
                        f"{AXES_PER_ELEMENT}"
                    )
                total = sum(per_el.values())
            else:
                total = len(data)
            if total != AXES_N:
                raise CanonError(f"{total} axes; expected canon.AXES_N ({AXES_N})")
            return data

    df = _read(table)
    if len(df) != AXES_N:
        raise CanonError(
            f"axis table has {len(df)} rows; expected canon.AXES_N ({AXES_N}). "
            f"Assert the count, do not trust the filename."
        )
    if "element" in df.columns:
        per_el = df["element"].value_counts().to_dict()
        if per_el != AXES_PER_ELEMENT:
            raise CanonError(
                f"axis per-element split {per_el} != canon.AXES_PER_ELEMENT "
                f"{AXES_PER_ELEMENT}"
            )
    return df


def assert_canonical_direction_table(table, *, n_reactions: int | None = None):
    """Refuse a direction table that cannot be trusted as a per-reaction annotator.

    Pass a path (parquet or tsv) or a DataFrame; returns the DataFrame. Unlike the
    significance idiom, `dir_tier > 0` is NOT a usability filter here -- every row
    is usable and ratio==1.0 (reversible) is a real value, so the checks are about
    completeness and the ratio being a well-formed conductance ratio:
      * the schema is canon.DIR_COLUMNS;
      * exactly one row per reaction (no missing, no duplicate MNXR);
      * ratio is never null -- a no-evidence reaction is ratio 1.0, never absent;
      * ratio is strictly positive (it is exp(dG'/RT));
      * if n_reactions is given, the row count matches (the base-graph basis).
    """
    if hasattr(table, "columns"):
        df = table
    else:
        p = Path(table)
        if not p.exists():
            raise CanonError(f"table does not exist: {p}")
        import pandas as pd
        df = pd.read_parquet(p) if p.suffix == ".parquet" else pd.read_csv(p, sep="\t")

    missing = [c for c in DIR_COLUMNS if c not in df.columns]
    if missing:
        raise CanonError(
            f"direction table missing canonical column(s) {missing}. "
            f"Expected canon.DIR_COLUMNS."
        )
    dups = df["mnxr"][df["mnxr"].duplicated()].unique().tolist()
    if dups:
        raise CanonError(
            f"duplicate MNXR row(s) {dups[:5]}{' ...' if len(dups) > 5 else ''}: the "
            f"annotator is one row per reaction, and a many-to-one curated collapse "
            f"must resolve, not duplicate."
        )
    if df["ratio"].isna().any():
        n = int(df["ratio"].isna().sum())
        raise CanonError(
            f"{n} row(s) carry a null ratio. A reaction with no evidence is ratio "
            f"1.0 (reversible), never null -- a null here means the default-reversible "
            f"limit was skipped, not that direction is missing."
        )
    if not (df["ratio"] > 0).all():
        raise CanonError(
            "ratio must be strictly positive: it is exp(dG'/RT), a conductance ratio, "
            "not a signed quantity."
        )
    if n_reactions is not None and len(df) != n_reactions:
        raise CanonError(
            f"table covers {len(df)} reactions; expected {n_reactions}. The annotator "
            f"is defined on the whole base graph -- a short table means reactions were "
            f"dropped instead of defaulted to reversible."
        )
    return df


def assert_canonical_reference(*, check_hash: bool = True):
    """Refuse an AAM+direction reference that is not the pinned, closed one.

    The schema, closure, MetaNetX-version, and self-consistency checks belong to the
    method module -- delegate them to ``reference.assert_reference`` rather than
    re-listing (a second copy is the drift this file exists to remove). Then add the
    one experiment-level check the method module cannot make: the reference must be the
    EXACT frozen universe this basis was validated against
    (``canon.REFERENCE_REAC_PROP_SHA256``), not merely a self-consistent rebuild
    against some other release.

    ``reference`` is the sibling method module (``main/fabfos/reference``); it is
    experiment-side, NOT the engine library, so importing it here is a consumer reading
    a method pin -- never the forbidden library->experiment edge. Returns the manifest.
    """
    import sys
    ref_dir = Path(__file__).resolve().parent / "reference"
    sys.path.insert(0, str(ref_dir))
    import reference as ref

    man = ref.assert_reference(check_hash=check_hash, require_closed=True)
    if man.get("reac_prop_sha256") != REFERENCE_REAC_PROP_SHA256:
        raise CanonError(
            f"the frozen reference pins universe {man.get('reac_prop_sha256')!r}, but "
            f"this basis expects canon.REFERENCE_REAC_PROP_SHA256. The reference is a "
            f"different MetaNetX universe than the one the basis was validated against; "
            f"rebuild the reference or re-pin the basis, do not trust it."
        )
    # Directed orientation folds a SECOND frozen input into the canonical answer: the
    # per-reaction direction ratios. Pin it the same way as the universe -- a silent
    # regenerate of direction.parquet would change every directed edge without touching the
    # graph, so the reac_prop hash alone cannot catch it.
    if _DIRECTED and check_hash:
        import hashlib
        if not REFERENCE_DIRECTION.exists():
            raise CanonError(
                f"CANONICAL_ORIENTATION is 'directed' but the direction table is missing: "
                f"{REFERENCE_DIRECTION}. The directed solve is scored against a null built "
                f"with these ratios; without it the canonical answer cannot be trusted."
            )
        got = hashlib.sha256(REFERENCE_DIRECTION.read_bytes()).hexdigest()
        if got != REFERENCE_DIRECTION_SHA256:
            raise CanonError(
                f"direction table {REFERENCE_DIRECTION} hashes {got!r}, but this basis "
                f"pins canon.REFERENCE_DIRECTION_SHA256. A different direction table means "
                f"different directed edges than the canonical directed null was built on; "
                f"rebuild the directed solve+null or re-pin, do not trust it."
            )
    return man

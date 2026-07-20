#!/usr/bin/env python3
"""Build src/fabfos/canon.py from the canonical fig-model canon.py.

Mechanical, and re-runnable: takes the newest canon verbatim, strips the two
things that tie it to this machine (the hardcoded DATA roots and the sys.path
import of the sibling reference.py), and routes every reference-data path
through the library manifest instead.
"""
import re
from pathlib import Path

SRC = Path("/home/tony/agentic_workspace/projects/scadc/fig-model/main/fabfos/canon.py")
# Repo-relative, NOT an absolute path to one checkout. This used to name
# .../projects/fabfos/scadc/src/fabfos/canon.py directly, so running the
# generator from any other worktree silently regenerated a DIFFERENT checkout's
# canon and left the current one untouched -- the exact class of drift this
# script exists to end. make_canon.py lives at transforms/build/_stage/.
DST = Path(__file__).resolve().parents[3] / "src" / "fabfos" / "canon.py"

text = SRC.read_text()

HEADER = '''
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
            f"no data library at [{root}]. Set ${_LIB_ENV}, or build it with\\n"
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
'''

# 1. Replace the whole hardcoded-roots + reference-import block. It runs from
#    the "Roots" banner to just before the "Solver orientation" banner.
start = text.index("# =====================================================================\n# Roots")
end = text.index("# =====================================================================\n# Solver orientation")
text = text[:start] + HEADER.strip() + "\n\n" + text[end:]

# 2. Drop every remaining assignment that would shadow a __getattr__ name --
#    a module-level binding always wins over __getattr__, so leaving one in
#    would silently keep the old absolute path.
shadowed = set(re.findall(r"^(\w+)\s*=", HEADER, re.M)) | set(
    re.search(r"_PATHS: dict\[str, str\] = \{(.*?)\n\}", HEADER, re.S).group(0).count("") * []
)
names = set(re.findall(r'^\s*"(\w+)":', HEADER, re.M))
names |= {"REFERENCE_NULL_DIR", "DIR_METACYC_PGDB", "DIR_ECOCYC_PGDB"}

out_lines, skipping = [], False
for line in text.split("\n"):
    m = re.match(r"^([A-Z_][A-Z0-9_]*)\s*=", line)
    if m and m.group(1) in names:
        indent = "# "
        out_lines.append(f"# [resolved through the library] {line.strip()[:100]}")
        skipping = line.rstrip().endswith(("(", "["))
        continue
    if skipping:
        out_lines.append("# " + line)
        if line.rstrip().endswith((")", "]")):
            skipping = False
        continue
    out_lines.append(line)
text = "\n".join(out_lines)

DST.write_text(text)
print(f"wrote {DST} ({len(text.splitlines())} lines)")
print("neutralised:", ", ".join(sorted(names)))

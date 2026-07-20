"""Locate the metasmith library that defines the FabFos transforms.

Resolution order:

1. ``$FABFOS_LIBRARY`` — explicit override (a metasmith library root).
2. A copy bundled inside this package at ``fabfos/_library`` — what a conda
   install ships (populated by ``dev.sh`` at build time).
3. The dev sibling submodule ``src/metasmith_libraries`` — when running from
   a source checkout of the FabFos repo (this package lives at ``src/fabfos``,
   so the submodule is a direct sibling).

The library root is the directory that contains ``data_types/``,
``resources/`` and ``transforms/``.
"""
import os
from pathlib import Path

# transform domains the FabFos pipeline draws from; mirrors the library layout.
DOMAINS = [
    "assembly",
    "fosmids",
    "functionalAnnotation",
    "logistics",
    "metabolicModelling",
    "metagenomics",
    "pangenome",
    "responseSurface",
    "transcriptomics",
    # ECSPr. The solve lives in two domains on purpose: both
    # ecsprUndirected/solve.py and ecsprDirected/solve_directed.py produce
    # ecspr::{reff,ieff}_axes_report, so co-locating them would let the planner
    # pick a lane by tiebreak. Which lane runs is a deliberate choice made by
    # the caller, not a planner accident -- see ECSPR_DOMAINS below.
    "ecspr",
    "ecsprDirected",
    "ecsprUndirected",
    "ecsprNetA",
    "ecsprNetB",
    # The benchmark lane PRODUCES no base graphs -- it CONSUMES
    # `ecspr::benchmark_base_graphs`, a staged item of the frozen benchmark
    # tree. (An earlier reading had it building them from X; X withholds the
    # atom mapping on purpose, and the atom-transit count is both the edge
    # weight and, via the w<=0 drop, the connectivity, so they are not
    # derivable from X. See transforms/build/benchmark/NOTES.md.)
    #
    # It is still mutually exclusive with netA/netB, for the opposite reason:
    # selecting it must DROP them, because on the benchmark the annotation
    # chain that feeds `ecspr::base_graphs` is upstream of the freeze and must
    # not be planned at all. A score has to measure the solver, not the lanes.
    "ecsprBenchmark",
]

# The ECSPr domains that must NOT both be offered to the planner in one run,
# keyed by the selector that chooses between them. Consumed by pipeline.py.
ECSPR_SOLVE_DOMAINS = {"directed": "ecsprDirected", "undirected": "ecsprUndirected"}
ECSPR_NETWORK_DOMAINS = {"A": "ecsprNetA", "B": "ecsprNetB", "benchmark": "ecsprBenchmark"}


def domains_for(solve: str | None = None, network: str | None = None) -> list[str]:
    """DOMAINS with the unselected ECSPr lanes removed.

    Offering both solve lanes at once makes ``ecspr::reff_axes_report`` have two
    producers, and the planner would resolve that by tiebreak rather than by
    intent. The same holds for the two base-graph lanes. Callers that do not ask
    for ECSPr get the full list minus every ECSPr domain.
    """
    drop: set[str] = set()
    if solve is None:
        drop |= set(ECSPR_SOLVE_DOMAINS.values()) | {"ecspr"}
    else:
        if solve not in ECSPR_SOLVE_DOMAINS:
            raise ValueError(f"unknown solve lane [{solve}], expected one of {sorted(ECSPR_SOLVE_DOMAINS)}")
        drop |= {d for k, d in ECSPR_SOLVE_DOMAINS.items() if k != solve}
    if network is None:
        drop |= set(ECSPR_NETWORK_DOMAINS.values())
    else:
        if network not in ECSPR_NETWORK_DOMAINS:
            raise ValueError(f"unknown network [{network}], expected one of {sorted(ECSPR_NETWORK_DOMAINS)}")
        drop |= {d for k, d in ECSPR_NETWORK_DOMAINS.items() if k != network}
    return [d for d in DOMAINS if d not in drop]

_MODULE = Path(__file__).resolve().parent


def _looks_like_library(root: Path) -> bool:
    return all((root / sub).is_dir() for sub in ("data_types", "resources", "transforms"))


def resolve_library_root() -> Path:
    override = os.environ.get("FABFOS_LIBRARY")
    if override:
        root = Path(override).expanduser().resolve()
        if not _looks_like_library(root):
            raise FileNotFoundError(
                f"FABFOS_LIBRARY=[{root}] is not a metasmith library "
                f"(missing data_types/ resources/ transforms/)"
            )
        return root

    bundled = _MODULE / "_library"
    if _looks_like_library(bundled):
        return bundled

    # src/fabfos/library.py -> _MODULE == src/fabfos, _MODULE.parent == src/
    dev_sibling = _MODULE.parent / "metasmith_libraries"
    if _looks_like_library(dev_sibling):
        return dev_sibling

    raise FileNotFoundError(
        "could not locate the FabFos metasmith library. Set FABFOS_LIBRARY, "
        "install the package with a bundled library, or run from a source "
        "checkout with the metasmith_libraries submodule present."
    )

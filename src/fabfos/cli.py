"""FabFos command line — a thin front end over the metasmith planner.

``fabfos`` builds a metasmith workflow that resolves fosmid inserts from
pooled reads and runs it through the chosen container runtime (apptainer by
default, which is what the cluster provides). Use ``--plan-only`` to just
resolve and print the DAG without executing.
"""
import argparse
import multiprocessing
import sys
from pathlib import Path

from metasmith.python_api import ContainerRuntime

from . import __version__, NAME, SHORT_SUMMARY
from .pipeline import FabFosInputs, generate_workflow, run_pipeline


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog=NAME, description=SHORT_SUMMARY)
    io = p.add_argument_group("inputs")
    io.add_argument("-r", "--reads", metavar="FASTQ", required=True,
                    help="forward (paired) or interleaved or single-end reads, fastq[.gz]")
    io.add_argument("-2", "--reverse", metavar="FASTQ", default=None,
                    help="reverse reads for paired-end input")
    io.add_argument("-i", "--interleaved", action="store_true", default=False,
                    help="--reads holds interleaved paired-end reads")
    io.add_argument("-o", "--output", metavar="PATH", required=True,
                    help="output directory")

    fos = p.add_argument_group("fosmid pool")
    fos.add_argument("-b", "--background", metavar="FASTA", default=None,
                     help="host background genome to filter out")
    fos.add_argument("--vector", metavar="FASTA", default=None,
                     help="vector backbone fasta; enables pool-size estimation")
    fos.add_argument("--endf", metavar="FASTA", default=None,
                     help="forward-junction end sequences fasta")
    fos.add_argument("--endr", metavar="FASTA", default=None,
                     help="reverse-junction end sequences fasta")
    fos.add_argument("--ends-facing", action="store_true", default=False,
                     help="end sequences face inward across both junctions")

    ec = p.add_argument_group("ecspr (metabolic graph prerequisite)")
    ec.add_argument("--ecspr", action="store_true", default=False,
                    help="also build per-fosmid bipartite metabolic graphs (stops before the axis/conductance step)")
    ec.add_argument("--base-graphs", metavar="DIR", default=None,
                    help="reference dir of base_{C,N,S,P}.pkl host graphs")
    ec.add_argument("--element-bipartite", metavar="DIR", default=None,
                    help="reference dir of mnx_bipartite_{C,N,S,P}.pkl universe graphs")
    ec.add_argument("--reaction-db", metavar="DIR", default=None,
                    help="reference dir with reactions.dmnd + bridge.tsv")

    run = p.add_argument_group("execution")
    run.add_argument("--runtime", choices=[r.value for r in ContainerRuntime],
                     default=ContainerRuntime.APPTAINER.value,
                     help="container runtime (default: apptainer)")
    run.add_argument("--solve-lane", choices=["directed", "undirected"], default=None,
                     help="which ECSPr solve lane to offer the planner; both produce "
                          "the same axes report, so leaving this unset offers neither "
                          "rather than letting the planner pick by tiebreak")
    run.add_argument("--network-lane", choices=["A", "B", "benchmark"], default=None,
                     help="which ECSPr base-graph lane to offer the planner")
    run.add_argument("-t", "--threads", type=int, default=multiprocessing.cpu_count(),
                     help="max threads per step")
    run.add_argument("--plan-only", action="store_true", default=False,
                     help="resolve and print the workflow DAG without executing")
    run.add_argument("--dag", metavar="PATH", default=None,
                     help="render the resolved DAG to PATH.svg. The render comes "
                          "from the SAME planner call the run makes, so it documents "
                          "what would actually execute, not a hand-drawn idea of it")
    run.add_argument("--provision-only", action="store_true", default=False,
                     help="create the per-tool mamba envs from the library *.env.yml specs, then exit")
    run.add_argument("--no-provision", action="store_true", default=False,
                     help="accepted and ignored: no runtime on the pinned engine needs "
                          "per-tool mamba envs, so nothing is auto-created. Use "
                          "--provision-only to stand them up explicitly")
    p.add_argument("-v", "--version", action="version", version=f"{NAME} {__version__}")

    # The METHOD version, distinct from the package version above. The package
    # is the CLI; the method is the composition (canon + library commit +
    # container digests + type contract + data-library index) that decides what
    # a number out of this pipeline means. A CLI bugfix is not a new method.
    meth = p.add_argument_group("method version")
    meth.add_argument("--method-version", action="store_true", default=False,
                      help="print the method id (version+hash) and exit")
    meth.add_argument("--describe-method", action="store_true", default=False,
                      help="print the full hashed method document and exit; "
                           "diff two of these to see WHICH component moved")
    meth.add_argument("--require-method", metavar="ID", default=None,
                      help="fail unless the live method matches ID "
                           "(full '0.3.0+abc1234' or bare '0.3.0')")
    return p


def _inputs_from_args(a: argparse.Namespace) -> FabFosInputs:
    return FabFosInputs(
        reads=Path(a.reads).resolve(),
        output=Path(a.output).resolve(),
        reverse=Path(a.reverse).resolve() if a.reverse else None,
        interleaved=a.interleaved,
        background=Path(a.background).resolve() if a.background else None,
        vector=Path(a.vector).resolve() if a.vector else None,
        end_forward=Path(a.endf).resolve() if a.endf else None,
        end_reverse=Path(a.endr).resolve() if a.endr else None,
        ends_facing=a.ends_facing,
        runtime=ContainerRuntime(a.runtime),
        threads=a.threads,
        solve_lane=a.solve_lane,
        network_lane=a.network_lane,
        ecspr=a.ecspr,
        base_graphs=Path(a.base_graphs).resolve() if a.base_graphs else None,
        element_bipartite=Path(a.element_bipartite).resolve() if a.element_bipartite else None,
        reaction_db=Path(a.reaction_db).resolve() if a.reaction_db else None,
    )


def _method_query(argv: list[str] | None) -> str | None:
    """Answer --method-version / --describe-method BEFORE the main parser.

    The main parser requires --reads and --output. Asking what method this is
    is a question about the installation, not about a run, so it must not
    require a runnable set of reads to answer.
    """
    import sys as _sys

    args = list(_sys.argv[1:] if argv is None else argv)
    for flag, key in (("--method-version", "version"), ("--describe-method", "describe")):
        if flag in args:
            return key
    return None


def _render_dag(task, base: Path) -> Path:
    """Render the resolved plan, refusing to draw an incomplete one.

    A renderer will happily draw a disconnected graph for a plan that never
    resolved, and the picture looks authoritative either way. If the plan is
    not ok, that is the thing worth reporting -- not a diagram of it.
    """
    if not getattr(task, "ok", False):
        raise RuntimeError(
            "refusing to render a DAG for a plan that did not resolve -- "
            "the drawing would look complete regardless. Fix the plan first."
        )
    base = base.resolve()
    base.parent.mkdir(parents=True, exist_ok=True)
    stem = base.with_suffix("") if base.suffix else base
    task.plan.RenderDAG(stem)
    out = stem.with_suffix(".svg")
    if not out.exists():
        raise RuntimeError(f"RenderDAG reported success but {out} is absent")
    return out


def main(argv: list[str] | None = None) -> int:
    query = _method_query(argv)
    if query is not None:
        from .method import describe_method

        desc = describe_method()
        if query == "version":
            print(desc.method_id)
            if desc.unresolved_containers:
                print(
                    "  NOT STAMPABLE -- unresolved containers: "
                    f"{', '.join(sorted(desc.unresolved_containers))}",
                    file=sys.stderr,
                )
                return 1
            return 0
        import yaml

        print(yaml.safe_dump(desc.to_dict(), sort_keys=False, default_flow_style=False))
        return 0

    args = _build_parser().parse_args(argv)

    if args.require_method is not None:
        from .method import MethodError, check_required

        try:
            check_required(args.require_method)
        except MethodError as e:
            print(f"fabfos: {e}", file=sys.stderr)
            return 1

    inp = _inputs_from_args(args)

    if args.provision_only:
        from .provision import provision_tool_environments
        from .library import resolve_library_root
        rep = provision_tool_environments(resolve_library_root())
        print(f"provision: created={rep.created} skipped={rep.skipped} failed={[n for n,_ in rep.failed]}")
        return 1 if rep.failed else 0

    if args.plan_only:
        inp.output.mkdir(parents=True, exist_ok=True)
        _agent, task = generate_workflow(inp, inp.output / "_fabfos")
        if not task.ok:
            print("workflow generation FAILED; planner hints:", file=sys.stderr)
            print(task.plan.RenderHints() if hasattr(task.plan, "RenderHints") else task, file=sys.stderr)
            return 1
        print(f"resolved workflow: {len(task.plan.steps)} steps")
        for i, step in enumerate(task.plan.steps):
            name = getattr(getattr(step, "transform", None), "name", None) or f"step{i}"
            print(f"  [{i}] {name}")
        if args.dag:
            print(f"DAG -> {_render_dag(task, Path(args.dag))}")
        return 0

    run_pipeline(inp, provision=not args.no_provision)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

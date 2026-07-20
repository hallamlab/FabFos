"""Gate and freeze benchmark v3.

Four gates, then a per-file sha256 manifest. A gate failure means the benchmark
is NOT frozen and the script exits nonzero -- the freeze is the reward for the
gates passing, never a separate step someone can run anyway.

WHY THIS IS A FORK AND NOT A FLAG ON 50_audit.py
------------------------------------------------
v1's `Y_BUILDERS` is a hardcoded list of v1's filenames. Pointed at v3 it would
scan four files that are not v3's builders, find nothing, and report PASS --
an audit that passes precisely because it audited nothing. The list is the thing
that has to change, so the file is the thing that has to fork.

GATE 4 IS NEW
-------------
v1 had three gates. v3 adds ANSWER-KEY PRECEDENCE: the key must be frozen and
hashed BEFORE any solve output exists, because a key written after the fact is
not an answer key. The gate refuses if a scored result is already sitting in the
tree, which is the one ordering error no amount of care in the builders can
catch.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _v3bench import BUILD_OUT, ROOT, X_DIR, Y_DIR, sha256_of  # noqa: E402

HERE = Path(__file__).resolve().parent

# v3's Y builders -- their INPUT paths must never touch ECSPr output.
Y_BUILDERS = ["20_build_panel_v3.py", "21_build_Y_v3.py",
              "10_build_targets.py", "_v3bench.py"]

# Tokens that name ECSPr OUTPUT (not input ground truth). `base_` is absent on
# purpose: base graphs are network CONSTRUCTION, which the panel's liveness gate
# legitimately reads. That line is deliberate -- do not "fix" it by adding them.
FORBIDDEN = [r"\bout/", r"panel_\w*\.parquet", r"analysis_report", r"_matrix_scored",
             r"potency_\w+\.tsv", r"specificity_\w+\.tsv", r"\breff_", r"ieff",
             r"sig_mix", r"sig_emp", r"sig_negbin",
             # The SOLVE output is named `observations.tsv`. The observation SET
             # is `gof_observations.tsv` / `lof_observations.tsv` and is
             # legitimate input, so the lookbehind is load-bearing: without it
             # this token fails all four builders for reading their own source.
             r"(?<![\w_])observations\.tsv"]


def gate_y_independence() -> list[str]:
    """No Y builder may READ an ECSPr-output path.

    A forbidden token counts only inside an actual file-read call on the same
    line, which is what keeps the gate from flagging its own contract prose --
    the docstrings name these tokens precisely so a reader knows the rule.

    KNOWN LIMIT, stated rather than hidden: a path assembled on one line and
    read on the next slips through. The gate is a tripwire against drift, not a
    proof. `observations.tsv` is in the list because that is the name the SOLVE
    output takes; the observation SET is `gof_observations.tsv` /
    `lof_observations.tsv` and is legitimate input.
    """
    reads = re.compile(r"read_csv|read_parquet|read_pickle|\bopen\s*\(|\.load|Path\s*\(")
    fails = []
    for name in Y_BUILDERS:
        path = HERE / name
        if not path.exists():
            fails.append(f"{name}: LISTED AS A Y BUILDER BUT ABSENT -- the audit "
                         f"would silently scan nothing")
            continue
        for i, ctx in enumerate(path.read_text().splitlines(), start=1):
            if not reads.search(ctx):
                continue
            for pat in FORBIDDEN:
                m = re.search(pat, ctx)
                if m:
                    fails.append(f"{name}:{i}: reads ECSPr output "
                                 f"`{m.group()}` -> {ctx.strip()[:70]}")
    return fails


def gate_x_sufficiency() -> list[str]:
    fails = []
    mets = pd.read_csv(X_DIR / "metabolites.tsv", sep="\t", dtype=str).fillna("")
    cond = pd.read_csv(Y_DIR / "conditions.tsv", sep="\t", dtype=str).fillna("")
    panel = pd.read_csv(Y_DIR / "panel_edges.tsv", sep="\t", dtype=str).fillna("")

    has_struct = dict(zip(mets.mnxm, mets.has_structure == "True"))
    y_mets = set(panel.product_mnxm) | set(panel.anchor_mnxm)
    y_mets |= {m for t in cond.target_mnxm for m in str(t).split(",") if m}

    missing = [m for m in y_mets if m and m not in has_struct]
    if missing:
        fails.append(f"{len(missing)} Y metabolites absent from X/metabolites.tsv: "
                     f"{missing[:5]}")

    structless = [m for m in y_mets if m in has_struct and not has_struct[m]]
    if structless:
        print(f"  note: {len(structless)} Y metabolites have no structure "
              f"(irreducible polymers/generics): {structless[:4]}")
    return fails


def gate_lof_uniformity() -> list[str]:
    fails = []
    cond = pd.read_csv(Y_DIR / "conditions.tsv", sep="\t", dtype=str).fillna("")
    exp = pd.read_csv(Y_DIR / "expectations.tsv", sep="\t", dtype=str).fillna("")

    lof = cond[cond.arm == "lof"]
    if set(lof.expected_dir) != {"-"}:
        fails.append(f"LOF arm carries non-'-' directions: "
                     f"{sorted(set(lof.expected_dir))}")
    if len(lof) != 166:
        # 166 is the observation count and every one must survive. This caught a
        # real defect: v3 writes the gene token as `argA:del` while pheno_edges
        # keys on `argA`, so the join matched 0 of 166 and about half then
        # resolved through the free-text fallback by luck -- which looks exactly
        # like ordinary coverage loss.
        fails.append(f"expected 166 LOF conditions, found {len(lof)} -- every "
                     f"deletion observation must carry a target")

    bad = exp[(exp.role == "target") & (exp["dir"] == "0")]
    if len(bad):
        fails.append(f"{len(bad)} expectation rows are target/dir=0 -- that code "
                     f"path must be empty, not merely unused")
    return fails


def gate_key_precedence() -> list[str]:
    """The key must be frozen BEFORE any solve output exists.

    A key authored after seeing results is not an answer key, and no amount of
    care inside the builders can detect that ordering error. This gate can.
    """
    fails = []
    for stray in sorted(ROOT.rglob("observations.tsv")) + sorted(ROOT.rglob("*_scored*")):
        fails.append(f"solve/scored output already present at {stray.relative_to(ROOT)} "
                     f"-- the answer key cannot be frozen after a result exists")
    return fails


def manifest() -> pd.DataFrame:
    rows = []
    for d, tag in ((X_DIR, "X"), (Y_DIR, "Y")):
        for f in sorted(d.glob("*.tsv")):
            df = pd.read_csv(f, sep="\t", dtype=str)
            rows.append(dict(artifact=tag, file=f.name, rows=len(df),
                             cols=df.shape[1], sha256=sha256_of(f)))
        for f in sorted(d.glob("*.json")):
            rows.append(dict(artifact=tag, file=f.name, rows="", cols="",
                             sha256=sha256_of(f)))
    return pd.DataFrame(rows)


def main() -> int:
    gates = [
        ("Y independence from ECSPr output", gate_y_independence,
         "no Y builder reads ECSPr output"),
        ("X sufficiency for Y", gate_x_sufficiency,
         "every Y metabolite resolves in X"),
        ("LOF uniformity", gate_lof_uniformity,
         "all 166 LOF dir=-, no target/dir=0 cell"),
        ("answer-key precedence", gate_key_precedence,
         "no solve output exists yet -- the key freezes first"),
    ]

    all_fails, results = [], {}
    for n, (title, fn, ok_msg) in enumerate(gates, start=1):
        print(f"{'' if n == 1 else chr(10)}GATE {n} -- {title}\n" + "=" * 60)
        f = fn()
        results[title] = not f
        print(f"  PASS -- {ok_msg}" if not f else "\n".join("  FAIL " + x for x in f))
        all_fails += f

    man = manifest()
    man.to_csv(BUILD_OUT / "MANIFEST.tsv", sep="\t", index=False)
    print(f"\nMANIFEST -- {len(man)} files\n" + "=" * 60)
    print(man.to_string(index=False))

    (BUILD_OUT / "_audit.json").write_text(json.dumps(dict(
        version="v3", gates=results, failures=all_fails, n_files=len(man),
    ), indent=2) + "\n")

    if all_fails:
        print(f"\n{len(all_fails)} GATE FAILURE(S) -- benchmark is NOT frozen")
        return 1
    print(f"\nALL GATES PASS -- benchmark v3 frozen at {ROOT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

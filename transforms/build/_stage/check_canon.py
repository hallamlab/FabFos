#!/usr/bin/env python3
"""Gate: canon resolves through the library, and fails loudly when it cannot.

The properties worth protecting are not "the paths are right" -- they are:

  1. Importing canon for a SCALAR must not require a library or metasmith.
     The paper worktrees do exactly this, and breaking it breaks them.
  2. Every data path must come from the manifest. A module-level assignment
     shadows __getattr__, so a reintroduced absolute path would keep every
     other test green while silently un-doing the migration. Check that the
     resolved path is under the library root.
  3. A missing key must raise CanonError naming the symbol -- never fall back.
     A fallback makes the migration untestable.
  4. canon's transcribed hashes must agree with the provenance records.
"""
from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
failures: list[str] = []


def check(label: str, ok: bool, detail: str = "") -> None:
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}" + (f" -- {detail}" if detail else ""))
    if not ok:
        failures.append(label)


def main() -> int:
    sys.path.insert(0, str(REPO / "src"))

    print("1. scalars import with no library present")
    env = {**os.environ, "FABFOS_REF": "/nonexistent", "PYTHONPATH": str(REPO / "src")}
    probe = "from fabfos import canon; print(canon.KMAX, canon.AXES_N, canon.CANONICAL_ORIENTATION)"
    proc = subprocess.run([sys.executable, "-c", probe], capture_output=True, text=True, env=env)
    check("canon.KMAX etc. with FABFOS_REF=/nonexistent", proc.returncode == 0,
          proc.stderr.strip().splitlines()[-1] if proc.returncode else proc.stdout.strip())

    from fabfos import canon

    print("\n2. every data path resolves UNDER the library root")
    root = canon.library_root().resolve()
    for name in sorted(canon._PATHS):
        try:
            p = getattr(canon, name).resolve()
        except canon.CanonError as e:
            check(name, False, str(e)[:90])
            continue
        under = root in p.parents or p == root
        check(f"{name} under <lib> and present", under and p.exists(),
              "" if under and p.exists() else f"{p}")

    print("\n3. a missing key raises CanonError naming the symbol")
    canon._PATHS["_GATE_PROBE"] = "derived/definitely-not-here"
    try:
        canon._GATE_PROBE
        check("missing key raises", False, "it did not raise -- a fallback has crept in")
    except canon.CanonError as e:
        check("missing key raises CanonError naming the symbol", "_GATE_PROBE" in str(e))
    finally:
        del canon._PATHS["_GATE_PROBE"]

    print("\n4. canon's transcribed hashes agree with the provenance records")
    import yaml
    pairs = [
        ("REFERENCE_DIRECTION_SHA256", "mnxref.direction"),
        ("REFERENCE_REAC_PROP_SHA256", "metanetx.reac_prop"),
    ]
    for symbol, item_id in pairs:
        if not hasattr(canon, symbol):
            print(f"  [SKIP] canon.{symbol} not defined")
            continue
        rec_path = REPO / "provenance" / "data" / f"{item_id}.yml"
        if not rec_path.exists():
            check(symbol, False, f"no record {item_id}.yml")
            continue
        rec = yaml.safe_load(rec_path.open())
        check(f"canon.{symbol} == {item_id}.sha256",
              getattr(canon, symbol) == rec.get("sha256"))

    print("\n5. the frozen null files canon names are all in the library")
    d = canon.REFERENCE_NULL_DIR
    missing = [f for f in canon.FROZEN_NULL_FILES if not (d / f).exists()]
    check(f"{len(canon.FROZEN_NULL_FILES)} frozen null files present", not missing,
          f"missing: {missing}")

    if failures:
        print(f"\n{len(failures)} FAILURE(S): {', '.join(failures)}")
        return 1
    print("\ncanon gate green")
    return 0


if __name__ == "__main__":
    sys.exit(main())

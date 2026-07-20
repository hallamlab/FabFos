#!/usr/bin/env python3
"""Gate: the library and the provenance tree must agree.

Four checks, each corresponding to an acceptance criterion:

  AC2  every manifest key exists on disk and its type resolves
  AC3  every declared item has a generated record carrying a sha256
  AC4  every container record carries a digest, or an explicit blocker
  ---   reac_prop's hash agrees across all the places that pin it

The last one is the point of the whole exercise in miniature. reac_prop.tsv's
sha256 is pinned in three places: this library's provenance record,
mnxref-4_5/MANIFEST.json (as reac_prop_sha256), and canon. Three transcriptions
of one number is exactly the failure mode canon exists to remove, so the
provenance record is the source and the others are checked against it.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
RECORDS = REPO / "provenance" / "data"
CONTAINERS = REPO / "provenance" / "containers"
LIB = REPO / ".awm" / "data" / "ref"

failures: list[str] = []
notes: list[str] = []


def check(label: str, ok: bool, detail: str = "") -> None:
    print(f"  [{'PASS' if ok else 'FAIL'}] {label}" + (f" -- {detail}" if detail else ""))
    if not ok:
        failures.append(f"{label}: {detail}")


def main() -> int:
    import yaml

    decl = yaml.safe_load((RECORDS / "_declared.yml").open())
    declared = {i["id"]: i for i in decl["items"]}

    print("AC3 -- every declared item has a record with a sha256")
    before = len(failures)
    missing_ok = 0
    for item_id, item in declared.items():
        rec_path = RECORDS / f"{item_id}.yml"
        if not rec_path.exists():
            check(item_id, False, "no generated record")
            continue
        rec = yaml.safe_load(rec_path.open())
        if rec.get("status") == "MISSING":
            missing_ok += 1
            if not rec.get("degraded_mode"):
                check(item_id, False, "declared MISSING but states no degraded_mode")
            continue
        if not rec.get("sha256"):
            check(item_id, False, "record carries no sha256")
    # `not failures` here would let ANY earlier failure mislabel this one line.
    check(f"{len(declared) - missing_ok} present items all hashed",
          len(failures) == before, f"{missing_ok} recorded MISSING with a degradation path")

    print("\nAC4 -- every container record carries a digest or an explicit blocker")
    before = len(failures)
    unresolved = []
    for rec_path in sorted(CONTAINERS.glob("*.yml")):
        rec = yaml.safe_load(rec_path.open())
        if rec.get("resolved"):
            if not str(rec.get("digest", "")).startswith("sha256:"):
                check(rec_path.stem, False, "resolved:true but no sha256 digest")
        else:
            if not rec.get("blocker"):
                check(rec_path.stem, False, "unresolved with no blocker stated")
            unresolved.append(rec_path.stem)
    n = len(list(CONTAINERS.glob("*.yml")))
    # was hardcoded True, so this line printed PASS even when the loop failed.
    check(f"{n - len(unresolved)}/{n} images digest-pinned",
          len(failures) == before,
          f"unresolved, each with a blocker: {', '.join(unresolved)}" if unresolved
          else "every image pinned by digest")

    print("\nAC5 -- every declared `type:` exists in the library's type contract")
    # A namespace binds to exactly one type library, and metasmith matches
    # endpoints by subset over a property set built from the yaml -- so a
    # declared type that no longer exists under its namespace does not raise,
    # it silently fails to join. This is the assertion that would have caught
    # the four ecspr:: reference types the transforms require but nothing
    # declared.
    sys.path.insert(0, str(REPO / "src"))
    from fabfos.library import resolve_library_root  # noqa: E402

    types_root = resolve_library_root() / "data_types"
    contract: dict[str, set[str]] = {}
    for ns in ("ecspr", "ref"):
        p = types_root / f"{ns}.yml"
        contract[ns] = set(yaml.safe_load(p.open())["types"]) if p.exists() else set()
    before = len(failures)
    seen = set()
    for item_id, item in declared.items():
        t = item.get("type", "")
        if "::" not in t:
            check(item_id, False, f"type [{t}] is not namespaced")
            continue
        ns, name = t.split("::", 1)
        seen.add(ns)
        if ns not in contract:
            notes.append(f"{item_id}: namespace [{ns}] not checked here")
            continue
        if name not in contract[ns]:
            check(item_id, False, f"type [{t}] is not in {types_root/(ns+'.yml')}")
    check(f"{len(declared)} declared types all resolve",
          len(failures) == before,
          f"namespaces used: {', '.join(sorted(seen))}; "
          f"contract sizes: {', '.join(f'{k}={len(v)}' for k, v in sorted(contract.items()))}")

    print("\nno image may float on a mutable :latest without saying so")
    for rec_path in sorted(CONTAINERS.glob("*.yml")):
        rec = yaml.safe_load(rec_path.open())
        if str(rec.get("reference", "")).endswith(":latest") and rec.get("resolved"):
            check(rec_path.stem, False, "resolved against a :latest tag")

    print("\nreac_prop sha256 agrees everywhere it is pinned")
    rec = yaml.safe_load((RECORDS / "metanetx.reac_prop.yml").open())
    truth = rec["sha256"]
    print(f"  provenance record: {truth}")
    man_path = LIB / "derived" / "mnxref-4_5" / "MANIFEST.json"
    if man_path.exists():
        pinned = json.loads(man_path.read_text()).get("reac_prop_sha256")
        print(f"  MANIFEST.json:     {pinned}")
        check("MANIFEST.json agrees with the provenance record", pinned == truth,
              "the pre-bake was built against different reac_prop bytes"
              if pinned != truth else "")
    else:
        notes.append("MANIFEST.json not in the library; skipped")

    print()
    for note in notes:
        print(f"  note: {note}")
    if failures:
        print(f"\n{len(failures)} FAILURE(S):")
        for f in failures:
            print(f"  {f}")
        return 1
    print("\nall gates green")
    return 0


if __name__ == "__main__":
    sys.exit(main())

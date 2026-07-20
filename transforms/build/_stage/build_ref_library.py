#!/usr/bin/env python3
"""Build `.awm/data/ref/` -- the metasmith DataInstanceLibrary -- from the
declaration table at `provenance/data/_declared.yml`.

Three phases, each independently re-runnable:

    --stage   resolve every declared source, verify it exists, and write the
              manifest with ABSOLUTE keys. Zero bytes move. This is where you
              find out the declaration is wrong, at no cost.
    --hash    sha256 every item and emit one generated record per item under
              provenance/data/. Pure read; slow; safe to background.
    --place   hardlink each item into its tiered destination, rewrite the
              manifest keys to relative, and verify by inode identity.

Why hardlinks. Every source root is on device 2128, so a link costs zero bytes
and shares the inode. The three live canon versions keep reading the same data
through their existing absolute paths -- there is no cutover moment and no
window in which two copies can diverge. Rollback is `rm -rf` on the new tree;
the originals are the same inodes and are untouched.

Why not symlinks: metasmith's ExecuteTransfers() defaults to
resolve_symlinks=False and only passes rsync -L when true, so a symlinked
library body arrives on a remote as dangling links. Silently. That is precisely
the failure this work exists to prevent.

Why not LocalizeContents(): it flattens every item to path.name at the library
root and disambiguates collisions as _1/_2. It would destroy the tiered tree
and silently merge the two ieff_axes_report.tsv.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
DECLARED = REPO / "provenance" / "data" / "_declared.yml"
RECORDS = REPO / "provenance" / "data"

# ---------------------------------------------------------------------------
# declaration loading
# ---------------------------------------------------------------------------


def load_declaration() -> dict:
    try:
        import yaml
    except ImportError:
        sys.exit("PyYAML is required: mamba run -n p312 python ...")
    with DECLARED.open() as fh:
        return yaml.safe_load(fh)


def expand_braces(pattern: str) -> list[str]:
    """`x_N{14,28}.tsv` -> [`x_N14.tsv`, `x_N28.tsv`]. One group, which is all
    the declaration uses; a second group would silently be ignored, so refuse."""
    matches = re.findall(r"\{([^}]*)\}", pattern)
    if not matches:
        return [pattern]
    if len(matches) > 1:
        raise ValueError(f"only one brace group is supported: {pattern}")
    return [pattern.replace("{" + matches[0] + "}", alt) for alt in matches[0].split(",")]


def resolve_members(item: dict, src: Path) -> list[tuple[Path, str]]:
    """Return [(absolute source file, path relative to the item's dest)]."""
    if src.is_file():
        return [(src, "")]

    selects = item.get("select")
    if selects:
        out = []
        for raw in selects:
            for pattern in expand_braces(raw):
                hits = sorted(src.glob(pattern))
                if not hits:
                    raise FileNotFoundError(f"{item['id']}: select matched nothing: {pattern}")
                out += [(h, h.relative_to(src).as_posix()) for h in hits if h.is_file()]
        return out

    excludes = item.get("exclude", [])
    out = []
    for path in sorted(src.rglob("*")):
        if not path.is_file():
            continue
        rel = path.relative_to(src)
        if any(rel.match(pat) or rel.as_posix().startswith(pat.rstrip("*/")) for pat in excludes):
            continue
        out += [(path, rel.as_posix())]
    return out


def build_plan(decl: dict) -> tuple[list[dict], list[dict]]:
    """(resolvable items, declared-but-missing items)."""
    root = Path(decl["defaults"]["data_root"])
    plan, missing = [], []

    for item in decl["items"]:
        if item.get("status") == "MISSING" or item.get("src") is None:
            missing.append(item)
            continue
        raw_src = str(item["src"])
        # `repo/...` resolves against the REPO, not the data root. Needed the
        # moment the repo itself produces a declared item (the benchmark's
        # target-resolution table is the first). An absolute path would pin the
        # declaration to one checkout, so it would resolve in the worktree that
        # wrote it and be MISSING everywhere else -- including after a merge.
        if raw_src.startswith("repo/"):
            src = REPO / raw_src[len("repo/"):]
        else:
            src = Path(raw_src)
            if not src.is_absolute():
                src = root / src
        if not src.exists():
            missing.append({**item, "_reason": f"declared source does not exist: {src}"})
            continue
        try:
            members = resolve_members(item, src)
        except FileNotFoundError as e:
            missing.append({**item, "_reason": str(e)})
            continue
        plan.append({**item, "_src": src, "_members": members})
    return plan, missing


# ---------------------------------------------------------------------------
# phases
# ---------------------------------------------------------------------------


def phase_stage(decl: dict, lib: Path) -> int:
    plan, missing = build_plan(decl)
    meta = lib / "_metadata"
    meta.mkdir(parents=True, exist_ok=True)

    manifest, n_files, n_bytes = {}, 0, 0
    for item in plan:
        for src, rel in item["_members"]:
            key = str(src)  # ABSOLUTE -- legal, and in production use elsewhere
            manifest[key] = item["type"]
            n_files += 1
            n_bytes += src.stat().st_size
        print(f"  {item['id']:28s} {len(item['_members']):4d} file(s)  {item['tier']}")

    (meta / "index.yml").write_text(
        "# staged with ABSOLUTE keys -- no bytes have moved yet\n"
        + "\n".join(f'? {k}\n: {v}' for k, v in sorted(manifest.items()))
        + "\n"
    )
    (meta / "stage.json").write_text(json.dumps(
        {"phase": "stage", "items": len(plan), "files": n_files, "bytes": n_bytes,
         "missing": [m["id"] for m in missing]}, indent=2) + "\n")

    print(f"\nstaged {len(plan)} items / {n_files} files / {n_bytes/1e9:.1f} GB (0 bytes moved)")
    if missing:
        print(f"\n{len(missing)} declared item(s) NOT resolvable:")
        for m in missing:
            print(f"  {m['id']:28s} {m.get('_reason', m.get('note', 'declared MISSING'))[:120]}")
    return 0


def _same_bytes(a: Path, b: Path) -> bool:
    return a.stat().st_size == b.stat().st_size and sha256_of(a) == sha256_of(b)


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def phase_hash(decl: dict, lib: Path) -> int:
    plan, missing = build_plan(decl)
    RECORDS.mkdir(parents=True, exist_ok=True)

    for item in plan:
        total, digests = 0, {}
        for src, rel in item["_members"]:
            total += src.stat().st_size
            digests[rel or src.name] = sha256_of(src)
            print(f"    {rel or src.name}", flush=True)

        # a single-file item gets its own hash; a bundle gets a hash OF THE
        # HASHES, so the item has one identity without pretending a directory
        # has a checksum.
        if len(digests) == 1:
            item_hash = next(iter(digests.values()))
        else:
            joined = "".join(f"{k}:{v}\n" for k, v in sorted(digests.items()))
            item_hash = "tree-" + hashlib.sha256(joined.encode()).hexdigest()

        rec = {k: v for k, v in item.items() if not k.startswith("_")}
        rec.update({"bytes": total, "sha256": item_hash, "n_files": len(digests),
                    "files": digests})
        _write_record(item["id"], rec)
        print(f"  {item['id']:28s} {item_hash[:24]}  {total/1e6:9.1f} MB", flush=True)

    for m in missing:
        rec = {k: v for k, v in m.items() if not k.startswith("_")}
        rec.update({"status": "MISSING", "sha256": None,
                    "reason": m.get("_reason", "declared MISSING")})
        _write_record(m["id"], rec)

    print(f"\nhashed {len(plan)} items; {len(missing)} recorded MISSING")
    return 0


def _write_record(item_id: str, rec: dict) -> None:
    import yaml
    path = RECORDS / f"{item_id}.yml"
    path.write_text(
        "# GENERATED by transforms/build/_stage/build_ref_library.py --hash\n"
        "# Source of truth is provenance/data/_declared.yml. Do not hand-edit.\n"
        + yaml.safe_dump(rec, sort_keys=False, default_flow_style=False, width=88)
    )


def phase_place(decl: dict, lib: Path) -> int:
    plan, _ = build_plan(decl)
    manifest, linked, reused, copied, failed = {}, 0, 0, [], []

    for item in plan:
        for src, rel in item["_members"]:
            dest = lib / item["dest"] / rel if rel else lib / item["dest"]
            dest.parent.mkdir(parents=True, exist_ok=True)

            if dest.exists():
                if dest.stat().st_ino == src.stat().st_ino:
                    reused += 1
                    manifest[dest.relative_to(lib).as_posix()] = item["type"]
                    continue
                if _same_bytes(src, dest):
                    reused += 1
                    manifest[dest.relative_to(lib).as_posix()] = item["type"]
                    continue
                failed.append((str(dest), "exists, different inode AND different bytes"))
                continue

            try:
                os.link(src, dest)
                linked += 1
            except OSError as e:
                # protected_hardlinks=1 forbids linking a file you neither own
                # nor can write. Some reference files are root-owned container
                # outputs. Copy those and verify by content instead -- the
                # inode-identity property is unavailable, so say so out loud.
                if e.errno != 1:
                    failed.append((str(dest), f"link failed: {e}"))
                    continue
                import shutil
                shutil.copy2(src, dest)
                if not _same_bytes(src, dest):
                    failed.append((str(dest), "copy fallback: content mismatch"))
                    continue
                copied.append((item["id"], dest.relative_to(lib).as_posix()))
                manifest[dest.relative_to(lib).as_posix()] = item["type"]
                continue

            # inode identity is O(1) and strictly stronger than a hash match
            if dest.stat().st_ino != src.stat().st_ino:
                failed.append((str(dest), "post-link inode mismatch"))
                continue
            manifest[dest.relative_to(lib).as_posix()] = item["type"]

        print(f"  {item['id']:28s} -> {item['dest']}")

    meta = lib / "_metadata"
    meta.mkdir(parents=True, exist_ok=True)
    (meta / "index.yml").write_text(
        "# keys are RELATIVE to the library root; bodies are hardlinks\n"
        + "\n".join(f'? {k}\n: {v}' for k, v in sorted(manifest.items()))
        + "\n"
    )

    print(f"\nlinked {linked}, already-present {reused}, copied {len(copied)}, "
          f"in manifest {len(manifest)}")
    if copied:
        print("\nCOPIED, not linked (root-owned; protected_hardlinks forbids the link).\n"
              "These are the only items that can drift from their source:")
        for item_id, rel in copied:
            print(f"  {item_id}: {rel}")
    if failed:
        print(f"\n{len(failed)} FAILURE(S):")
        for path, why in failed:
            print(f"  {path}: {why}")
        return 1
    print("all placed files verified by inode identity against their sources")
    return 0


def phase_index(decl: dict, lib_path: Path) -> int:
    """Write a REAL metasmith index over the already-placed tree.

    The hand-rolled index.yml the stage phase writes is a placeholder for human
    eyes. metasmith's own format carries a schema, per-item parent lineage and a
    bundled copy of every type library, and only `Save()` produces it correctly
    -- so build it through the API rather than emitting YAML that merely looks
    right.
    """
    sys.path.insert(0, str(REPO / "src" / "metasmith" / "src"))
    from metasmith.models.libraries import DataInstanceLibrary  # noqa: E402

    plan, _ = build_plan(decl)
    lib = DataInstanceLibrary(location=lib_path)
    # A namespace binds to exactly one type library. These used to be read from
    # a repo-local data_types/, which meant staged items were typed against one
    # file while the transforms resolved against the library's own -- and since
    # metasmith matches endpoints by subset over a property set built from the
    # yaml, two same-named types with different text are different endpoints.
    # The failure is a wrong join, not an exception. Resolve through the library.
    sys.path.insert(0, str(REPO / "src"))
    from fabfos.library import resolve_library_root  # noqa: E402

    types_root = resolve_library_root() / "data_types"
    for ns in ("ecspr", "ref"):
        lib.AddTypeLibrary(types_root / f"{ns}.yml", namespace=ns)

    by_id = {}
    for item in plan:
        for src, rel in item["_members"]:
            dest = Path(item["dest"]) / rel if rel else Path(item["dest"])
            if not (lib_path / dest).exists():
                print(f"  MISSING ON DISK, skipped: {dest}")
                continue
            lib.AddItem(dest, item["type"])
            by_id.setdefault(item["id"], []).append(dest)

    # declared `inputs` become real parent lineage, so the library itself
    # records that uniref50.dmnd came from uniref50.fasta.gz -- provenance the
    # planner can traverse, not just prose in a record.
    for item in plan:
        for dep in item.get("inputs", []):
            if dep not in by_id:
                continue
            for child in by_id.get(item["id"], []):
                lib.AddParentsTo(child, [lib.Get(p) for p in by_id[dep]])

    lib.Save(update_types=True)
    print(f"\nindexed {len(lib.manifest)} items across {len(lib.types)} type namespaces")
    return 0


def phase_verify(lib_path: Path) -> int:
    """AC2: the library loads with integrity checking on."""
    sys.path.insert(0, str(REPO / "src" / "metasmith" / "src"))
    from metasmith.models.libraries import DataInstanceLibrary  # noqa: E402

    lib = DataInstanceLibrary.Load(lib_path, check_integrity=True)
    types = sorted({t for t in lib.manifest.values()})
    print(f"loaded {len(lib.manifest)} items, all present on disk")
    print(f"{len(types)} distinct types, all resolving:")
    for t in types:
        n = sum(1 for v in lib.manifest.values() if v == t)
        print(f"  {t:34s} {n:4d}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", action="store_true")
    ap.add_argument("--hash", action="store_true")
    ap.add_argument("--place", action="store_true")
    ap.add_argument("--index", action="store_true")
    ap.add_argument("--verify", action="store_true")
    ap.add_argument("--library", default=str(REPO / ".awm" / "data" / "ref"))
    args = ap.parse_args()

    if not (args.stage or args.hash or args.place or args.index or args.verify):
        ap.error("pick a phase: --stage, --hash, --place, --index or --verify")

    decl = load_declaration()
    lib = Path(args.library)

    if args.stage:
        print("=== stage ===")
        if (rc := phase_stage(decl, lib)):
            return rc
    if args.hash:
        print("=== hash ===")
        if (rc := phase_hash(decl, lib)):
            return rc
    if args.place:
        print("=== place ===")
        if (rc := phase_place(decl, lib)):
            return rc
    if args.index:
        print("=== index ===")
        if (rc := phase_index(decl, lib)):
            return rc
    if args.verify:
        print("=== verify ===")
        if (rc := phase_verify(lib)):
            return rc
    return 0


if __name__ == "__main__":
    sys.exit(main())

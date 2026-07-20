#!/usr/bin/env python3
"""Resolve every `*.oci` container reference to an immutable registry digest.

The library pins images by tag (`docker://staphb/bbtools:39.49`). Tags are
mutable: the same reference pulls different bits six months from now, on
another machine. This walks the library's container resources, asks the
registry what the tag currently points at, and writes one provenance record
per image under `provenance/containers/`.

An image that cannot be resolved is not an error -- it is the finding. An
untagged reference, or a `:latest` tag that exists only in this machine's
local docker daemon, cannot work on a fresh machine at all. Those get a
record with `resolved: false` and a `blocker` explaining why.

Usage:
    python resolve_container_digests.py <library_root> <provenance_dir>
"""
import json
import subprocess
import sys
from pathlib import Path

TIMEOUT = 120


def parse_ref(raw: str) -> str:
    """`docker://quay.io/x/y:1.0` -> `quay.io/x/y:1.0`."""
    return raw.strip().removeprefix("docker://").strip()


def resolve(ref: str) -> tuple[str | None, str | None]:
    """Return (digest, error). Digest is the manifest descriptor digest."""
    try:
        proc = subprocess.run(
            ["docker", "manifest", "inspect", "--verbose", ref],
            capture_output=True, text=True, timeout=TIMEOUT,
        )
    except subprocess.TimeoutExpired:
        return None, f"registry query timed out after {TIMEOUT}s"
    if proc.returncode != 0:
        return None, proc.stderr.strip().splitlines()[-1] if proc.stderr.strip() else "unknown error"
    try:
        payload = json.loads(proc.stdout)
    except json.JSONDecodeError as e:
        return None, f"could not parse manifest: {e}"
    # a multi-arch tag yields a list, one entry per platform; they share a
    # manifest-list digest, so any entry's Descriptor.digest is the pin.
    entry = payload[0] if isinstance(payload, list) else payload
    digest = entry.get("Descriptor", {}).get("digest")
    if not digest:
        return None, "manifest carried no Descriptor.digest"
    return digest, None


def classify(ref: str) -> str | None:
    """Pre-flight blockers that make a registry query pointless."""
    name = ref.split("/")[-1]
    if ":" not in name:
        return (
            "no tag -- the reference floats on the registry's implicit :latest "
            "and pins nothing"
        )
    if name.endswith(":latest") and "/" not in ref:
        return (
            "bare `:latest` on no registry -- this image exists only in the "
            "local docker daemon and is unobtainable on a fresh machine"
        )
    return None


def main() -> int:
    lib_root = Path(sys.argv[1])
    out_dir = Path(sys.argv[2])
    out_dir.mkdir(parents=True, exist_ok=True)

    oci_files = sorted((lib_root / "resources" / "containers").glob("*.oci"))
    if not oci_files:
        print(f"no *.oci under {lib_root}/resources/containers", file=sys.stderr)
        return 1

    resolved, blocked, curated = 0, [], []
    for path in oci_files:
        name = path.stem
        raw = path.read_text()
        ref = parse_ref(raw)

        # This writer emits a MINIMAL record and hardcodes `origin:
        # public-registry` / `built_by_us: false`. For an image we build
        # ourselves that would silently drop the hand-authored provenance --
        # the recipe path, the build assertions, how the tag was resolved --
        # which is the whole reason those records exist. Refuse to touch them;
        # a human edits those by hand.
        existing = out_dir / f"{name}.yml"
        if existing.exists() and "built_by_us: true" in existing.read_text():
            curated.append(name)
            print(f"{name:28s} SKIPPED (hand-authored, built_by_us: true)", flush=True)
            continue

        blocker = classify(ref)
        digest, error = (None, blocker) if blocker else resolve(ref)
        if digest:
            resolved += 1
            status = digest[:19] + "..."
        else:
            blocked.append((name, error))
            status = f"BLOCKED: {error}"
        print(f"{name:28s} {status}", flush=True)

        lines = [
            f"id: containers.{name}",
            f"name: {name}",
            "kind: container",
            f"resource: resources/containers/{name}.oci",
            f"reference: {ref}",
            f"resolved: {'true' if digest else 'false'}",
        ]
        if digest:
            lines += [
                f"digest: {digest}",
                f"pinned_reference: {ref.split(':')[0] if ref.count(':') else ref}@{digest}",
            ]
        else:
            lines += [f"blocker: >-\n  {error}"]
        lines += [
            "origin: public-registry",
            "built_by_us: false",
            "",
        ]
        (out_dir / f"{name}.yml").write_text("\n".join(lines))

    print(f"\nresolved {resolved}/{len(oci_files) - len(curated)} (of the records this script owns)")
    for name, why in blocked:
        print(f"  BLOCKED {name}: {why}")
    if curated:
        print(f"  left alone, hand-authored: {', '.join(curated)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

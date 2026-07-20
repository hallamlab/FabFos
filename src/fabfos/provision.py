"""Provision the per-tool conda/mamba environments a MAMBA run needs.

Under ``Runtime.MAMBA`` each tool is launched as ``mamba run -n <env> <cmd>``
(see metasmith ``env/environment.py``). Unlike the container runtimes — where
``docker run`` / the apptainer pull step materializes the image on demand — a
mamba env does **not** come into being on first use: it must already exist.

The unified ``<tool>.env.yml`` resources in the library carry a full conda
spec (``conda: {name, channels, dependencies}``) alongside the container
pointer. This module reads those specs, groups them by env name (several tools
share one env, e.g. the core bio tools all live in ``fabfos-bio``), and ensures
each named env exists with the union of the declared dependencies — so a
freshly-installed FabFos can stand up its own tool environments with no
containers involved.

Idempotent by design: an env that already exists is left alone unless
``force`` is set, so the common case (envs already provisioned) is instant and
a run is never blocked waiting on a solver.
"""
from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path

import yaml


@dataclass
class ToolEnvSpec:
    name: str
    channels: list[str] = field(default_factory=list)
    dependencies: list[str] = field(default_factory=list)
    sources: list[str] = field(default_factory=list)  # env.yml files that fed this env


def _parse_env_yml(path: Path) -> tuple[str, list[str], list[str]] | None:
    """Return (env_name, channels, dependencies) for a ``*.env.yml`` resource.

    Accepts both the full-block form (``conda: {name, channels, dependencies}``)
    and the shorthand where ``conda`` is a bare string naming the env (in which
    case there are no deps to install — the env is expected to pre-exist or be
    provisioned elsewhere). Returns ``None`` for a legacy bare-URI ``.oci`` or a
    doc without a conda section.
    """
    try:
        doc = yaml.safe_load(path.read_text())
    except yaml.YAMLError:
        return None
    if not isinstance(doc, dict):
        return None
    conda = doc.get("conda")
    if isinstance(conda, str) and conda.strip():
        return conda.strip(), [], []
    if isinstance(conda, dict) and conda.get("name"):
        name = str(conda["name"]).strip()
        channels = [str(c) for c in (conda.get("channels") or [])]
        deps = [str(d) for d in (conda.get("dependencies") or [])]
        return name, channels, deps
    return None


def collect_tool_environments(lib_root: Path) -> dict[str, ToolEnvSpec]:
    """Scan the library's container resources for ``*.env.yml`` and group the
    conda specs by env name, unioning channels + dependencies."""
    envs: dict[str, ToolEnvSpec] = {}
    containers = lib_root / "resources" / "containers"
    if not containers.is_dir():
        return envs
    for p in sorted(containers.glob("*.env.yml")):
        parsed = _parse_env_yml(p)
        if parsed is None:
            continue
        name, channels, deps = parsed
        spec = envs.setdefault(name, ToolEnvSpec(name=name))
        for c in channels:
            if c not in spec.channels:
                spec.channels.append(c)
        for d in deps:
            if d not in spec.dependencies:
                spec.dependencies.append(d)
        spec.sources.append(p.name)
    return envs


def _existing_envs() -> set[str]:
    mamba = shutil.which("mamba") or shutil.which("conda")
    if not mamba:
        return set()
    try:
        out = subprocess.run([mamba, "env", "list"], capture_output=True, text=True, check=True).stdout
    except subprocess.CalledProcessError:
        return set()
    names: set[str] = set()
    for line in out.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        names.add(line.split()[0])
    return names


@dataclass
class ProvisionReport:
    created: list[str] = field(default_factory=list)
    skipped: list[str] = field(default_factory=list)   # already existed
    failed: list[tuple[str, str]] = field(default_factory=list)


def provision_tool_environments(
    lib_root: Path,
    *,
    force: bool = False,
    dry_run: bool = False,
    log=print,
) -> ProvisionReport:
    """Ensure every conda env named by the library's ``*.env.yml`` resources
    exists. Creates missing envs from the union of their declared deps; leaves
    existing envs untouched unless ``force``. Never raises for a single env's
    failure — it is collected in the report so one bad solve does not sink the
    whole run.
    """
    report = ProvisionReport()
    envs = collect_tool_environments(lib_root)
    if not envs:
        log("provision: no *.env.yml tool environments found; nothing to do")
        return report

    mamba = shutil.which("mamba") or shutil.which("conda")
    if not mamba and not dry_run:
        raise RuntimeError("provision: neither `mamba` nor `conda` found on PATH")

    existing = _existing_envs()
    for name, spec in envs.items():
        if name in existing and not force:
            log(f"provision: env [{name}] already present (from {spec.sources}) — skipping")
            report.skipped.append(name)
            continue
        if not spec.dependencies:
            log(f"provision: env [{name}] names no dependencies (shorthand pointer); "
                f"cannot create it from spec — ensure it exists out of band")
            report.failed.append((name, "no dependencies declared in env.yml specs"))
            continue
        cmd = [mamba, "create", "-y", "-n", name]
        for c in spec.channels:
            cmd += ["-c", c]
        cmd += list(spec.dependencies)
        log(f"provision: creating env [{name}] with {spec.dependencies} "
            f"(channels={spec.channels}) [from {spec.sources}]")
        if dry_run:
            log("  (dry-run) " + " ".join(cmd))
            report.created.append(name)
            continue
        try:
            subprocess.run(cmd, check=True)
            report.created.append(name)
        except subprocess.CalledProcessError as e:
            log(f"provision: FAILED to create env [{name}]: {e}")
            report.failed.append((name, str(e)))
    return report

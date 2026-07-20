import os, sys
import stat
from pathlib import Path
import yaml

HERE = Path(os.path.realpath(__file__)).parent
sys.path = list({str(HERE.joinpath("../").absolute())} | set(sys.path))

# constants from the package
from setup import USER, NAME, __version__ as VERSION, ENTRY_POINTS, SHORT_SUMMARY  # type: ignore

# ------------------------------------------------------------------
# dependencies (from envs/base.yml). conda recipes can't carry pip deps.
with open(HERE.joinpath("../envs/base.yml")) as y:
    raw_deps = yaml.safe_load(y)


def _parse_deps(level: list, compiled: str, depth: int):
    tabs = "  " * depth
    for item in level:
        if not isinstance(item, str) or item in {"pip"}:
            continue
        compiled += f"{tabs}- {item}\n"
    return compiled[:-1]


reqs = _parse_deps(raw_deps["dependencies"], "", 2)
python_dep = [d for d in raw_deps["dependencies"] if isinstance(d, str) and d.startswith("python=")]
python_ver = _parse_deps(python_dep or ["python=3.12"], "", 2)

# ------------------------------------------------------------------
# entry points
entry_points = "".join(f"{'  '*2}- {e}\n" for e in ENTRY_POINTS)[:-1]

# ------------------------------------------------------------------
# render recipe
with open(HERE.joinpath("meta_template.yaml")) as f:
    template = f.read()
for k, v in {
    "USER": USER,
    "NAME": NAME,
    "SHORT_SUMMARY": SHORT_SUMMARY,
    "VERSION": VERSION,
    "ENTRY": entry_points,
    "REQUIREMENTS": reqs,
    "PYTHON": python_ver,
}.items():
    template = template.replace(f"<{k}>", v)
with open(HERE.joinpath("meta.yaml"), "w") as f:
    f.write(template)

build_file = HERE.joinpath("call_build.sh")
channels = " ".join(f"-c {ch}" for ch in raw_deps["channels"])
with open(build_file, "w") as f:
    f.write(
        'HERE=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )\n'
        f"conda mambabuild {channels} --output-folder $HERE/../conda_build $HERE/\n"
    )
os.chmod(build_file, os.stat(build_file).st_mode | stat.S_IEXEC)
print(f"wrote {HERE/'meta.yaml'} and {build_file}")

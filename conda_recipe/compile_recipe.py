import os, sys
import stat
from pathlib import Path
import yaml

HERE = Path(os.path.realpath(__file__)).parent
sys.path = list(set([
    str(HERE.joinpath("../").absolute())
]+sys.path))

# import constants from setup.py
from setup import NAME, VERSION, ENTRY_POINTS


# ======================================================
# parse dependencies
with open(HERE.joinpath(f"../envs/base.yml")) as y:
    raw_deps = yaml.safe_load(y)
def _parse_deps(level: list, compiled: str, depth: int):
    tabs_space = "  "*depth
    for item in level:
        if isinstance(item, str):
            compiled += f"{tabs_space}- {item}\n"
        else:
            k, v = list(item.items())[0]
            compiled += f"{tabs_space}- {k}:\n"
            compiled = _parse_deps(v, compiled, depth+1)
    return compiled
reqs = _parse_deps(raw_deps["dependencies"], "", 2)[:-1] # remove trailing \n


# ======================================================
# entry points

entry_points = ""
for e in ENTRY_POINTS:
    tabs_space = "  "*2
    entry_points += f"{tabs_space}- {e}\n"
entry_points = entry_points[:-1] # remove trailing \n


# ======================================================
# path to tar archive of source code

dist_path = Path(os.path.abspath(HERE.joinpath("../dist")))
tar_path = [dist_path.joinpath(f) for f in os.listdir(dist_path) if VERSION in f and ".tar.gz" in f][0]


# ======================================================
# generate recipe files

with open(HERE.joinpath("meta_template.yaml")) as f:
    template = "".join(f.readlines())
meta_values = {
    "NAME": NAME,
    "VERSION": VERSION,
    "ENTRY": entry_points,
    "REQUIREMENTS": reqs,
    "TAR": f"file://{tar_path}"
}
for k, v in meta_values.items():
    template = template.replace(f"<{k}>", v)
with open(HERE.joinpath("meta.yaml"), "w") as f:
    f.write(template)

build_file = HERE.joinpath("call_build.sh")
with open(build_file, "w") as f:
    channels = " ".join(f"-c {ch}" for ch in raw_deps["channels"])
    _here = 'HERE=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )'
    f.write(f"""\
        {_here}
        conda mambabuild {channels} --output-folder $HERE/../conda_build $HERE/
    """.replace("    ", ""))
st = os.stat(build_file)
os.chmod(build_file, st.st_mode | stat.S_IEXEC)

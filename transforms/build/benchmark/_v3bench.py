"""Shared paths and helpers for the benchmark v3 answer-key build.

Everything resolves through `canon`, which resolves through the data library.
No absolute path appears here -- that is the whole point of the migration, and
v1's own `_bench.py` is the counter-example: it pins

    BENCH = DATA / "ecspr/benchmark/v1"

so the scorer and every builder were welded to v1 and could not be pointed at
another version without editing the module.

The v3 tree is SELF-CONTAINED: X, the contract shape, the ground truth, the
baseline and v1's provenance all live inside it, so nothing here reaches
sideways into a sibling version.
"""
from __future__ import annotations

import hashlib
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO / "src"))

from fabfos import canon  # noqa: E402

# ---- the v3 tree, through the library ----
ROOT = Path(canon.BENCH_V3_ROOT)
X_DIR = Path(canon.BENCH_V3_X)
OBS_DIR = Path(canon.BENCH_V3_OBSERVATIONS)
GT_DIR = Path(canon.BENCH_V3_GROUND_TRUTH)
DEC_DIR = Path(canon.BENCH_V3_DECISIONS)
CONTRACT_DIR = Path(canon.BENCH_V3_CONTRACT)

# The answer key is the thing being BUILT, so it is addressed as a destination
# rather than resolved through canon -- canon.BENCH_V3_Y deliberately raises
# until the key exists, which is correct for readers and useless for a writer.
#
# It is built INTO THE REPO, not into the library, for the same reason
# target_resolution.tsv is: the library is built from declarations by the
# staging script, and a file written into it directly is one that no record
# describes and no hash covers. Placement then hardlinks this into
# `validation/benchmark/v3/Y`, where canon.BENCH_V3_Y finds it.
BUILD_OUT = Path(__file__).resolve().parent / "v3_build"
Y_DIR = BUILD_OUT / "Y"

ELEMENTS = ["C", "N", "S", "P"]
FACETS = ["netA_iML1515", "netA_iECDH10B", "netB_iML1515", "netB_iECDH10B"]

# Currency metabolites, as v1's phenotype builder defined them. A reaction's
# "target" is its dedicated product, so the ubiquitous cofactors have to come
# out or every condition would target ATP. Acetyl-CoA is deliberately KEPT:
# it is a real biosynthetic endpoint, not bulk currency.
CURRENCY_NAMES = {
    "h+", "h2o", "water", "co2", "carbon dioxide", "atp", "adp", "amp",
    "phosphate", "diphosphate", "triphosphate", "nad+", "nadh", "nadp+",
    "nadph", "coa", "coenzyme a", "o2", "oxygen", "nh4+", "ammonium",
    "so4(2-)", "sulfate", "h", "pi", "ppi", "pppi", "fad", "fadh2",
}


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 22), b""):
            h.update(chunk)
    return h.hexdigest()


def write_tsv(df, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False)
    print(f"  wrote {path.relative_to(ROOT) if ROOT in path.parents else path}"
          f"  {len(df)} rows")


# ---------------------------------------------------------------------------
# name normalisation
# ---------------------------------------------------------------------------

_GREEK = {
    "alpha": "a", "beta": "b", "gamma": "g", "delta": "d",
    "α": "a", "β": "b", "γ": "g", "δ": "d",
}
_STRIP_PUNCT = re.compile(r"[\s\-_,'’\"()\[\]{}]+")


def norm_name(s: str) -> str:
    """Fold a chemical name to a comparison key.

    Deliberately AGGRESSIVE about punctuation and whitespace and deliberately
    CONSERVATIVE about stereochemistry: `L-` and `D-` are kept, because
    L-serine and D-serine are different metabolites and folding them together
    would silently retarget a condition. The greek-letter fold is safe because
    the alternatives are spelling variants of one name, not distinct species.
    """
    s = s.strip().lower()
    for k, v in _GREEK.items():
        s = s.replace(k, v)
    s = _STRIP_PUNCT.sub("", s)
    return s


# Cells that name more than one metabolite. The normalized LOF column collapses
# these to a single label, which throws away information the sparse
# expectations format can express perfectly well as several target rows.
# Every separator REQUIRES surrounding whitespace except `/` and `;`, and that
# is not fussiness -- chemical names contain the bare characters. A plain `,`
# would split `2,3-dihydroxybenzoate` into two non-metabolites, and a plain `+`
# would split `NAD+` and `h+`. Both would resolve to nothing and be reported as
# coverage misses, so the damage would look like missing data rather than a
# parser bug.
_SPLIT = re.compile(r"\s+\+\s+|\s*/\s*|\s*;\s*|\s+and\s+|\s+or\s+|,\s+")


def split_components(raw: str) -> list[str]:
    """`L-isoleucine + L-valine` -> two targets. `FMN/FAD` -> two targets.

    A compound cell is a condition with several legitimate targets, not an
    unparseable one. Splitting keeps them all; the normalized column keeps one.
    """
    parts = [p.strip() for p in _SPLIT.split(raw) if p.strip()]
    # a lone "none"/"n/a" style cell is not a target at all
    return [p for p in parts if norm_name(p) not in {"none", "na", "nan", ""}]

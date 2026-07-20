#!/bin/bash
# FabFos dev/build automation.
#
# FabFos is a thin metasmith front end: the planner/executor comes from the
# `metasmith` conda package, and the fosmid pipeline definition (transforms,
# data types, container/conda env resources) is the metasmith library bundled
# into the wheel as `fabfos/_library`.
set -e
# this script sits at the repo root; the package is in ./src/fabfos (src-layout)
REPO=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
HERE="$REPO"
LIB_SRC="$REPO/src/metasmith_libraries"
LIB_DST="$REPO/src/fabfos/_library"

case $1 in
    --ibase) # create the dev conda env
        mamba env create --no-default-packages -f "$HERE/envs/base.yml"
    ;;
    -b|--bundle-library) # copy the metasmith library into the package for shipping
        echo "bundling metasmith library: $LIB_SRC -> $LIB_DST"
        rm -rf "$LIB_DST"
        mkdir -p "$LIB_DST"
        # only the pieces the planner loads at runtime
        for sub in data_types resources transforms; do
            cp -r "$LIB_SRC/$sub" "$LIB_DST/$sub"
        done
    ;;
    -bp|--build-pip) # build the wheel/sdist (bundle first)
        "$HERE/dev.sh" --bundle-library
        cd "$HERE"
        rm -rf build dist *.egg-info
        python -m build
    ;;
    -bc|--build-conda) # compile + build the conda package
        "$HERE/dev.sh" --bundle-library
        python "$HERE/conda_recipe/compile_recipe.py"
        "$HERE/conda_recipe/call_build.sh"
    ;;
    -r|--run) # run the CLI from source (dev): ./dev.sh -r --plan-only ...
        shift
        PYTHONPATH="$REPO/src:$REPO/src/metasmith/src" \
        FABFOS_LIBRARY="${FABFOS_LIBRARY:-$LIB_SRC}" \
        python -m fabfos "$@"
    ;;
    *)
        echo "usage: dev.sh [--ibase|-b|-bp|-bc|-r ...]"
    ;;
esac

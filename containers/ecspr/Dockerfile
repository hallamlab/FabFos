ARG CONDA_ENV=for_container

# ECSPr numeric core: SMW/Woodbury effective-conductance solve, per-gene ablation,
# and the full-mixture-SF significance scorer. Consumed by the fabfos metasmith
# library (transforms/ecspr/*) as `containers::ecspr.oci`.
#
# Why this image exists rather than reusing the pre-existing one: that image is
# missing pins this lane needs, and its torch is built for a CUDA the target
# cluster's V100s do not take. See load/env.yml for the numpy note -- in
# particular, the "numpy 1.26 cannot unpickle the base graphs" rationale that gets
# repeated for this image is false, and is not why we are here.

# Base is PINNED by digest-bearing tag rather than :latest. The sibling recipes in
# this repo track :latest and have silently rotted against it -- mamba 2.x dropped
# `--no-default-packages`, so they no longer build. An image whose base moves under
# it is not a reproducible image.
FROM condaforge/miniforge3:25.3.0-3 AS build-env
ARG CONDA_ENV

# avoid permission issues when manipulating added files under /opt
RUN old_umask=`umask` \
    && umask 0000 \
    && umask $old_umask

ADD ./load/env.yml /opt/env.yml
RUN --mount=type=cache,target=/opt/conda/pkgs \
    mamba env create -n ${CONDA_ENV} -f /opt/env.yml
ENV PATH=/opt/conda/envs/${CONDA_ENV}/bin:$PATH

# torch, pinned to the cu128 wheel index. NON-NEGOTIABLE: the default CUDA build
# breaks the V100 driver on the target cluster. The solve and the ablation are
# CPU-intractable on the dense maps, so this is not a nice-to-have fallback --
# without a working GPU build those two steps do not run at all.
RUN --mount=type=cache,target=/root/.cache/pip \
    pip install --index-url https://download.pytorch.org/whl/cu128 torch

# Fail the BUILD, not a queued GPU job six hours later, if the pins regressed.
RUN python -c "\
import numpy, scipy, sklearn, pandas, networkx, pyarrow, torch; \
assert int(numpy.__version__.split('.')[0]) >= 2, f'numpy {numpy.__version__} < 2'; \
assert torch.version.cuda == '12.8', f'torch CUDA {torch.version.cuda} != 12.8'; \
print('ecspr image OK: numpy', numpy.__version__, '| torch', torch.__version__, \
      '| cuda', torch.version.cuda, '| nx', networkx.__version__)"

# Singularity uses tini, but raises warnings because the -s flag is not used.
# This used to be `ADD ./load/tini /tini` against a vendored binary that is
# gitignored and absent, so the recipe could not build from a fresh clone at
# all -- which made the `recipe:` field in provenance/containers/ecspr.yml a
# false claim. Fetch a pinned release and verify it, so the recipe is
# self-contained and the bytes are the ones recorded here.
ADD https://github.com/krallin/tini/releases/download/v0.19.0/tini-static-amd64 /tini
RUN echo "c5b0666b4cb676901f90dfcb37106783c5fe2077b04590973b885950611b30ee  /tini" | sha256sum -c - \
    && chmod +x /tini
ENTRYPOINT ["/tini", "-s", "-g", "--"]

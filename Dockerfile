ARG CONDA_ENV=fabfos

FROM condaforge/mambaforge as build-env
# https://mamba.readthedocs.io/en/latest/user_guide/mamba.html

# scope var from global
ARG CONDA_ENV
COPY ./envs/base.yml /opt/env.yml
# Create conda environment:
# name must match what is in conda.yml
RUN mamba env create --no-default-packages -n $CONDA_ENV -f /opt/env.yml\
    && mamba clean -afy

# use a smaller runtime image
# jammy is ver. 22.04 LTS
# https://wiki.ubuntu.com/Releases
FROM ubuntu:jammy
# scope var from global
ARG CONDA_ENV
COPY --from=build-env /opt/conda/envs/${CONDA_ENV} /opt/conda/envs/${CONDA_ENV}
ENV PATH /opt/conda/envs/${CONDA_ENV}/bin:$PATH

COPY ./src/fabfos /app/fabfos
RUN echo "python -m fabfos \$@" >/usr/bin/fabfos && chmod +x /usr/bin/fabfos
ENV PYTHONPATH /app:$PYTHONPATH

# Singularity uses tini, but raises warnings
# we set it up here correctly for singularity
ENV TINI_VERSION v0.19.0
ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /tini
RUN chmod +x /tini
    
# singularity doesn't use the -s flag, and that causes warnings
ENTRYPOINT ["/tini", "-s", "--"]

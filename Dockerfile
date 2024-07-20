ARG CONDA_ENV=fabfos

FROM condaforge/miniforge3
# https://mamba.readthedocs.io/en/latest/user_guide/mamba.html

# scope var from global
ARG CONDA_ENV
COPY ./envs/base.yml /opt/env.yml
RUN --mount=type=cache,target=/opt/conda/pkgs \
    mamba env create -n ${CONDA_ENV} --no-default-packages -f /opt/env.yml

# # use a smaller runtime image
# # jammy is ver. 22.04 LTS
# # https://wiki.ubuntu.com/Releases
# FROM ubuntu:jammy
# # scope var from global
# ARG CONDA_ENV
# COPY --from=build-env /opt/conda/envs/${CONDA_ENV} /opt/conda/envs/${CONDA_ENV}
ENV PATH /opt/conda/envs/${CONDA_ENV}/bin:$PATH

ADD ./src/fabfos /app/fabfos
RUN echo "python -m fabfos \$@" >/usr/bin/fabfos && chmod +x /usr/bin/fabfos
ENV PYTHONPATH /app:$PYTHONPATH

# Singularity uses tini, but raises warnings
# we set it up here correctly for singularity
ADD ./lib/tini /tini
RUN chmod +x /tini
    
# singularity doesn't use the -s flag, and that causes warnings.
# -g kills process group on ctrl+C
ENTRYPOINT ["/tini", "-s", "-g", "--"]

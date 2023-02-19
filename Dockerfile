FROM condaforge/mambaforge as build-env
# https://mamba.readthedocs.io/en/latest/user_guide/mamba.html

# Create conda environment:
# name must match what is in conda.yml
ENV CONDA_ENV "fabfos"
COPY ./envs/conda.yml /opt/
RUN mamba env create --no-default-packages -f /opt/conda.yml

# Singularity uses tini, but raises warnings
# we set it up here correctly for singularity
ENV TINI_VERSION v0.19.0
ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /tini
RUN chmod +x /tini

# use a smaller runtime image
FROM ubuntu:22.10
ENV CONDA_ENV "fabfos"
COPY --from=build-env /opt/conda/envs/${CONDA_ENV} /opt/conda/envs/${CONDA_ENV}
COPY --from=build-env /tini /tini
ENV PATH /opt/conda/envs/${CONDA_ENV}/bin:$PATH

COPY ./src /app
ENV PATH /app:$PATH

# singularity doesn't use the -s flag, and that causes warnings
ENTRYPOINT ["/tini", "-s", "--"]

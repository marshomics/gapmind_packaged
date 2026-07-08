# GapMind container: pinned conda env + patched PaperBLAST, ready for the
# prebuilt analysis path. Build with containers/build.sh (produces a Docker image
# and, for HPC, a Singularity/Apptainer .sif).
FROM condaforge/miniforge3:24.3.0-0

LABEL org.opencontainers.image.title="gapmind" \
      org.opencontainers.image.description="GapMind (amino-acid + carbon) analysis, prebuilt-DB path" \
      org.opencontainers.image.source="https://github.com/morgannprice/PaperBLAST"

# Non-conda helpers needed at build/run time.
RUN apt-get update && apt-get install -y --no-install-recommends \
        git ca-certificates procps xz-utils \
    && rm -rf /var/lib/apt/lists/*

# Pinned conda environment (installed into base; tools land on /opt/conda/bin).
COPY environment.yml /tmp/environment.yml
RUN conda env update -n base -f /tmp/environment.yml && conda clean -afy
ENV PATH=/opt/conda/bin:$PATH

# PaperBLAST/GapMind code, pinned. Override at build time:
#   docker build --build-arg PAPERBLAST_REF=<sha> ...
ARG PAPERBLAST_REF=master
COPY setup_paperblast.sh /usr/local/bin/setup_paperblast.sh
RUN bash /usr/local/bin/setup_paperblast.sh /opt/PaperBLAST "${PAPERBLAST_REF}"

# The workflow's bin/ scripts read $GAPMIND_DIR; put the code's bin/ on PATH.
ENV GAPMIND_DIR=/opt/PaperBLAST
ENV PATH=/opt/PaperBLAST/bin:$PATH
# Deterministic locale/threads for reproducible runs.
ENV LC_ALL=C.UTF-8 LANG=C.UTF-8

WORKDIR /work
CMD ["bash"]

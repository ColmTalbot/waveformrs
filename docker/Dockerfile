# syntax=docker/dockerfile:1
#
# Linux environment for testing and benchmarking waveformrs against the
# conda-forge build of lalsimulation: a pinned stable Rust toolchain and a
# micromamba environment with python-lalsimulation 6.2.0 (docker/environment.yml).
# The build context is docker/; the repository is not copied in, docker/run.sh
# mounts it at /work.
#
# Multi-platform (linux/arm64 and linux/amd64): build one platform at a time,
#     docker/run.sh --platform arm64 build
# See docs/linux.md.

ARG RUST_VERSION=1.90.0
ARG MICROMAMBA_VERSION=1.5.10

FROM mambaorg/micromamba:${MICROMAMBA_VERSION}-bookworm-slim AS micromamba

FROM rust:${RUST_VERSION}-slim-bookworm

# curl for the conda downloads; gcc/binutils (objdump) come with the image
RUN apt-get update \
    && apt-get install -y --no-install-recommends curl \
    && rm -rf /var/lib/apt/lists/*
RUN rustup component add clippy rustfmt

COPY --from=micromamba /bin/micromamba /usr/local/bin/micromamba
ENV MAMBA_ROOT_PREFIX=/opt/conda
COPY environment.yml conda-install.sh /tmp/conda/
RUN --mount=type=cache,target=/var/cache/conda,sharing=locked \
    CONDA_DOWNLOAD_CACHE=/var/cache/conda sh /tmp/conda/conda-install.sh /tmp/conda/environment.yml \
    && micromamba clean -a -y \
    && rm -rf /tmp/conda

# The conda environment goes last on PATH so that its linker (ld_impl) does
# not shadow the system one that rustc's linker driver uses.
ENV PATH=/usr/local/cargo/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/opt/conda/envs/lal/bin \
    OMP_NUM_THREADS=1 \
    OPENBLAS_NUM_THREADS=1 \
    CARGO_TARGET_DIR=/target \
    LAL_DATA_PATH=/work/lal_data \
    PYTHONDONTWRITEBYTECODE=1
WORKDIR /work

#!/bin/sh
# Create the conda environment "lal" from an environment file without
# hammering the conda-forge CDN: solve first, then download the packages one
# at a time (at most one request every ${CONDA_DOWNLOAD_INTERVAL:-1} s), then
# install from the package cache.
#
#     sh conda-install.sh environment.yml
#
# Downloads are kept in ${CONDA_DOWNLOAD_CACHE} if set (a BuildKit cache
# mount in docker/Dockerfile), so a failed or repeated build does not fetch
# them again.
set -eu

env_file=$1
interval=${CONDA_DOWNLOAD_INTERVAL:-1}
pkgs=${MAMBA_ROOT_PREFIX:?}/pkgs
downloads=${CONDA_DOWNLOAD_CACHE:-$pkgs}
# a single package cache (the offline install insists that every cache exists)
export CONDA_PKGS_DIRS="$pkgs"
plan=$(mktemp)
mkdir -p "$pkgs" "$downloads"

# The solve downloads the repodata of conda-forge (two requests).
micromamba create -y -n lal -f "$env_file" --dry-run --json >"$plan"

grep -o '"url": *"[^"]*"' "$plan" | sed 's/.*"\(http[^"]*\)"/\1/' | sort -u |
    while read -r url; do
        name=${url##*/}
        if [ ! -s "$downloads/$name" ]; then
            echo "fetching $name"
            curl -fsSL --retry 3 --retry-delay 15 --retry-all-errors -o "$downloads/$name.part" "$url"
            mv "$downloads/$name.part" "$downloads/$name"
            sleep "$interval"
        fi
        if [ "$downloads" != "$pkgs" ]; then
            cp "$downloads/$name" "$pkgs/$name"
        fi
    done

# Everything is cached now: no further downloads.
micromamba create -y -n lal -f "$env_file" --offline
micromamba env export -n lal --explicit >"$MAMBA_ROOT_PREFIX/lal-explicit.txt"
rm -f "$plan"

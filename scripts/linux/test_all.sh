#!/usr/bin/env bash
# The test matrix of docs/linux.md, run from the host: cargo test --release
# against the committed (osx-arm64) reference data and against reference data
# from each container's own lalsimulation, on linux/arm64 and linux/amd64.
#
#     scripts/linux/test_all.sh [arm64] [amd64]
set -uo pipefail
run=$(dirname "${BASH_SOURCE[0]}")/../../docker/run.sh
[[ $# -eq 0 ]] && set -- arm64 amd64
status=0
for platform in "$@"; do
    echo "== linux/$platform, committed osx-arm64 reference data"
    "$run" --platform "$platform" test || status=1
    echo "== linux/$platform, reference data from the container's lalsimulation"
    "$run" --platform "$platform" reference || status=1
done
exit $status

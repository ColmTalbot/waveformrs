#!/usr/bin/env bash
# The test matrix of docs/linux.md, run from the host. On each platform:
# cargo test --release against the committed (osx-arm64) reference data,
# against reference data from the container's own lalsimulation, and
# without default features (no `fft`); then the build and lint checks.
#
#     scripts/linux/test_all.sh [arm64] [amd64]
set -uo pipefail
run=$(dirname "${BASH_SOURCE[0]}")/../../docker/run.sh
[[ $# -eq 0 ]] && set -- arm64 amd64
status=0
for platform in "$@"; do
    echo "== linux/$platform, committed osx-arm64 reference data"
    "$run" --platform "$platform" test --no-fail-fast -- --nocapture || status=1
    echo "== linux/$platform, reference data from the container's lalsimulation"
    "$run" --platform "$platform" reference --no-fail-fast -- --nocapture || status=1
    echo "== linux/$platform, --no-default-features"
    # (not the doctests: those of src/conversions need the fft feature)
    "$run" --platform "$platform" test --no-default-features --lib --tests --no-fail-fast || status=1
    echo "== linux/$platform, build warnings, clippy and rustfmt"
    "$run" --platform "$platform" check || status=1
done
echo "status $status"
exit $status

#!/usr/bin/env bash
# Generate the reference data of tests/data with the lalsimulation of this
# machine (into /target/reference/<platform> by default, not tests/data) and
# run the test suite against them:
#
#     scripts/linux/reference.sh [--generate-only] [cargo test arguments]
#
# The tests read the data from WAVEFORMRS_REFERENCE_DATA and then apply the
# tolerances of the reference platform (tests/common/mod.rs).
# scripts/linux/compare_reference.py compares two such directories (e.g.
# this platform's data with the committed osx-arm64 data).
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.."

platform=${WAVEFORMRS_PLATFORM:-linux-$(uname -m)}
out=${WAVEFORMRS_REFERENCE_DATA:-/target/reference/$platform}
generate_only=
if [[ ${1:-} == --generate-only ]]; then
    generate_only=1
    shift
fi
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
if [[ ! -f $out/spintaylor.json ]]; then
    mkdir -p "$out"
    suites=()
    if [[ ! -f lal_data/NRSur7dq4_v1.0.h5 ]]; then
        # the NRSur7dq4 suites need the data file (their tests skip without it)
        suites=(--suites $(python -c 'import sys; sys.path.insert(0, "scripts")
import generate_reference_data as g
print(" ".join(s for s in g.SUITES if not s.startswith("nrsur7dq4")))'))
    fi
    LAL_DATA_PATH=lal_data python scripts/generate_reference_data.py --outdir "$out" \
        ${suites[@]+"${suites[@]}"}
    python scripts/generate_spintaylor_reference.py --outdir "$out"
fi
[[ -n $generate_only ]] && exit 0
WAVEFORMRS_REFERENCE_DATA=$out cargo test --release "$@"

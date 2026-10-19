#!/usr/bin/env bash
# Re-run every measurement used by the technical report, one after another
# (run on an otherwise idle machine), then rebuild the report.
#
#     bash report/measure_all.sh
set -euo pipefail
cd "$(dirname "$0")/.."
export OMP_NUM_THREADS=1
export LAL_DATA_PATH="${LAL_DATA_PATH:-$PWD/lal_data}"
export NRSUR7DQ4_DATA="${NRSUR7DQ4_DATA:-$PWD/lal_data/NRSur7dq4_v1.0.h5}"
PY="${PYTHON:-python}"
cargo build --release --examples
"$PY" scripts/validate_against_lalsimulation.py --n-draws 200 --outdir report/data \
    --approximants TaylorF2 IMRPhenomD IMRPhenomXAS IMRPhenomXHM IMRPhenomXP IMRPhenomXP_NNLO \
    IMRPhenomXP_MSA IMRPhenomXPHM IMRPhenomXPHM_NNLO IMRPhenomXPHM_MSA IMRPhenomXPNR
# IMRPhenomXPNR with the parity fix against lalsimulation 6.2.0 + lalsuite
# !2593, if built (scripts/build_patched_lalsimulation.sh)
if [ -n "${PATCHED_LALSIMULATION:-}" ]; then
    "$PY" scripts/validate_against_lalsimulation.py --n-draws 200 --outdir report/data \
        --approximants IMRPhenomXPNR_parityfix --lalsimulation-library "$PATCHED_LALSIMULATION"
fi
"$PY" scripts/validate_nrsur7dq4.py --n-draws 200 --outdir report/data
"$PY" report/lal_benchmark.py
"$PY" report/nrsur_benchmark.py
"$PY" report/serial_benchmark.py --n-draws 300
"$PY" report/long_benchmark.py --memory
"$PY" report/conversions_benchmark.py
cargo run --quiet --release --example batch_benchmark -- --json report/data/batch_scaling.json
# the Python bindings (installed with `maturin develop --release` in python/)
"$PY" python/benchmarks/benchmark.py --json report/data/python_timings.json
"$PY" report/build_report.py

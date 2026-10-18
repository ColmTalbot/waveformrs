#!/usr/bin/env bash
# Timings of waveformrs and of the container's lalsimulation, summarised
# under the key <label> of report/data/linux_timings_<platform>.json.
# Run inside the container (docker/run.sh [--platform ...] bench).
#
#     scripts/linux/bench.sh [--label NAME] [--serial-draws N] [--skip-lal]
#
# Runs report/lal_benchmark.py (examples/benchmark.rs against
# SimInspiralChooseFDWaveformSequence), report/serial_benchmark.py
# (examples/serial_benchmark.rs) and examples/grid_benchmark.rs on one
# thread, and examples/batch_benchmark.rs on 1 to all threads. Use the
# options of docker/run.sh (e.g. --target-cpu native) to time
# another build of the crate; --skip-lal then times only the crate (single
# calls, grid and batch).
#
# Timings under emulation (linux/amd64 on Apple silicon) are not meaningful
# as absolute numbers.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.."

platform=${WAVEFORMRS_PLATFORM:-linux-$(uname -m)}
label=$platform${WAVEFORMRS_VARIANT:-}
serial_draws=300
skip_lal=
while [[ $# -gt 0 ]]; do
    case "$1" in
    --label) label=$2; shift 2 ;;
    --serial-draws) serial_draws=$2; shift 2 ;;
    --skip-lal) skip_lal=1; shift ;;
    *) echo "unknown argument $1" >&2; exit 2 ;;
    esac
done
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
out=${BENCH_DIR:-/target/bench/$label}
mkdir -p "$out"

cargo build --quiet --release --examples
if [[ -z $skip_lal ]]; then
    python report/lal_benchmark.py --output "$out/timings.json"
    python report/serial_benchmark.py --n-draws "$serial_draws" --output "$out/serial_timings.json"
else
    NRSUR7DQ4_DATA=lal_data/NRSur7dq4_v1.0.h5 cargo run --quiet --release --example benchmark -- --json \
        >"$out/benchmark.json"
fi
cargo run --quiet --release --example grid_benchmark | tee "$out/grid.txt"
# throughput on 1, 2, 4, 8 and all threads (examples/batch_benchmark.rs)
NRSUR7DQ4_DATA=lal_data/NRSur7dq4_v1.0.h5 cargo run --quiet --release --example batch_benchmark \
    | tee "$out/batch.txt"
python scripts/linux/summarise.py bench "$out" --label "$label" \
    --output "report/data/linux_timings_${platform#linux-}.json"

#!/usr/bin/env bash
# Random-draw comparison of waveformrs with the lalsimulation of this
# (Linux) machine, summarised under the key <label> of
# report/data/linux_validation_<platform>.json.
# Run inside the container (docker/run.sh [--platform ...] validate).
#
#     scripts/linux/validate.sh [--n-draws N] [--nrsur-draws N] [--label NAME]
#
# --label names the configuration in the summary (default: the platform;
# with the options of docker/run.sh, e.g. "linux-arm64-fma-off" for --fma off).
# The per-draw records go to $VALIDATION_DIR (default /target/validation/<label>).
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.."

n_draws=200
nrsur_draws=200
platform=${WAVEFORMRS_PLATFORM:-linux-$(uname -m)}
label=$platform${WAVEFORMRS_VARIANT:-}
while [[ $# -gt 0 ]]; do
    case "$1" in
    --n-draws) n_draws=$2; shift 2 ;;
    --nrsur-draws) nrsur_draws=$2; shift 2 ;;
    --label) label=$2; shift 2 ;;
    *) echo "unknown argument $1" >&2; exit 2 ;;
    esac
done
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
outdir=${VALIDATION_DIR:-/target/validation/$label}
mkdir -p "$outdir"

python scripts/validate_against_lalsimulation.py --n-draws "$n_draws" --outdir "$outdir" \
    --approximants TaylorF2 IMRPhenomD IMRPhenomXAS IMRPhenomXHM \
    IMRPhenomXP_NNLO IMRPhenomXP_MSA IMRPhenomXP \
    IMRPhenomXPHM_NNLO IMRPhenomXPHM_MSA IMRPhenomXPHM IMRPhenomXPNR
if [[ $nrsur_draws -gt 0 && -f lal_data/NRSur7dq4_v1.0.h5 ]]; then
    NRSUR7DQ4_DATA=lal_data/NRSur7dq4_v1.0.h5 LAL_DATA_PATH=lal_data \
        python scripts/validate_nrsur7dq4.py --n-draws "$nrsur_draws" --outdir "$outdir"
fi
python scripts/linux/summarise.py validation "$outdir" --label "$label" \
    --output "report/data/linux_validation_${platform#linux-}.json"

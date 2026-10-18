#!/usr/bin/env bash
# Build and use the Linux containers of docker/Dockerfile. See docs/linux.md.
#
#     docker/run.sh [options] <command> [args...]
#
# Options
#   --platform arm64|amd64   container platform (default: the Docker host's;
#                            amd64 runs emulated on Apple silicon)
#   --fma on|off             build with --cfg waveformrs_fma="on"|"off",
#                            overriding the platform's FMA mirroring
#   --target-cpu CPU         build with -C target-cpu=CPU (e.g. native,
#                            x86-64-v3)
#   Builds with options use their own target directory under /target.
#
# Commands
#   build             build the image waveformrs-linux:<platform>
#   test [args]       cargo test --release [args]
#   check             cargo build (no warnings, also --no-default-features),
#                     clippy and fmt for this platform
#   validate [args]   scripts/linux/validate.sh: random draws against the
#                     container's lalsimulation (OMP_NUM_THREADS=1)
#   bench [args]      scripts/linux/bench.sh: timings against lalsimulation
#   reference [args]  scripts/linux/reference.sh: reference data from the
#                     container's lalsimulation, and cargo test against them
#   shell             an interactive shell
#   <anything else>   run that command in the container
#
# The repository is mounted read-write at /work. Build products go to the
# volume waveformrs-target-<platform> (mounted at /target), so they never mix
# with the host's target/. The NRSur7dq4 data directory (default: lal_data/ of
# the checkout or of the main checkout of a git worktree; override with
# WAVEFORMRS_LAL_DATA) is mounted read-only at /work/lal_data.
set -euo pipefail

repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
case "$(docker version --format '{{.Server.Arch}}' 2>/dev/null)" in
amd64 | x86_64) platform=amd64 ;;
*) platform=arm64 ;;
esac
rustflags=${RUSTFLAGS:-}
variant=

while [[ $# -gt 0 ]]; do
    case "$1" in
    --platform)
        platform=${2#linux/}
        shift 2
        ;;
    --fma)
        [[ $2 == on || $2 == off ]] || { echo "--fma on|off" >&2; exit 2; }
        rustflags="$rustflags --cfg waveformrs_fma=\"$2\""
        variant="$variant-fma-$2"
        shift 2
        ;;
    --target-cpu)
        rustflags="$rustflags -C target-cpu=$2"
        variant="$variant-$2"
        shift 2
        ;;
    -h | --help)
        sed -n '2,32p' "${BASH_SOURCE[0]}"
        exit 0
        ;;
    *) break ;;
    esac
done
case "$platform" in
arm64 | aarch64) platform=arm64 ;;
amd64 | x86_64) platform=amd64 ;;
*)
    echo "unknown platform $platform (arm64 or amd64)" >&2
    exit 2
    ;;
esac
image=${WAVEFORMRS_IMAGE:-waveformrs-linux}:$platform
command=${1:-shell}
[[ $# -gt 0 ]] && shift

if [[ $command == build ]]; then
    exec docker build --platform "linux/$platform" -t "$image" "$@" "$repo/docker"
fi

# The NRSur7dq4 data: the checkout's own lal_data/, else the main checkout's.
lal_data=${WAVEFORMRS_LAL_DATA:-}
if [[ -z $lal_data ]]; then
    common=$(git -C "$repo" rev-parse --path-format=absolute --git-common-dir 2>/dev/null || true)
    for candidate in "$repo/lal_data" "${common:+$common/../lal_data}"; do
        if [[ -n $candidate && -f $candidate/NRSur7dq4_v1.0.h5 ]]; then
            lal_data=$(cd "$candidate" && pwd)
            break
        fi
    done
fi
mounts=(-v "$repo:/work"
    -v "waveformrs-target-$platform:/target"
    -v "waveformrs-cargo-registry:/usr/local/cargo/registry")
if [[ -n $lal_data ]]; then
    mounts+=(-v "$lal_data:/work/lal_data:ro")
else
    echo "note: NRSur7dq4_v1.0.h5 not found; the NRSur7dq4 tests will be skipped" >&2
fi

tty=()
[[ -t 0 && -t 1 ]] && tty=(-it)
run() {
    docker run --rm ${tty[@]+"${tty[@]}"} --platform "linux/$platform" "${mounts[@]}" \
        -e WAVEFORMRS_PLATFORM="linux-$platform" \
        -e WAVEFORMRS_VARIANT="$variant" \
        -e RUSTFLAGS="${rustflags# }" \
        -e CARGO_TARGET_DIR="/target/build${variant}" \
        -e PYTHONUNBUFFERED=1 \
        "$image" "$@"
}

case "$command" in
shell) run bash ;;
test) run cargo test --release "$@" ;;
# (clippy 1.90 lints that the toolchain of the reference checks, 1.86, lacks are
# allowed)
check) run bash -c 'set -o pipefail; { cargo build --all-targets && cargo build --all-targets --no-default-features; } 2>&1 | tee /tmp/build.log \
    && ! grep -q "^warning" /tmp/build.log \
    && cargo clippy --all-targets -- -D warnings -A clippy::neg_multiply -A clippy::manual_is_multiple_of \
    && cargo fmt --check' ;;
validate) run bash scripts/linux/validate.sh "$@" ;;
bench) run bash scripts/linux/bench.sh "$@" ;;
reference) run bash scripts/linux/reference.sh "$@" ;;
*) run "$command" "$@" ;;
esac

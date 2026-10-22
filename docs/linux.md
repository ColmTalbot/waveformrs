# waveformrs on Linux

The committed reference data in `tests/data` come from the conda-forge
**osx-arm64** build of lalsimulation 6.2.0. The crate does not copy that
build, or any other, operation by operation: it uses Rust's elementary
functions and fuses multiply-adds wherever the hardware has the
instruction (see "Floating-point decisions" below). This page covers Linux. It describes the Docker containers
that test and benchmark the crate against the conda-forge **linux-aarch64**
and **linux-64** builds of the same lalsimulation, what agrees with what and
how closely, and what was changed for Linux.

All numbers below were measured on an Apple M2 Ultra (16 performance and 8
efficiency cores) with Docker Desktop 29.1.3: a linux/aarch64 VM with 24 CPUs
and 7.6 GB, and linux/amd64 under Rosetta emulation. The data are in
`report/data/linux_*.json`.

## Quick start

```sh
docker/run.sh --platform arm64 build      # once per platform (amd64: emulated)
docker/run.sh --platform arm64 test       # cargo test --release, committed data
docker/run.sh --platform arm64 reference  # tests/data regenerated with the container's LAL, then cargo test
docker/run.sh --platform arm64 validate   # random draws against the container's LAL
docker/run.sh --platform arm64 bench      # timings against the container's LAL
docker/run.sh --platform arm64 check      # build warnings, clippy, rustfmt
docker/run.sh --platform arm64 shell
scripts/linux/test_all.sh                 # the whole test matrix on both platforms
```

Options for every command: `--fma on|off` builds with
`--cfg waveformrs_fma="on"|"off"`, and `--target-cpu CPU` builds with
`-C target-cpu=CPU`. Each variant gets its own target directory.

* **Image** (`docker/Dockerfile`, build context `docker/`): `rust:1.90.0`
  (with clippy and rustfmt), and micromamba 1.5.10 with
  `docker/environment.yml`: python 3.12, python-lalsimulation 6.2.0, liblal
  7.7.0, GSL 2.8, OpenBLAS 0.3.34, numpy, scipy, h5py and matplotlib-base.
  These are the versions of the osx-arm64 reference environment. The exact
  package list is in `/opt/conda/lal-explicit.txt` in the image.
  `docker/conda-install.sh` solves first. It then downloads the ~175
  packages one at a time, at most one per second, into a BuildKit cache, so
  a rebuild fetches nothing again.
* **Mounts** (`docker/run.sh`): the repository is mounted at `/work`. Build
  products go to the volume `waveformrs-target-<platform>` at `/target` and
  never mix with the host's `target/`. The cargo registry is a shared volume.
  `lal_data/` is mounted read-only at `/work/lal_data`. It is taken from the
  checkout, or from the main checkout of a git worktree, and
  `WAVEFORMRS_LAL_DATA` overrides it. `OMP_NUM_THREADS=1` and
  `OPENBLAS_NUM_THREADS=1` are set.
* **Scripts** (`scripts/linux/`): `validate.sh`, `bench.sh`, `reference.sh`
  and `test_all.sh`. `summarise.py` writes `linux_validation_*.json` and
  `linux_timings_*.json`. `compare_reference.py` compares two directories of
  reference data and writes `linux_reference_comparison.json`.

## Test matrix

`scripts/linux/test_all.sh` (results at 4520bfd):

| | macOS arm64 (host) | linux/arm64 | linux/amd64 (emulated) |
|---|---|---|---|
| `cargo test --release`, committed data | pass | pass | pass |
| against the container's own LAL (`reference`) | — | pass | pass |
| `--no-default-features` (lib and tests) | pass | pass | pass |
| build warnings / clippy / rustfmt | clean | clean¹ | clean¹ |
| `-C target-cpu=x86-64-v3` | — | — | not testable² |

¹ `docker/run.sh check` uses clippy 1.90. It allows `neg_multiply` and
`manual_is_multiple_of`, two lints that the host toolchain (1.86) does not
have.

² Any build with AVX enabled gives wrong results under Rosetta. That covers
`-C target-feature=+avx`, which `+fma` and `x86-64-v3` imply. With `+avx`
alone, the crate's arithmetic is unchanged (no FMA), and still the
IMRPhenomXAS duration test is off by 6%. An `+sse4.1` build
passes. These builds need real x86-64 hardware to be validated.

The doctests of `src/conversions` fail with `--no-default-features` on every
platform, because they use items behind the `fft` feature. `test_all.sh`
therefore runs only the library and integration tests in that
configuration.

### How the tests decide their tolerances

The tolerances are the same on every platform. `tests/common::tolerances(suite,
strict)` looks the suite up in `AMPLIFYING_SUITE_TOLERANCES`, which covers
only the suites that amplify last-bit differences:

* the IMRPhenomXHM (3,±2) mode, the l = 3 modes of XPHM, and every XHM/XPHM
  polarisation (relative 1e-3 to 2e-1; the (3,±2) modes differ more on
  the Ubuntu runners of GitHub Actions than in these containers);
* the MSA angles;
* the SpinTaylor angles of IMRPhenomXP_NRTidalv2 and v3 (relative 1e-6);
* the SpinTaylor orbit and angles;
* the long-waveform and conversion tests of these models.

Each value is about ten times the largest difference seen on macOS, Linux
arm64 and Linux amd64, against either the osx-arm64 or the Linux LAL. Those
differences are the size of the differences between LAL's own builds (next
section). Every other suite has its strict tolerance.

Unit tests with hard-coded LAL values use `src/reference.rs`: a few ulp (up
to 1e4 ulp for the QNM fits, whose denominators amplify rounding).

`WAVEFORMRS_REFERENCE_DATA=<dir>` makes the tests read the reference suites
from `<dir>` instead of `tests/data`. `scripts/linux/reference.sh` generates
`<dir>` with the container's LAL. The long-waveform and conversion data
always come from `tests/data`, and so do suites missing from `<dir>`, such
as the `imrphenomxpnr_parityfix*` suites, which need a lalsimulation built
with lalsuite !2593 (`scripts/build_patched_lalsimulation.sh`).

`tests/batch.rs` (bit for bit across thread counts) passes unchanged on
every platform.

## How far lalsimulation itself moves between builds

`compare_reference.py` on the reference suites, LAL against LAL. Each entry
is the largest relative error, or the absolute error for angles:

| suite | linux-aarch64 vs osx-arm64 | linux-64 vs osx-arm64 | linux-64 vs linux-aarch64 |
|---|---|---|---|
| IMRPhenomXHM modes | 1.1e-3 | 2.2e-3 | 2.1e-3 |
| ChooseFDModes (XHM, XPHM) | 7.0e-4 | 1.2e-3 | 1.2e-3 |
| IMRPhenomXPHM (SpinTaylor) | 2.2e-5 | 9.0e-6 | 3.0e-5 |
| IMRPhenomXP (MSA) | 3.6e-9 | 3.2e-9 | 3.3e-9 |
| SpinTaylor orbit | 7.5e-10 | 1.3e-9 | 5.4e-10 |
| SpinTaylor angles (rad) | 9.0e-11 | 1.8e-11 | 9.2e-11 |
| TaylorF2 / IMRPhenomXAS | 7e-12 / 2e-12 | 7e-12 / 3e-12 | 7e-12 / 2e-12 |
| NRSur7dq4 | 3.1e-14 | 7.0e-13 | 7.0e-13 |

LAL's own Linux builds therefore differ from its macOS build by up to 0.2%
in the (3,2) mode of IMRPhenomXHM. That mode's intermediate phase takes a
finite-difference second derivative with a step of 1e-7. The crate uses the
analytic derivatives instead (CHANGELOG), so its own builds agree to ~1e-12
in IMRPhenomXHM, and it differs from each LAL build by about as much as
the LAL builds differ from each other.

## Floating-point decisions

The crate used to copy LAL's arithmetic per platform: clang's or GCC's
choice of fused multiply-adds, a port of glibc's `cbrt` on Linux, and the C
library's `cexp`, `cpow`, `csqrt` and Darwin's `__sincos_stret` through FFI.
With the causes of the differences understood (the previous section), it now
uses the best available operations instead, and the tests allow for the
differences on every platform.

* **Fused multiply-adds** (`src/spintaylor/fp.rs`): `fma`, `fms` and `fnms`
  fuse wherever the target has the instruction (aarch64, and x86-64 with the
  `fma` target feature, e.g. `-C target-cpu=native` on Haswell and later).
  A fused product and sum is at least as accurate as a separately rounded
  one. Baseline x86-64 rounds them separately, as `f64::mul_add` would be a
  slow software call there. The `--cfg waveformrs_fma` override and
  `docker/run.sh --fma` are gone.
* **Cube roots**: `f64::cbrt` everywhere, and `cbrt` instead of
  `powf(1.0 / 3.0)`. On linux-gnu this is the correctly rounded `cbrt` of
  Rust's `compiler-builtins`; on macOS it is Apple's. glibc's, which LAL
  calls on Linux, differs from the correctly rounded value in the last bit
  for 57% of arguments.
* **Complex functions** (`src/complex.rs`): `exp(i phi)` from `sin_cos`, the
  principal square root by C99's cancellation-free algorithm, and the
  principal cube root as `cbrt(|z|) exp(i arg(z) / 3)` instead of
  `exp(log(z) / 3)`. None of them goes through FFI any more.
* **NRSur7dq4**: `f64::sin_cos`, `sin` and `cos` rather than the calls that
  matched clang's combination into `__sincos_stret` on Darwin.

## Accuracy against the container's lalsimulation

`scripts/validate_against_lalsimulation.py` and
`scripts/validate_nrsur7dq4.py` (`docker/run.sh validate`), with
OMP_NUM_THREADS=1. The arm64 runs use 200 draws per model and the emulated
amd64 runs use 100 (NRSur7dq4: 200 and 50). The table gives the largest and
median of max|h_rs − h_lal| / max|h_lal|. Data:
`report/data/linux_validation_{arm64,amd64}.json`. On osx-arm64 every
model agrees with LAL to ~1e-14 or better, except TaylorF2, IMRPhenomD and
IMRPhenomXAS, which are at ~1e-11 everywhere because they evaluate powers of
the frequency differently from LAL.

| model | linux/arm64 max | median | linux/amd64 max | median |
|---|---|---|---|---|
| TaylorF2 | 2.9e-11 | 8.0e-13 | 2.9e-11 | 8.8e-13 |
| IMRPhenomD | 3.6e-11 | 4.5e-13 | 2.9e-11 | 5.4e-13 |
| IMRPhenomXAS | 4.6e-11 | 5.7e-13 | 4.4e-11 | 4.0e-13 |
| IMRPhenomXHM | 9.9e-6 | 1.9e-7 | 2.7e-5 | 2.5e-7 |
| IMRPhenomXP, NNLO | 3.0e-12 | 5.5e-14 | 5.6e-13 | 3.7e-14 |
| IMRPhenomXP, MSA | 5.8e-8 | 1.2e-13 | 5.3e-9 | 7.9e-14 |
| IMRPhenomXP, SpinTaylor | 1.1e-5 | 2.0e-12 | 5.6e-6 | 1.8e-12 |
| IMRPhenomXPHM, NNLO | 3.4e-5 | 3.6e-7 | 7.6e-5 | 3.0e-7 |
| IMRPhenomXPHM, MSA | 7.5e-5 | 4.1e-7 | 1.1e-4 | 4.3e-7 |
| IMRPhenomXPHM, SpinTaylor | 1.1e-4 | 6.3e-7 | 6.9e-5 | 6.1e-7 |
| IMRPhenomXPNR | 5.1e-5 | 4.9e-7 | 1.5e-5 | 4.6e-7 |
| NRSur7dq4 | 1.5e-13 | 1.3e-14 | 1.0e-15 | 5.5e-16 |

The IMRPhenomXPNR row was measured when the crate could still reproduce
lalsimulation 6.2.0's IMRPhenomXPNR (without the parity fix of lalsuite
!2593). The crate's IMRPhenomXPNR now always includes the fix, so
`validate.sh` no longer compares it: the containers have no lalsimulation
built with the fix.

All mismatches are below 4e-8. Before the crate stopped copying LAL's
arithmetic, linux-64 agreed with LAL to ~1e-15 for most models (it then
mirrored that build's unfused arithmetic and glibc's `cbrt`), and the
higher-mode models on linux-aarch64 to ~1e-8 (median). The higher-mode and
SpinTaylor models now differ from both Linux builds of LAL by ~1e-7
(median) to ~1e-4 (max), the size of LAL's own differences between its
builds, and are closer than before to the committed osx-arm64 data (e.g.
IMRPhenomXPHM 7e-4 -> 4e-4, IMRPhenomXP_NRTidal 6e-8 -> 5e-9 on linux/arm64).

The crate's own builds differ from each other (`examples/lal_validation` on
the same 20 binaries per model on macOS arm64, Linux arm64 and Linux
amd64) by at most ~1e-12 for TaylorF2, IMRPhenomD, IMRPhenomXAS and
IMRPhenomXP with NNLO angles (TaylorF2 and IMRPhenomD are identical on the
two Linux builds), ~1e-11 to ~1e-9 for IMRPhenomXP with SpinTaylor and MSA
angles, and ~1e-6 (median) to ~2e-5 (max) for the higher-mode models.

## Performance on Linux arm64

`docker/run.sh --platform arm64 bench`, single thread, on the same M2 Ultra
as the macOS numbers. Data: `report/data/linux_timings_arm64.json`; the
macOS columns are from `report/data/timings.json`.

**One call** (`examples/benchmark.rs` against
`SimInspiralChooseFDWaveformSequence`, 5,921 frequencies, best of 10
batches, µs):

| model | Rust, Linux | LAL, Linux | speed-up | Rust, macOS | LAL, macOS |
|---|---|---|---|---|---|
| TaylorF2 | 164 | 215 | 1.3× | 110 | 172 |
| IMRPhenomD | 177 | 295 | 1.7× | 141 | 306 |
| IMRPhenomXAS | 177 | 197 | 1.1× | 103 | 208 |
| IMRPhenomXHM | 1,161 | 1,705 | 1.5× | 1,035 | 1,734 |
| IMRPhenomXP (102) | 404 | 588 | 1.5× | 350 | 648 |
| IMRPhenomXP (223) | 1,231 | 1,540 | 1.3× | 1,130 | 1,549 |
| IMRPhenomXP (320) | 864 | 2,587 | 3.0× | 790 | 2,916 |
| IMRPhenomXPHM (102) | 2,028 | 4,035 | 2.0× | 2,018 | 4,271 |
| IMRPhenomXPHM (223) | 5,274 | 8,805 | 1.7× | 5,003 | 8,664 |
| IMRPhenomXPHM (320) | 3,525 | 7,590 | 2.2× | 3,526 | 8,521 |
| NRSur7dq4 | 1,105 | 5,593 | 5.1× | 1,011 | 5,249 |

**Serial use** (`report/serial_benchmark.py`, 300 random binaries,
8,033 frequencies, mean µs per waveform): the grid paths are 1.9× (XHM) to
7.3× (TaylorF2) faster than LAL. In the per-frequency paths, TaylorF2 is
1.4×, IMRPhenomD 1.7×, IMRPhenomXAS 1.05×, IMRPhenomXHM 1.3×, IMRPhenomXP
1.2–2.2×, IMRPhenomXPHM 1.7–1.9× and NRSur7dq4 5.0× faster.

**Threads** (`examples/batch_benchmark.rs`, 24 vCPUs): 1 to 24 threads give
13–15× the single-thread throughput for every model; 8 threads give
6.5–7.3×. IMRPhenomXAS runs at 150,000 waveforms per second on 24 threads,
and NRSur7dq4 at 12,800.

**Where Linux differs.** LAL takes about as long on both systems (within
10%; TaylorF2 25% longer). So do the crate's higher-mode, precessing and
NRSur7dq4 models (within 15%). Only the per-frequency
paths of TaylorF2, IMRPhenomD and IMRPhenomXAS are 1.3–1.7× slower than on
macOS. `perf` (in a privileged container) shows why: in IMRPhenomXAS about
half of the time is in glibc's libm, mostly
`sincos` (16%), `cbrt` with its `frexp`/`scalbn` (15%), `atan` (5%) and
`log`/`exp` (6%). glibc's `sin`/`cos` are ~1.5× slower than Apple's and its
`cbrt` ~3×, while its `exp`, `log` and `pow` are faster, which is why LAL,
built on `pow`, does not slow down. The crate calls the corresponding
functions, now Rust's correctly rounded `cbrt` (an out-of-line call) rather
than an inlined port of glibc's: these per-frequency paths are 6–18% slower
for it than with the port (`examples/benchmark.rs`, best of two runs:
TaylorF2 166 → 185 µs, IMRPhenomD 188 → 223 µs, IMRPhenomXAS 186 → 198
µs), while the other models are within ±6%. For many waveforms on the same
frequencies, the grid paths avoid libm for the phase and are unaffected
(TaylorF2, IMRPhenomD and IMRPhenomXAS on grids: 42, 110 and 118 µs per
waveform in `examples/grid_benchmark.rs`, as on macOS).

No allocator effects showed up in the profiles of these models:
`malloc`/`free` are not among the top symbols. The other models were not
profiled. `-C target-cpu=native` made no measurable difference on arm64,
where NEON is baseline.

## x86-64

* `f64::mul_add` without the `fma` target feature calls compiler-builtins'
  software `fma`, and `round_ties_even` without SSE4.1 calls `roundeven`.
  Both are slow, and they kept `src/fastmath.rs` from vectorising.
  `fastmath` now has variants for targets without FMA, which a test checks
  on every target against the C library (within 1 ulp):
  * the Cody–Waite reduction splits the quadrant number so that every
    product with the 33-bit parts of π/2 is exact;
  * the polynomials multiply and add separately;
  * rounding uses the 1.5·2^52 shift.

  The QNM fit sums go through `fp::fma`, so they are unfused on baseline
  x86-64.
* `-C target-cpu=x86-64-v3` (`docker/run.sh --platform amd64 --target-cpu
  x86-64-v3`) is an opt-in. It turns on hardware FMA (in `fp.rs` as well
  as in `fastmath`) and AVX2 for the vectorised grid loops. The FMA
  variants pass their tests everywhere, because aarch64 uses them. Whole
  AVX builds can't be validated under Rosetta, though (note ² above).
* The linux/amd64 container runs under Rosetta. Its timings are not
  meaningful as absolute numbers, so this page reports correctness only
  there. Rosetta reports AVX2 and FMA, so glibc and OpenBLAS pick their
  dispatch variants as for that CPU (OpenBLAS picks Nehalem).

## Open ends

* **Faster per-frequency paths on Linux.** glibc's `sincos` and `atan` are
  slower than Apple's. A vectorised `sin_cos`, which the grid paths already
  use, would speed up the per-frequency paths of TaylorF2, D and XAS, at
  an accuracy of a few ulp instead of under one.
* **Real x86-64 hardware.** Nothing here was measured on it: not the
  timings, not OpenBLAS's Haswell+ kernels, and not the AVX/FMA builds
  that Rosetta miscomputes.

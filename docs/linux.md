# waveformrs on Linux

The crate reproduces the conda-forge **osx-arm64** build of lalsimulation
6.2.0 to rounding, for most models bit for bit, and that build produced the
committed reference data
in `tests/data`. This page covers Linux. It describes the Docker containers
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
| `cargo test --release`, committed data | pass, strict | pass, cross-build tolerances | pass, cross-build tolerances |
| against the container's own LAL (`reference`) | — | pass, cross-build tolerances | pass, cross-build tolerances |
| `--no-default-features` (lib and tests) | pass | pass | pass |
| build warnings / clippy / rustfmt | clean | clean¹ | clean¹ |
| `-C target-cpu=x86-64-v3` | — | — | not testable² |

¹ `docker/run.sh check` uses clippy 1.90. It allows `neg_multiply` and
`manual_is_multiple_of`, two lints that the host toolchain (1.86) does not
have.

² Any build with AVX enabled gives wrong results under Rosetta. That covers
`-C target-feature=+avx`, which `+fma` and `x86-64-v3` imply. With `+avx`
alone, the crate's arithmetic is unchanged (no FMA, unfused `fp.rs`), and
still the IMRPhenomXAS duration test is off by 6%. An `+sse4.1` build
passes. These builds need real x86-64 hardware to be validated.

The doctests of `src/conversions` fail with `--no-default-features` on every
platform, because they use items behind the `fft` feature. `test_all.sh`
therefore runs only the library and integration tests in that
configuration.

### How the tests decide their tolerances

The crate mirrors one build of LAL operation by operation: osx-arm64 (clang
contractions, Apple's libm). `tests/common::mirrored_build()` is true only on
macOS arm64, and the strict tolerances (and the bit-for-bit checks) apply
only there.

Everywhere else, `tests/common::tolerances(suite, strict)` looks the suite up
in `CROSS_BUILD_TOLERANCES`. That table covers only the suites that amplify
last-bit differences:

* the IMRPhenomXHM (3,±2) mode, the l = 3 modes of XPHM, and every XHM/XPHM
  polarisation (relative 1e-3 to 2e-1; the (3,±2) modes differ more on
  the Ubuntu runners of GitHub Actions than in these containers);
* the MSA angles;
* the SpinTaylor angles of IMRPhenomXP_NRTidalv2 and v3 (relative 1e-6);
* the SpinTaylor orbit and angles;
* the long-waveform and conversion tests of these models.

Each value is about ten times the largest difference seen on Linux arm64 and
amd64, against either the osx-arm64 or the Linux LAL. Those differences are
the size of the differences between LAL's own builds (next section). Every
other suite keeps its strict tolerance on every platform.

Unit tests with hard-coded LAL values use `src/reference.rs`: bit for bit on
macOS arm64, a few ulp elsewhere (up to 1e4 ulp for the QNM fits, whose
denominators amplify rounding).

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
finite-difference second derivative with a step of 1e-7.

## Floating-point decisions

### Fused multiply-adds (`src/spintaylor/fp.rs`)

`FUSED_MULTIPLY_ADD` depends only on the architecture. The choice was made
from the disassembly of the conda-forge binaries and from the random-draw
validation:

* **aarch64, macOS: fused.** Clang (`-ffp-contract=on`) builds LAL and GSL,
  and the crate mirrors it bit for bit.
* **aarch64, Linux: fused.** GCC 14.4 builds LAL and GSL with
  `-ffp-contract=fast`, which fuses every product that clang fuses inside an
  expression, with the same choice of operand. It also fuses products and
  sums that sit in different statements, so `liblalsimulation.so` has about
  11% more fused instructions than the osx-arm64 dylib (46,824 against
  42,090). Mirroring clang's subset is not bit for bit, but it is much
  closer than not fusing (table below). OpenBLAS's neoversen1 kernels fuse
  exactly as mirrored, so the GSL LU solves are bit-identical.
* **x86-64: not fused.** The conda-forge linux-64 builds use
  `-march=nocona`, and `liblalsimulation.so`, `liblal.so` and `libgsl.so`
  contain no FMA instruction. The IMRPhenomX remnant and PN fits and the
  GSL LU solves agree bit for bit. The LU solves were checked against the
  Nehalem kernels that OpenBLAS selects under Rosetta. On Haswell and later,
  OpenBLAS's kernels fuse and the LU solves agree only to rounding.

`--cfg waveformrs_fma="on"|"off"` (`docker/run.sh --fma on|off`) overrides
the choice, for example for a LAL built with clang on Linux or with
`-march=haswell`.

On Linux arm64, median relative errors against the container's LAL over the
same draws:

| | FMA mirrored (default) | not fused |
|---|---|---|
| IMRPhenomXHM | 9.7e-9 | 4.5e-7 |
| IMRPhenomXPHM, NNLO / MSA / SpinTaylor | 1.4e-8 / 3.9e-8 / 3.7e-7 | 5.1e-7 / 4.8e-7 / 6.1e-7 |
| NRSur7dq4 | 1.3e-14 | 9.5e-14 |

### glibc's `cbrt` (`src/cmath.rs`)

On linux-gnu, Rust links its `compiler-builtins` ahead of libm.
`compiler-builtins` defines a correctly rounded `cbrt`, so `f64::cbrt` is
not the glibc `cbrt` that LAL calls. The two differ in the last bit for 57%
of arguments. That was the largest source of differences on Linux.

The crate now uses `cmath::LibmCbrt::libm_cbrt` wherever LAL calls `cbrt`.
On linux-gnu this is an inlined port of glibc's algorithm, with the
polynomial fused on aarch64. It is checked bit for bit against glibc's own
`cbrt`, through its versioned symbol, on 2 million arguments on both
architectures. It is also 10–20% faster than the libm call in the
per-frequency paths of TaylorF2, IMRPhenomD and IMRPhenomXAS. On linux-64,
LAL's IMRPhenomXHM modes then agree bit for bit, and IMRPhenomXP to ~4e-16.

## Accuracy against the container's lalsimulation

`scripts/validate_against_lalsimulation.py` and
`scripts/validate_nrsur7dq4.py` (`docker/run.sh validate`), with
OMP_NUM_THREADS=1. The arm64 runs use 200 draws per model and the emulated
amd64 runs use 100 (NRSur7dq4: 200 and 50). The table gives the largest and
median of max|h_rs − h_lal| / max|h_lal|. Data:
`report/data/linux_validation_{arm64,amd64}.json`. For comparison,
osx-arm64 agrees bit for bit or to ~1e-15 for all models but TaylorF2,
IMRPhenomD and IMRPhenomXAS, which are at ~1e-11 everywhere because they
evaluate powers of the frequency differently from LAL.

| model | linux/arm64 max | median | linux/amd64 max | median |
|---|---|---|---|---|
| TaylorF2 | 4.9e-11 | 7.9e-13 | 6.2e-11 | 9.8e-13 |
| IMRPhenomD | 4.1e-11 | 6.4e-13 | 5.7e-11 | 6.1e-13 |
| IMRPhenomXAS | 6.2e-11 | 1.0e-12 | 1.0e-10 | 7.0e-13 |
| IMRPhenomXHM | 1.6e-5 | 9.7e-9 | 3.8e-7 | 0 (56/100 identical) |
| IMRPhenomXP, NNLO | 2.9e-12 | 3.5e-14 | 1.6e-15 | 3.6e-16 |
| IMRPhenomXP, MSA | 9.4e-9 | 6.8e-14 | 1.8e-15 | 3.7e-16 |
| IMRPhenomXP, SpinTaylor | 9.3e-7 | 1.1e-12 | 2.1e-15 | 3.7e-16 |
| IMRPhenomXPHM, NNLO | 8.8e-6 | 1.4e-8 | 3.4e-7 | 3.8e-16 |
| IMRPhenomXPHM, MSA | 2.0e-5 | 3.9e-8 | 5.5e-7 | 3.9e-16 |
| IMRPhenomXPHM, SpinTaylor | 1.3e-4 | 3.7e-7 | 4.8e-7 | 3.8e-16 |
| NRSur7dq4 | 1.5e-13 | 1.3e-14 | 1.0e-15 | 5.5e-16 |

All mismatches are below 1e-8 on arm64 and below 1e-12 on amd64.

* **linux-64:** apart from the frequency powers of TaylorF2, D and XAS, the
  crate agrees with LAL at the rounding level, as it does on macOS. The
  remaining (3,2)-mode differences (up to 4e-7) come from a few last-bit
  differences that have not been traced.
* **linux-aarch64:** the cross-statement contractions of GCC are not
  mirrored. The higher-mode models agree to ~1e-8 typically and to 1e-4 at
  worst, which is the size of LAL's own difference between its builds.

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
built on `pow`, does not slow down. The crate calls the functions LAL
calls, so that it rounds as LAL does. Two changes help:

* the inlined glibc `cbrt` (above) cut these per-frequency paths by 10–20%
  (TaylorF2 185 → 164 µs, IMRPhenomD 219 → 177 µs, IMRPhenomXAS 205 → 177
  µs);
* for many waveforms on the same frequencies, the grid paths avoid libm
  for the phase (TaylorF2, IMRPhenomD and IMRPhenomXAS on grids: 42, 119
  and 118 µs per waveform in `examples/grid_benchmark.rs`, as on macOS).

No allocator effects showed up in the profiles of these models:
`malloc`/`free` are not among the top symbols. The other models were not
profiled. `-C target-cpu=native` made no measurable difference on arm64,
where NEON is baseline; this was measured before the `cbrt` change.

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

  The QNM fit sums now go through `fp::fma`, so they are unfused on x86-64
  like LAL's.
* `-C target-cpu=x86-64-v3` (`docker/run.sh --platform amd64 --target-cpu
  x86-64-v3`) is an opt-in. It turns on hardware FMA and AVX2 for
  `fastmath`'s vectorised grid loops, while `fp.rs` and `cmath.rs` still
  mirror LAL's and glibc's unfused arithmetic. The FMA variants of
  `fastmath` pass their tests everywhere, because aarch64 uses them. Whole
  AVX builds can't be validated under Rosetta, though (note ² above).
* The linux/amd64 container runs under Rosetta. Its timings are not
  meaningful as absolute numbers, so this page reports correctness only
  there. Rosetta reports AVX2 and FMA, so glibc and OpenBLAS pick their
  dispatch variants as for that CPU (OpenBLAS picks Nehalem).

## Open ends

* **Bit-for-bit agreement with GCC builds.** On linux-aarch64 this would
  need a mirror of GCC's cross-statement contractions, a second `fp.rs`
  contraction rule derived from the disassembly. On linux-64 a few
  differences remain to be traced: up to 4e-7 in the (3,2) mode, and
  SpinTaylor angles that differ by 1e-6 rad after 900 s.
* **Faster per-frequency paths on Linux.** glibc's `sincos` and `atan` are
  slower than Apple's. The per-frequency paths of TaylorF2, D and XAS call
  them as LAL does. A vectorised `sin_cos`, which the grid paths already
  use, would speed them up at the cost of last-bit differences from LAL.
* **Real x86-64 hardware.** Nothing here was measured on it: not the
  timings, not OpenBLAS's Haswell+ kernels, and not the AVX/FMA builds
  that Rosetta miscomputes.

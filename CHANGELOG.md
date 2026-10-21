# Changelog

All notable changes to `waveformrs` are recorded here.

## Unreleased

### License

- **waveformrs is licensed under the GNU General Public License version 3
  or later** (`LICENSE`), as the code it is translated from requires:
  LALSuite (GPL-2.0-or-later) and GSL (GPL-3.0-or-later), with parts of
  glibc (LGPL-2.1-or-later) and fdlibm. `NOTICE.md` lists the original files
  and their authors; the library sources carry SPDX headers, and the crates
  and the Python package declare `GPL-3.0-or-later`. The crate had no
  license before.

### Breaking changes

- **`phi_c` is the orbital reference phase in every `Waveform` method.**
  `Waveform::phase`, `waveform_single_frequency`, `waveform`,
  `waveform_modes`, `waveform_on_grid` and `waveform_modes_on_grid` now take
  the orbital reference phase, LALSimulation's `phi_ref`, for TaylorF2,
  IMRPhenomD and IMRPhenomXAS alike: the (2,2) mode is shifted by
  `2 phi_c`. Previously TaylorF2 used the orbital phase while IMRPhenomD and
  IMRPhenomXAS used the gravitational-wave phase, so reproducing
  `SimInspiralChooseFDWaveform(phi_ref, ...)` needed `phi_c = 2 phi_ref` for
  those two models.

  To migrate, pass `phi_ref` where you passed `2.0 * phi_ref` to IMRPhenomD
  or IMRPhenomXAS (or halve the value you passed). TaylorF2 is unchanged.
  The waveforms are bit for bit those of the old calls with twice the
  phase, so the agreement with LALSimulation is unchanged.

### Added

- **A common model interface**, `waveformrs::model` (re-exported from the
  crate root):
  - `Binary`, the parameters of `SimInspiralChooseFDWaveform` (component
    masses in solar masses, spin vectors at `f_ref` with `L` along `z`,
    distance in Mpc, inclination, orbital `phi_ref`, `f_ref` in Hz).
    `imrphenomx::PrecessingBinary` (and `nrsur::PrecessingBinary`) is now a
    type alias of `Binary`, so existing code is unaffected.
  - `FrequencyDomainModel`, implemented by `TaylorF2`, `IMRPhenomD`,
    `IMRPhenomXAS`, `IMRPhenomXHM`, `IMRPhenomXP` and `IMRPhenomXPHM`:
    `from_binary` (with LAL's checks, including the rejection of in-plane
    spins by aligned-spin approximants and `f_ref = 0`), `polarizations`
    (`SimInspiralChooseFDWaveformSequence`), `polarizations_into`
    (caller-owned buffers) and `polarizations_on_grid` (the fast grid
    paths).
  - `TimeDomainModel`, implemented by `NRSur7dq4` (polarizations and modes,
    with a reusable workspace).
  - `Approximant`, with LAL's names (`Display`, `FromStr`), and
    `choose_fd_model` / `choose_td_model`, which build a boxed model from an
    approximant, a `Binary` and `WaveformOptions` (Euler angles, mode
    array).
  - `waveformrs::Error`, wrapping `IMRPhenomXError` and `NRSurError`.
- **Frequency-domain modes** in the conventions of
  `SimInspiralChooseFDModes`: `FrequencyDomainModel::modes` returns
  `FrequencyModes`, the modes `h_lm` at `+f` and `-f` with the frame in
  which to combine them (`FrequencyModes::polarizations` rebuilds the
  polarizations). IMRPhenomXHM returns its modes (`m < 0` at positive,
  `m > 0` at negative frequencies, including `phi_ref`); IMRPhenomXPHM
  returns the J-frame modes of `XLALSimIMRPhenomXPHMOneMode`, including
  LAL's per-mode range of the SpinTaylor angles; TaylorF2, IMRPhenomD and
  IMRPhenomXAS return their (2,+-2) modes and IMRPhenomXP its J-frame (2,m)
  modes, consistent with their polarizations. `choose_fd_modes` reproduces
  `SimInspiralChooseFDModes` on uniform grids for IMRPhenomXHM and
  IMRPhenomXPHM. Tested against new reference data
  (`tests/data/choose_fd_modes.json`, `tests/modes.rs`).
- **A switch for the automatic changes of the precession prescription** of
  IMRPhenomXP and IMRPhenomXPHM: `with_version_fallback(false)` (and
  `WaveformOptions::version_fallback`). By default the models still change
  the Euler angles as LALSimulation does (SpinTaylor to NNLO for in-plane
  spins below 1e-7, SpinTaylor to MSA when the SpinTaylor angles cannot be
  computed, MSA to NNLO when the MSA system cannot be initialised); with the
  switch off they return the new
  `IMRPhenomXError::PrecessionVersionUnavailable` instead (polarizations,
  grid and uniform paths, and modes).
- **Batched evaluation**, `waveformrs::batch::Batch`: many binaries (built
  with any `from_binary`-like function) or prepared models on one
  `FrequencyGrid` (`polarizations_on_grid`, `models_on_grid`, into
  caller-provided row-major buffers), at arbitrary frequencies
  (`polarizations`, returning a `Vec`), and NRSur7dq4 in the time domain
  (`time_domain_polarizations`, `time_domain_models`, one workspace per
  thread). `Batch::new(threads)` splits the work into chunks taken by
  `std::thread::scope` threads (no new dependency); one thread is exactly
  a serial loop and the results do not depend on the number of threads
  (`tests/batch.rs`). `examples/batch_benchmark.rs` reports the throughput
  for 1, 2, 4, 8 and all threads.
- **IMRPhenomXPNR** (`imrphenomx::IMRPhenomXPNR`, `Approximant::IMRPhenomXPNR`):
  IMRPhenomXPHM with LAL's PNR flags (NR-tuned precession angles with the
  frequency map of the higher modes, the tuned co-precessing (2,2) mode,
  its antisymmetric part, SpinTaylor angles 330, final spin 7 and the
  122019 release of the higher modes), including LAL's own switches
  (plain XPHM for in-plane spins below 1e-7, NNLO angles below 1e-6, MSA
  angles when the SpinTaylor integration fails, ...). Polarizations on
  frequency sequences and uniform grids, the J-frame modes of
  `SimInspiralChooseFDModes` (`choose_fd_modes`) and of the single-mode
  functions (`FrequencyDomainModel::modes`), `simulate_td`, and
  `with_version_fallback(false)`. Tested against LALSimulation 6.2.0
  (`tests/xpnr.rs`, worst relative difference ~1e-15).
- **The IMRPhenomXPNR parity fix of lalsuite master** (merge request
  !2593, merged after lalsimulation 6.2.1), on by default:
  `IMRPhenomXPNR::with_parity_fix(false)` and
  `WaveformOptions::pnr_parity_fix = false` (`with_pnr_parity_fix`, also
  set by `WaveformOptions::lalsimulation_defaults()`) reproduce
  lalsimulation 6.2.0 and 6.2.1. The phase of the antisymmetric (2,2)
  co-precessing mode is referred to `-alpha(f_ref) - pi / 2` and the
  azimuth of `J` instead of `alpha_offset` and the polarization angle,
  which depend on the line of sight, and the spins at the end of the
  SpinTaylor integration are rotated with the components of `L`. The
  waveforms change by up to a few per cent; the mirror image of a binary
  now gives `(h_+, -h_x)` from the mirrored direction, as reflection
  symmetry requires (violated at 2e-3 to 1e-1 before). Tested against
  lalsimulation 6.2.0 built with the patch
  (`scripts/build_patched_lalsimulation.sh`, suites
  `imrphenomxpnr_parityfix*`, ~1e-15) and LAL's own unit test.
- **A grid path for IMRPhenomXPNR** (`IMRPhenomXPNR::polarizations_on_grid`
  and `FrequencyDomainModel::polarizations_on_grid`, previously the
  per-frequency path): the set-up of the per-frequency path, then for each
  co-precessing mode its PNR angles (with the spline coefficients computed
  once per waveform), the twisting-up factors, the mode from the powers of
  the grid and the antisymmetric (2,2) mode in passes over the grid.
  Agrees with the per-frequency path to ~1e-13 (`tests/grid.rs`); 2.2 times
  faster on 8033 frequencies (`examples/grid_benchmark.rs`). The moving
  average of the antisymmetric amplitude ratio no longer recomputes windows
  whose inputs repeat those of the previous window (the same values, on
  both paths).
- **IMRPhenomXHM release 122019** (`with_release_version(122019)`), the
  original higher-mode amplitudes used by IMRPhenomXPNR.
- **Python bindings** (`python/`, the package `waveformrs` built with
  maturin; PyO3 0.29 and rust-numpy 0.29, stable ABI for CPython >= 3.9):
  the approximants by LAL's names, `WaveformOptions` as keyword arguments
  or an object, `polarizations` (`SimInspiralChooseFDWaveformSequence`),
  `FrequencyGrid` and `polarizations_on_grid`, `modes` and
  `choose_fd_modes` (two-sided, as LAL's `SphHarmFrequencySeries`),
  `td_polarizations` / `td_modes` (`SimInspiralChooseTDWaveform` /
  `ChooseTDModes`, NRSur7dq4 with a workspace per thread), `simulate_td`
  and `simulate_fd`, and batches (`batch_polarizations` into
  `(n_binaries, n_frequencies)` numpy arrays, `batch_td_polarizations`)
  with `threads=`, per-binary errors or `strict=True`. The GIL is released
  while computing; `Error` maps to `waveformrs.WaveformError` subclasses
  that also derive from `ValueError`, `RuntimeError` or `OSError`. Tested
  against lalsimulation's Python API (`python/tests`); see
  `python/README.md`. `pnr_parity_fix` is one of the options.
- **NRTidal binary-neutron-star models**: the approximants
  `IMRPhenomD_NRTidalv2`, `IMRPhenomXAS_NRTidalv2`, `IMRPhenomXAS_NRTidalv3`,
  `IMRPhenomXP_NRTidalv2` and `IMRPhenomXP_NRTidalv3`, every model of
  LALSimulation 6.2.0 that adds NRTidalv2 or NRTidalv3 to a model of the
  crate (IMRPhenomXHM, IMRPhenomXPHM, IMRPhenomXPNR, TaylorF2 and NRSur7dq4
  have none in LAL). The tidal deformabilities are
  `WaveformOptions::lambda_1` and `lambda_2` (LAL's `lambda1`, `lambda2`,
  also `lambda_1=` / `lambda_2=` in Python); as in LAL, the other
  approximants reject non-zero values, and the quadrupole-monopole
  parameters follow from the deformabilities
  (`XLALSimInspiralSetQuadMonParamsFromLambdas`).
  - `waveformrs::nrtidal`: `LALSimNRTunedTides.c` (the NRTidalv2 and
    NRTidalv3 phase, the amplitude correction, the merger frequencies and
    the Planck window, with NRTidalv3 held at its first minimum after
    merger on the grid as in LAL), `XLALSimInspiralEOSQfromLambda` and the
    3.5PN spin-squared and spin-cubed terms.
  - `IMRPhenomD::with_nrtidal_v2`: the quadrupole-monopole terms in the
    post-Newtonian phase of IMRPhenomD (`Subtract3PNSS` keeps the
    black-hole values, as in LAL) and the NRTidalv2 corrections of
    `XLALSimIMRPhenomDNRTidalFrequencySequence`.
  - `IMRPhenomXAS::with_nrtidal` and `IMRPhenomXP::with_nrtidal`: the tidal
    corrections of `IMRPhenomXASGenerateFD` and `IMRPhenomXPGenerateFD`
    (the frequency sequences of LAL), including the time shift at the
    merger frequency (or at the last frequency of a grid that ends below
    it) and the tidal reference phase; the SpinTaylor angles integrate the
    tidal and quadrupole-monopole terms. On uniform grids,
    `FdWaveform::polarizations_uniform` reproduces
    `SimInspiralChooseFDWaveform`, which builds IMRPhenomXAS_NRTidal with
    IMRPhenomXHM and IMRPhenomXP_NRTidal with IMRPhenomXPHM restricted to
    the (2,2) mode (multibanding off), and IMRPhenomD_NRTidalv2 up to 1.3
    times the merger frequency. As in LAL, `f_ref = 0` is `f_min` in
    `SimInspiralTD` for IMRPhenomXP_NRTidalv2 but not for
    IMRPhenomXP_NRTidalv3.
  - Against LALSimulation 6.2.0 (`tests/model.rs`, suites
    `imrphenomd_nrtidalv2`, `imrphenomxas_nrtidalv2/v3`,
    `imrphenomxp_nrtidalv2/v3` and `nrtidal_uniform` of
    `scripts/generate_reference_data.py`, binary neutron stars through
    merger and a black hole with a neutron star): maximum relative errors
    of 1.8e-11 (IMRPhenomD_NRTidalv2, whose phase reaches ~1e5 rad),
    1.2e-12 (IMRPhenomXAS) and 1.2e-12 (IMRPhenomXP, all three
    precession versions), with mismatches below 1e-15. The start of the
    SpinTaylor angles is sensitive to the last bit of the
    quadrupole-monopole parameter (a change of one ulp in the
    deformability moves the phase at `f_min` by ~1e-8), so it is computed
    with the fused multiply-adds of LAL's arm64 build.
  - `docs/lal-nrtidalv3-integer-division.md`: LAL's NRTidalv3 phase
    derivative contains `(4/3)`, an integer division (= 1); it has no
    effect on LAL's waveforms and is kept here to match.

- **Continuous integration** (`.github/workflows/ci.yml`): build warnings,
  clippy and rustfmt; `cargo test --release` on macOS arm64 (strict
  tolerances) and Linux x86-64 and arm64 (the cross-build tolerances of
  `tests/common`), with and without the `fft` feature; and the Python tests
  against lalsimulation 6.2.0 from conda-forge on macOS and Linux. The
  NRSur7dq4 data are downloaded from Zenodo and cached. The SpinTaylor
  suites of IMRPhenomXP_NRTidalv2/v3 have cross-build tolerances of 1e-6
  (6e-8 seen on Linux arm64).

### Fixed

- **IMRPhenomXHM returns an error instead of panicking** through
  `FrequencyDomainModel` when `f_ref = 0` and the first frequency is
  negative (it becomes the reference frequency, which LAL refuses). The
  inherent `IMRPhenomXHM::polarizations` and `polarizations_on_grid` still
  panic, as documented. No model panics on frequencies LAL would refuse
  (`tests/model.rs`).
- **The crate package leaves out** `tests/data`, `report/`, `scripts/` and
  `lal_data/`, which made the source distribution of the Python bindings
  14.6 MB instead of 0.5 MB.

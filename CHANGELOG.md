# Changelog

All notable changes to `waveformrs` are recorded here.

## Unreleased

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

# waveformrs

Gravitational-wave waveform models for compact binary coalescences, written
in Rust, with Python bindings. The models reproduce
[LALSimulation](https://git.ligo.org/lscsoft/lalsuite) and use its names,
parameters and conventions. Every model is checked against LALSimulation in
the test suite.

## Models

Frequency domain:

- **TaylorF2**: post-Newtonian inspiral.
- **IMRPhenomD**: aligned-spin inspiral, merger and ringdown.
- **IMRPhenomXAS**: aligned spins, dominant (2,2) mode.
- **IMRPhenomXHM**: aligned spins with higher modes.
- **IMRPhenomXP**: precessing spins, dominant mode. Choice of NNLO, MSA or
  SpinTaylor precession angles.
- **IMRPhenomXPHM**: precessing spins with higher modes, with the same
  choice of precession angles.
- **IMRPhenomXPNR**: precessing spins tuned to numerical relativity. It
  includes the parity fix from the development version of LALSuite, so it
  differs from the IMRPhenomXPNR of the LALSimulation releases.

Binary neutron stars, with the NRTidal tidal corrections:

- **IMRPhenomD_NRTidalv2**: IMRPhenomD with NRTidalv2.
- **IMRPhenomXAS_NRTidalv2** and **IMRPhenomXAS_NRTidalv3**: aligned spins,
  dominant mode, with NRTidalv2 or NRTidalv3.
- **IMRPhenomXP_NRTidalv2** and **IMRPhenomXP_NRTidalv3**: precessing spins,
  dominant mode, with NRTidalv2 or NRTidalv3 and the same choice of
  precession angles as IMRPhenomXP.

Time domain:

- **NRSur7dq4**: numerical-relativity surrogate for precessing binaries. It
  needs the surrogate data file from LALSuite's data repository.

## Features

- **One interface for every model.** Choose a model by its LALSimulation
  name, describe the binary once, and get the plus and cross polarizations.
- **Individual modes.** Get the spherical-harmonic modes as well as the
  polarizations.
- **Fixed frequency grids.** Evaluate many waveforms quickly on the same
  frequencies, as samplers do.
- **Batches.** Evaluate many binaries at once across several threads.
- **Time and frequency domains.** Convert between them, following
  LALSimulation's conditioning.
- **Options.** Settings matching LALSimulation's, including the choice of
  precession angles, a restricted set of modes, whether a model may
  switch precession versions automatically, and the tidal deformabilities
  of neutron stars.
- **Python bindings.** The same models and features from Python with numpy
  arrays, including batches across threads. See
  [python/README.md](python/README.md).
- **Platforms.** Tested on macOS and Linux for every pull request on GitHub.

## Getting started

```sh
cargo test --release           # the Rust crate and its tests
cd python && maturin develop --release   # the Python package
```

## More information

- [CHANGELOG.md](CHANGELOG.md): what has changed.
- [python/README.md](python/README.md): the Python package.
- [docs/linux.md](docs/linux.md): running on Linux.
- [docs/lal-quirks.md](docs/lal-quirks.md): where waveformrs differs from
  LALSimulation, and the oddities of LALSimulation it keeps.
- [report/waveformrs-report.html](report/waveformrs-report.html): the
  technical report, with validation and performance details.

## License

waveformrs is free software, licensed under the GNU General Public License
version 3 or later ([LICENSE](LICENSE)). Much of it is translated from
LALSuite and the GNU Scientific Library, whose licenses require this;
[NOTICE.md](NOTICE.md) lists the original code and its authors.

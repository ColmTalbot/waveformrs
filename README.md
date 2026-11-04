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
  SpinTaylor precession angles, plus a variant of the SpinTaylor angles
  that changes smoothly with the masses and spins (LAL's version jumps in
  small steps, which upsets samplers).
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
  needs the surrogate data file from LALSuite's data repository. It can be
  sampled on any time grid, such as each detector's own samples, and made
  to change smoothly with the masses and spins (LAL's version jumps in
  small steps).

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
  arrays, including batches across threads. Free-threaded Python is
  supported. See
  [python/README.md](python/README.md).
- **Platforms.** Tested on macOS and Linux for every pull request on GitHub.

## Installation

The crate is on crates.io and the Python package on PyPI:

```sh
cargo add waveformrs           # the Rust crate
pip install waveformrs         # the Python package
```

## Usage

IMRPhenomXPNR in the frequency domain and NRSur7dq4 in the time domain, for
the same binary. NRSur7dq4 needs its data file, `NRSur7dq4_v1.0.h5` from
[zenodo.org/records/14999310](https://zenodo.org/records/14999310): set
`NRSUR7DQ4_DATA` to the file or `LAL_DATA_PATH` to its directory.

In Rust:

```rust
use waveformrs::{
    choose_fd_model, choose_td_model, Approximant, Binary, FrequencyDomainModel, WaveformOptions,
};

let binary = Binary {
    mass_1: 36.0,                // solar masses
    mass_2: 29.0,
    spin_1: [0.4, 0.2, 0.3],     // dimensionless, at f_ref
    spin_2: [-0.2, 0.3, -0.1],
    luminosity_distance: 500.0,  // Mpc
    inclination: 0.7,            // rad
    phi_ref: 1.2,                // rad
    f_ref: 20.0,                 // Hz
};
let options = WaveformOptions::default();

// IMRPhenomXPNR at the frequencies of your choice (Hz)
let frequencies: Vec<f64> = (160..=8192).map(|k| k as f64 / 8.0).collect();
let model = choose_fd_model(Approximant::IMRPhenomXPNR, &binary, &options)?;
let (hplus, hcross) = model.polarizations(&frequencies)?;

// NRSur7dq4 sampled at 4096 Hz from 20 Hz
let model = choose_td_model(Approximant::NRSur7dq4, &binary, &options)?;
let h = model.polarizations(1.0 / 4096.0, 20.0)?;
// h.hplus, h.hcross; h.epoch is the time of the first sample from the peak
```

In Python:

```python
import numpy as np
import waveformrs as wf

binary = dict(
    mass_1=36.0, mass_2=29.0,
    spin_1=(0.4, 0.2, 0.3), spin_2=(-0.2, 0.3, -0.1),
    luminosity_distance=500.0, inclination=0.7, phi_ref=1.2, f_ref=20.0,
)

frequencies = np.arange(160, 8193) / 8.0
hplus, hcross = wf.polarizations("IMRPhenomXPNR", frequencies, **binary)

h = wf.td_polarizations("NRSur7dq4", 1 / 4096, 20.0, **binary)
hplus, hcross = h
h.times                         # seconds from the peak
```

Each model takes options matching LALSimulation's, such as
`nrsur_continuous` for the smoothly changing NRSur7dq4.

## Getting started from a checkout

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

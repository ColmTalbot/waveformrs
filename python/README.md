# waveformrs for Python

Python bindings of the Rust crate `waveformrs` (the parent directory; see its
`CHANGELOG.md` and `report/`), gravitational-wave models that
reproduce LALSimulation 6.2.0 to rounding: TaylorF2, IMRPhenomD, IMRPhenomXAS,
IMRPhenomXHM, IMRPhenomXP, IMRPhenomXPHM and IMRPhenomXPNR in the frequency
domain, the binary-neutron-star models IMRPhenomD_NRTidalv2,
IMRPhenomXAS_NRTidalv2/v3 and IMRPhenomXP_NRTidalv2/v3, and NRSur7dq4 in the
time domain, named as in LALSimulation.

The bindings use [PyO3](https://pyo3.rs) and
[rust-numpy](https://github.com/PyO3/rust-numpy). Every function releases the
GIL while it computes and returns numpy arrays; the batch functions spread
many binaries over threads and write straight into numpy-owned memory.

## Installation

From this directory (`python/`), in the environment where the package should
go (a virtualenv or an activated conda environment; for an environment that
is not activated, set `CONDA_PREFIX` to it):

```sh
pip install maturin           # >= 1.9.4
maturin develop --release     # build and install in editable mode
```

`maturin build --release` builds a wheel into `target/wheels/` instead, and
`pip install .` builds and installs one in a step. The extension is built
against Python's stable ABI (`abi3`), so one wheel serves CPython 3.9 and
later on a platform.

For this repository's development environment:

```sh
cd python
CONDA_PREFIX=/Users/colmtalbot/mambaforge/envs/waveformrs \
    /Users/colmtalbot/mambaforge/envs/waveformrs/bin/maturin develop --release
```

## Usage

```python
import numpy as np
import waveformrs as wf

wf.approximants()
# ['TaylorF2', 'IMRPhenomD', 'IMRPhenomXAS', 'IMRPhenomXHM', 'IMRPhenomXP',
#  'IMRPhenomXPHM', 'IMRPhenomXPNR', 'NRSur7dq4', 'IMRPhenomD_NRTidalv2',
#  'IMRPhenomXAS_NRTidalv2', 'IMRPhenomXAS_NRTidalv3', 'IMRPhenomXP_NRTidalv2',
#  'IMRPhenomXP_NRTidalv3']

binary = dict(
    mass_1=36.0,                 # solar masses
    mass_2=29.0,
    spin_1=(0.4, 0.2, 0.3),      # dimensionless (x, y, z) at f_ref, L along z
    spin_2=(-0.2, 0.3, -0.1),
    luminosity_distance=500.0,   # Mpc
    inclination=0.7,             # rad
    phi_ref=1.2,                 # orbital reference phase, rad
    f_ref=20.0,                  # Hz (0: LAL's default of the model)
)

# SimInspiralChooseFDWaveformSequence
frequencies = np.linspace(20.0, 1024.0, 8033)
hplus, hcross = wf.polarizations("IMRPhenomXPHM", frequencies, **binary, precession_version=223)
```

The binary is given by keyword arguments with the names of the Rust `Binary`:
`mass_1` and `mass_2` are required; the spins default to zero, the distance
to 1 Mpc, the angles and `f_ref` to zero. `f_ref = 0` selects LAL's default
reference frequency of each model (the first frequency for the IMRPhenom
models, coalescence for TaylorF2).

### Options

The settings of LAL's `LALDict` are `WaveformOptions`, given as keyword
arguments of any waveform function or as an object (`options=`); keywords
override the object.

| option | meaning | default |
|---|---|---|
| `precession_version` | Euler angles of IMRPhenomXP/XPHM, LAL's `PhenomXPrecVersion`: `102` (`"NNLO"`), `223` (`"MSA"`), `320` (`"SpinTaylor"`) | `320` |
| `mode_array` | LAL's `ModeArray`: `[(l, m), ...]`, or `None` for all modes | `None` |
| `version_fallback` | whether IMRPhenomXP/XPHM/XPNR may change the precession version as LAL does; `False` raises `PrecessionVersionUnavailableError` instead | `True` |
| `pnr_parity_fix` | whether IMRPhenomXPNR applies the fix of lalsuite merge request !2593 to its antisymmetric (2,2) mode, restoring the mirror symmetry of the waveform (as LAL master); `False` reproduces LALSimulation 6.2.0 and 6.2.1 | `True` |
| `lambda_1`, `lambda_2` | LAL's `lambda1`, `lambda2`: the dimensionless tidal deformabilities of the bodies of `mass_1` and `mass_2`, for the NRTidal approximants (the others reject non-zero values) | `0.0` |

```python
options = wf.WaveformOptions(precession_version="MSA", mode_array=[(2, 2), (3, 3)])
hp, hc = wf.polarizations("IMRPhenomXPHM", frequencies, **binary, options=options)
hp, hc = wf.polarizations("IMRPhenomXPHM", frequencies, **binary, options=options, version_fallback=False)

wf.WaveformOptions.lalsimulation_defaults()   # LAL 6.2.0 with an empty LALDict: MSA angles, no XPNR parity fix
options.replace(mode_array=None)              # a changed copy
options.as_dict()                             # {'precession_version': 223, ...}
```

Note that the default precession version is the crate's (SpinTaylor, 320),
not LAL 6.2.0's (MSA, 223), and that LAL's multibanding is not implemented:
the models reproduce LAL with `PhenomXHMThresholdMband = PhenomXPHMThresholdMband = 0`.

### Frequency grids

For many waveforms on the same frequencies, a `FrequencyGrid` precomputes the
powers of the frequency that the models need. The grid paths agree with
`polarizations` to about 1e-13 relative rather than bit for bit.

```python
grid = wf.FrequencyGrid.uniform(20.0, 1024.0, 0.125)    # or wf.FrequencyGrid(frequencies)
aligned = dict(binary, spin_1=(0.0, 0.0, 0.3), spin_2=(0.0, 0.0, -0.1))   # XAS takes aligned spins only
hp, hc = wf.polarizations_on_grid("IMRPhenomXAS", grid, **aligned)
grid.frequencies, len(grid)
```

### Modes

`modes` returns the modes `h_lm` without the spherical harmonics at `+f` and
`-f`, in the conventions of `SimInspiralChooseFDModes`, and `choose_fd_modes`
reproduces `SimInspiralChooseFDModes` itself (IMRPhenomXHM, IMRPhenomXPHM,
IMRPhenomXPNR) on its uniform grid. Both return a `FrequencyModes`, laid out
as LAL's `SphHarmFrequencySeries`: each mode on the two-sided frequencies
`-f[::-1], f` (with `f = 0` once).

```python
two_sided, hlm = wf.modes("IMRPhenomXHM", np.linspace(20.0, 512.0, 1000), mass_1=30.0, mass_2=20.0)
hlm[2, -2]                   # h_{2,-2} on the frequencies two_sided

result = wf.choose_fd_modes("IMRPhenomXPHM", 0.25, 20.0, 1024.0, **binary)
result.frequencies, result.modes, result.theta, result.phi, result.polarization_angle
hp, hc = result.polarizations()        # summed with the spherical harmonics (edit result.modes to leave modes out)
```

### Time domain and conversions

NRSur7dq4 needs its data file `NRSur7dq4_v1.0.h5` (zenodo.org/records/14999310):
set `NRSUR7DQ4_DATA` to the file or `LAL_DATA_PATH` to its directory before
the first call (`wf.nrsur7dq4_data_file()` shows what is found).

```python
# SimInspiralChooseTDWaveform (also for the frequency-domain models but TaylorF2,
# through XLALSimInspiralTDFromFD as in LAL)
h = wf.td_polarizations("NRSur7dq4", 1 / 4096, 20.0, **binary)
hplus, hcross = h
h.epoch, h.delta_t, h.times         # first sample (s, relative to the peak), sampling, sample times

modes = wf.td_modes("NRSur7dq4", 1 / 4096, 20.0, **binary)   # SimInspiralChooseTDModes
times, hlm = modes

# SimInspiralTD and SimInspiralFD: conditioned waveforms of any approximant
td = wf.simulate_td("IMRPhenomXPHM", 1 / 4096, 20.0, **binary)
fd = wf.simulate_fd("NRSur7dq4", 0.0, 20.0, 2048.0, **binary)    # delta_f = 0: automatic
fd.hplus, fd.frequencies, fd.delta_f, fd.epoch
```

`td_polarizations` keeps an NRSur7dq4 workspace per thread: repeated calls
allocate nothing, and calls that change only the inclination, phase, distance
or sampling skip the surrogate evaluation.

### Batches

`batch_polarizations` evaluates many binaries on one grid and returns
`(n_binaries, n_frequencies)` arrays. The binaries are a dict of arrays (or
scalars, broadcast), a numpy structured array, a pandas `DataFrame` or an
`(n, 12)` table (`wf.parameter_table`); the spins are `(n, 3)` arrays or their
components `spin_1x`, `spin_1y`, ... The results do not depend on `threads`
(0 for all cores) and each row equals the single `polarizations_on_grid`
call (`exact=True`: the single `polarizations` call) bit for bit.

```python
n = 1000
rng = np.random.default_rng(1)
parameters = dict(
    mass_1=rng.uniform(20, 50, n),
    mass_2=rng.uniform(10, 20, n),
    spin_1=rng.uniform(-0.5, 0.5, (n, 3)),
    spin_2=rng.uniform(-0.5, 0.5, (n, 3)),
    luminosity_distance=400.0,
    inclination=rng.uniform(0, np.pi, n),
    phi_ref=0.0,
    f_ref=20.0,
)
result = wf.batch_polarizations("IMRPhenomXPHM", grid, parameters, threads=8)
hplus, hcross = result                    # complex128, (1000, 8033)
result.ok                                 # False for the binaries that failed (their rows are zero)
result.errors                             # {index: exception}

wf.batch_polarizations("IMRPhenomXPHM", grid, parameters, strict=True)   # raise the first error

# NRSur7dq4 (or any approximant) in the time domain: one workspace per thread
heavy = dict(parameters, mass_1=rng.uniform(50, 70, n), mass_2=rng.uniform(25, 35, n))
td = wf.batch_td_polarizations("NRSur7dq4", 1 / 4096, 20.0, heavy, threads=8)
td[0].hplus, td[0].epoch, td.ok, td.errors      # td[i] is None where td.ok[i] is False
```

Python threads also run in parallel, since the GIL is released:

```python
from concurrent.futures import ThreadPoolExecutor

with ThreadPoolExecutor(4) as pool:
    waveforms = list(pool.map(lambda m: wf.polarizations("IMRPhenomXHM", frequencies, mass_1=m, mass_2=20.0),
                              np.linspace(20, 60, 100)))
```

### Errors

The crate's errors are raised as subclasses of `waveformrs.WaveformError`
(in `waveformrs.exceptions`) that also derive from the matching built-in
exception:

| exception | built-in base | raised for |
|---|---|---|
| `InvalidParameterError` | `ValueError` | parameters LAL rejects (masses, spins, frequencies, in-plane spins for an aligned-spin model, mode arrays, sampling) |
| `UnsupportedError` | `ValueError` | requests the approximant does not support (a mode array for a (2,2) model, `ChooseFDModes` for IMRPhenomD, ...) |
| `UnknownApproximantError` | `ValueError` | names that are not approximants (case-sensitive, as LAL) |
| `ModelError` | `RuntimeError` | failures of the model (SpinTaylor or PNR angles, NRSur7dq4) |
| `PrecessionVersionUnavailableError` | `RuntimeError` (via `ModelError`) | `version_fallback=False` and the requested angles cannot be used |
| `DataFileError`, `DataFileNotFoundError` | `OSError`, `FileNotFoundError` | the NRSur7dq4 data |

Missing or unknown keyword arguments raise `TypeError`.

### Reference

| function | LALSimulation | returns |
|---|---|---|
| `polarizations(approximant, frequencies, **binary)` | `SimInspiralChooseFDWaveformSequence` | `(hplus, hcross)` |
| `polarizations_on_grid(approximant, grid, **binary)` | (fast grid paths) | `(hplus, hcross)` |
| `modes(approximant, frequencies, **binary)` | conventions of `SimInspiralChooseFDModes` | `FrequencyModes` |
| `choose_fd_modes(approximant, delta_f, f_min, f_max, **binary)` | `SimInspiralChooseFDModes` | `FrequencyModes` |
| `td_polarizations(approximant, delta_t, f_min, **binary)` | `SimInspiralChooseTDWaveform` | `TimeSeriesPolarizations` |
| `td_modes(approximant, delta_t, f_min, **binary)` | `SimInspiralChooseTDModes` (NRSur7dq4) | `TimeSeriesModes` |
| `simulate_td(approximant, delta_t, f_min, **binary)` | `SimInspiralTD` | `TimeSeriesPolarizations` |
| `simulate_fd(approximant, delta_f, f_min, f_max, **binary)` | `SimInspiralFD` | `FrequencySeriesPolarizations` |
| `batch_polarizations(approximant, grid, parameters, threads=1, strict=False, exact=False)` | | `BatchPolarizations` |
| `batch_td_polarizations(approximant, delta_t, f_min, parameters, threads=1, strict=False)` | | `BatchTimeSeries` |

The result classes unpack like tuples of their main arrays (`hplus, hcross =
...`, `frequencies, modes = ...`, `times, modes = ...`) and keep the sampling
metadata as attributes. Docstrings and the type stub `_waveformrs.pyi` give
the details.

## Tests

The tests compare with lalsimulation's Python API (python-lalsimulation
6.2.0), at the crate's tolerances: a relative error (max |difference| / max
|LAL|) of 1e-10 for TaylorF2 and IMRPhenomD and 1e-9 for the IMRPhenomX
models, and a mismatch of 1e-12, on macOS arm64; elsewhere the crate's
cross-build tolerances. They also check that batches do not depend on the
number of threads and equal the single calls, the exceptions, and that the
GIL is released.

```sh
cd python
OMP_NUM_THREADS=1 NRSUR7DQ4_DATA=/path/to/NRSur7dq4_v1.0.h5 python -m pytest tests
```

The NRSur7dq4 tests are skipped without the data (found through
`NRSUR7DQ4_DATA`, `LAL_DATA_PATH` or `../lal_data/`); the comparisons with LAL
are skipped without lalsimulation. LAL runs with `OMP_NUM_THREADS=1` (the tests
set it), since the OpenMP loop of IMRPhenomXP races in the MSA and PNR angles.

## Benchmark

`benchmarks/benchmark.py` compares the Python-level call overhead and the
throughput with lalsimulation's Python API, single-threaded and batched:

```sh
OMP_NUM_THREADS=1 NRSUR7DQ4_DATA=... LAL_DATA_PATH=... python benchmarks/benchmark.py [--quick] [--json out.json]
```

Indicative results (Apple M2 Ultra, 24 cores, with other jobs running;
lalsimulation 6.2.0 single-threaded with multibanding off; median of 7 runs):

| call overhead, 4 frequencies | waveformrs | lalsimulation |
|---|---|---|
| TaylorF2 | 1.9 µs | 2.7 µs |
| IMRPhenomD | 2.3 µs | 45 µs |
| IMRPhenomXAS | 2.9 µs | 6.6 µs |

| one waveform, 8033 frequencies, one thread | `polarizations` | `polarizations_on_grid` | lalsimulation | speed-up |
|---|---|---|---|---|
| TaylorF2 | 154 µs | 52 µs | 335 µs | 2.2× |
| IMRPhenomD | 193 µs | 117 µs | 500 µs | 2.6× |
| IMRPhenomXAS | 257 µs | 141 µs | 471 µs | 1.8× |
| IMRPhenomXHM | 1.46 ms | 1.03 ms | 2.50 ms | 1.7× |
| IMRPhenomXP (320) | 1.12 ms | 0.72 ms | 3.12 ms | 2.8× |
| IMRPhenomXPHM (102) | 3.68 ms | 1.25 ms | 7.90 ms | 2.1× |
| IMRPhenomXPHM (223) | 9.13 ms | 3.90 ms | 15.7 ms | 1.7× |
| IMRPhenomXPHM (320) | 4.83 ms | 2.46 ms | 11.1 ms | 2.3× |
| NRSur7dq4 (time domain, 4096 Hz from 20 Hz) | 1.00 ms | – | 5.67 ms | 5.7× |

| batch of 64, waveforms per second | lalsimulation loop | 1 thread | 4 threads | 8 threads | 24 threads |
|---|---|---|---|---|---|
| IMRPhenomXAS | 2,034 | 6,491 | 18,235 | 25,111 | 32,823 |
| IMRPhenomXHM | 401 | 939 | 3,220 | 6,174 | 10,612 |
| IMRPhenomXPHM (223) | 64 | 265 | 938 | 1,837 | 3,582 |
| IMRPhenomXPHM (320) | 88 | 406 | 1,444 | 2,776 | 5,094 |
| NRSur7dq4 (time domain) | 170 | 994 | 3,645 | 6,790 | 13,100 |

The speed-ups are those of the Rust models (see the crate's report); the
Python layer adds about 2 µs per call. A batch on one thread is faster than
single calls because it uses the grid paths.

## Layout

* `Cargo.toml`: the extension crate `waveformrs-python` (library
  `_waveformrs`), depending on the crate with `waveformrs = { path = ".." }`
  (the root is not a workspace), PyO3 0.29 with `abi3-py39` and numpy 0.29,
  built with `lto = "fat"` as the crate recommends.
* `src/`: the bindings (`lib.rs` functions and `FrequencyGrid`, `params.rs`
  keyword arguments and `WaveformOptions`, `batch.rs`, `results.rs`,
  `errors.rs`, `models.rs`).
* `python/waveformrs/`: the Python package (`__init__.py`, the exceptions,
  the batch helpers, the stub `_waveformrs.pyi` and `py.typed`); maturin
  puts the extension module `waveformrs/_waveformrs*.so` next to them.
* `tests/`, `benchmarks/`.

A new option of the crate's `WaveformOptions` is bound by an arm of
`apply_option` and an entry of `option_items` in `src/params.rs` (plus a
getter and a line of the stub).

#!/usr/bin/env python
"""Python-level call overhead and throughput of waveformrs against
lalsimulation's Python API, single-threaded and batched.

Usage::

    OMP_NUM_THREADS=1 python benchmarks/benchmark.py [--quick] [--threads 1 4 8 0]
        [--json results.json]

* overhead: one call on a handful of frequencies, where the time is the cost
  of the Python call, the argument conversion and the model set-up;
* throughput (one thread): one waveform on 8033 frequencies (20-1024 Hz in
  steps of 1/8 Hz) with ``polarizations`` (``SimInspiralChooseFDWaveformSequence``)
  and ``polarizations_on_grid``, against ``SimInspiralChooseFDWaveformSequence``;
  NRSur7dq4 in the time domain (``SimInspiralChooseTDWaveform``, 4096 Hz from
  20 Hz);
* batches: ``batch_polarizations`` / ``batch_td_polarizations`` of 64 binaries
  with several thread counts, against a Python loop over lalsimulation.

LALSimulation runs with ``OMP_NUM_THREADS=1`` (its OpenMP loops race in the
MSA and PNR angles) and with multibanding off for IMRPhenomXHM/XPHM, which is
what waveformrs reproduces. The timings are the median of several repeats and
are only indicative: other processes on the machine affect them.
"""

import argparse
import json
import os
import platform
import statistics
import time

os.environ.setdefault("OMP_NUM_THREADS", "1")

import lal  # noqa: E402
import lalsimulation as lalsim  # noqa: E402
import numpy as np  # noqa: E402

import waveformrs as wf  # noqa: E402

BINARY = dict(
    mass_1=36.0,
    mass_2=29.0,
    spin_1=(0.4, 0.2, 0.3),
    spin_2=(-0.2, 0.3, -0.1),
    luminosity_distance=500.0,
    inclination=0.7,
    phi_ref=1.2,
    f_ref=20.0,
)
ALIGNED = dict(BINARY, spin_1=(0.0, 0.0, 0.3), spin_2=(0.0, 0.0, -0.1))
NRSUR = dict(BINARY, mass_1=40.0, mass_2=20.0)

# label, approximant, binary, options
FD_MODELS = [
    ("TaylorF2", "TaylorF2", ALIGNED, {}),
    ("IMRPhenomD", "IMRPhenomD", ALIGNED, {}),
    ("IMRPhenomXAS", "IMRPhenomXAS", ALIGNED, {}),
    ("IMRPhenomXHM", "IMRPhenomXHM", ALIGNED, {}),
    ("IMRPhenomXP (320)", "IMRPhenomXP", BINARY, dict(precession_version=320)),
    ("IMRPhenomXPHM (102)", "IMRPhenomXPHM", BINARY, dict(precession_version=102)),
    ("IMRPhenomXPHM (223)", "IMRPhenomXPHM", BINARY, dict(precession_version=223)),
    ("IMRPhenomXPHM (320)", "IMRPhenomXPHM", BINARY, dict(precession_version=320)),
    # with the parity fix of lalsuite !2593, which lalsimulation 6.2.0 lacks
    # (the cost is the same)
    ("IMRPhenomXPNR", "IMRPhenomXPNR", BINARY, {}),
]


def median_time(function, repeats, number=1):
    """Median over `repeats` of the mean time of `number` calls (s)."""
    function()
    times = []
    for _ in range(repeats):
        start = time.perf_counter()
        for _ in range(number):
            function()
        times.append((time.perf_counter() - start) / number)
    return statistics.median(times)


def lal_params(options):
    params = lal.CreateDict()
    lalsim.SimInspiralWaveformParamsInsertPhenomXHMThresholdMband(params, 0)
    lalsim.SimInspiralWaveformParamsInsertPhenomXPHMThresholdMband(params, 0)
    if "precession_version" in options:
        lalsim.SimInspiralWaveformParamsInsertPhenomXPrecVersion(
            params, options["precession_version"]
        )
    return params


def lal_sequence(approximant, binary, options, frequencies):
    """A closure calling SimInspiralChooseFDWaveformSequence as a Python user
    would (converting the parameters and the frequencies on each call)."""
    approximant = lalsim.GetApproximantFromString(approximant)

    def call():
        vector = lal.CreateREAL8Vector(len(frequencies))
        vector.data = frequencies
        return lalsim.SimInspiralChooseFDWaveformSequence(
            binary["phi_ref"],
            binary["mass_1"] * lal.MSUN_SI,
            binary["mass_2"] * lal.MSUN_SI,
            *binary["spin_1"],
            *binary["spin_2"],
            binary["f_ref"],
            binary["luminosity_distance"] * 1e6 * lal.PC_SI,
            binary["inclination"],
            lal_params(options),
            approximant,
            vector,
        )

    return call


def lal_nrsur(binary, delta_t, f_min):
    def call():
        return lalsim.SimInspiralChooseTDWaveform(
            binary["mass_1"] * lal.MSUN_SI,
            binary["mass_2"] * lal.MSUN_SI,
            *binary["spin_1"],
            *binary["spin_2"],
            binary["luminosity_distance"] * 1e6 * lal.PC_SI,
            binary["inclination"],
            binary["phi_ref"],
            0.0,
            0.0,
            0.0,
            delta_t,
            f_min,
            binary["f_ref"],
            None,
            lalsim.NRSur7dq4,
        )

    return call


def draws(binary, n, seed=7):
    """`n` binaries around `binary` (masses and angles varied) as a dict."""
    rng = np.random.default_rng(seed)
    precessing = binary["spin_1"][0] != 0.0
    spins = rng.uniform(-0.5, 0.5, (2, n, 3))
    if not precessing:
        spins[..., :2] = 0.0
    return dict(
        mass_1=binary["mass_1"] * rng.uniform(0.9, 1.1, n),
        mass_2=binary["mass_2"] * rng.uniform(0.9, 1.1, n),
        spin_1=spins[0],
        spin_2=spins[1],
        luminosity_distance=binary["luminosity_distance"],
        inclination=rng.uniform(0.0, np.pi, n),
        phi_ref=rng.uniform(0.0, 2 * np.pi, n),
        f_ref=binary["f_ref"],
    )


def row(parameters, i):
    return {k: (v[i] if np.ndim(v) > 0 else v) for k, v in parameters.items()}


def fmt(seconds):
    if seconds < 1e-3:
        return f"{seconds * 1e6:.1f} µs"
    return f"{seconds * 1e3:.2f} ms"


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--quick", action="store_true", help="fewer repeats and binaries")
    parser.add_argument("--threads", type=int, nargs="+", default=[1, 4, 8, 0])
    parser.add_argument("--json", help="write the results to this file")
    args = parser.parse_args()
    repeats = 3 if args.quick else 7
    n_batch = 16 if args.quick else 64
    results = dict(
        machine=f"{platform.machine()} {platform.system()} {platform.release()}",
        cpus=os.cpu_count(),
        lalsimulation=lalsim.__version__,
        waveformrs=wf.__version__,
        overhead={},
        throughput={},
        batch={},
    )
    print(f"waveformrs {wf.__version__}, lalsimulation {lalsim.__version__}, "
          f"{results['machine']}, {results['cpus']} CPUs; median of {repeats} repeats\n")

    # call overhead: a few frequencies
    few = np.array([20.0, 40.0, 80.0, 160.0])
    print("Call overhead (4 frequencies)")
    print(f"{'model':<22}{'waveformrs':>14}{'lalsimulation':>16}")
    for label, approximant, binary, options in FD_MODELS[:3]:
        ours = median_time(
            lambda: wf.polarizations(approximant, few, **binary, **options), repeats, 2000
        )
        theirs = median_time(lal_sequence(approximant, binary, options, few), repeats, 2000)
        results["overhead"][label] = dict(waveformrs=ours, lalsimulation=theirs)
        print(f"{label:<22}{fmt(ours):>14}{fmt(theirs):>16}")

    # one waveform, one thread
    grid = wf.FrequencyGrid.uniform(20.0, 1024.0, 0.125)
    frequencies = grid.frequencies
    print(f"\nOne waveform, one thread ({len(grid)} frequencies)")
    print(f"{'model':<22}{'polarizations':>15}{'on grid':>12}{'lalsimulation':>16}{'speed-up':>10}")
    for label, approximant, binary, options in FD_MODELS:
        number = 20 if "HM" not in approximant else 5
        ours = median_time(
            lambda: wf.polarizations(approximant, frequencies, **binary, **options),
            repeats,
            number,
        )
        on_grid = median_time(
            lambda: wf.polarizations_on_grid(approximant, grid, **binary, **options),
            repeats,
            number,
        )
        theirs = median_time(
            lal_sequence(approximant, binary, options, frequencies), repeats, number
        )
        results["throughput"][label] = dict(
            waveformrs=ours, on_grid=on_grid, lalsimulation=theirs
        )
        print(
            f"{label:<22}{fmt(ours):>15}{fmt(on_grid):>12}{fmt(theirs):>16}"
            f"{theirs / ours:>9.1f}x"
        )
    nrsur = wf.nrsur7dq4_data_file() is not None
    if nrsur:
        delta_t, f_min = 1.0 / 4096, 20.0
        # a different binary each call, so that the cached modes of the
        # workspace are not reused
        masses = iter(np.linspace(39.0, 41.0, 100_000))
        ours = median_time(
            lambda: wf.td_polarizations(
                "NRSur7dq4", delta_t, f_min, **dict(NRSUR, mass_1=next(masses))
            ),
            repeats,
            20,
        )
        theirs = median_time(lal_nrsur(NRSUR, delta_t, f_min), repeats, 20)
        results["throughput"]["NRSur7dq4 (time domain)"] = dict(
            waveformrs=ours, lalsimulation=theirs
        )
        print(
            f"{'NRSur7dq4 (TD)':<22}{fmt(ours):>15}{'-':>12}{fmt(theirs):>16}"
            f"{theirs / ours:>9.1f}x"
        )

    # batches
    print(f"\nBatches of {n_batch} binaries: waveforms per second")
    columns = "".join(f"{f'{t} thr' if t else 'all':>10}" for t in args.threads)
    print(f"{'model':<22}{'lal loop':>10}{columns}")
    batch_models = [m for m in FD_MODELS if m[0] in (
        "IMRPhenomXAS", "IMRPhenomXHM", "IMRPhenomXPHM (223)", "IMRPhenomXPHM (320)",
        "IMRPhenomXPNR")]
    for label, approximant, binary, options in batch_models:
        parameters = draws(binary, n_batch)
        calls = [
            lal_sequence(approximant, row(parameters, i), options, frequencies)
            for i in range(n_batch)
        ]
        lal_rate = n_batch / median_time(lambda: [c() for c in calls], max(2, repeats // 2))
        rates = {}
        for threads in args.threads:
            seconds = median_time(
                lambda: wf.batch_polarizations(
                    approximant, grid, parameters, threads=threads, **options
                ),
                max(2, repeats // 2),
            )
            rates[threads] = n_batch / seconds
        results["batch"][label] = dict(lalsimulation=lal_rate, waveformrs=rates)
        print(f"{label:<22}{lal_rate:>10.0f}" + "".join(f"{rates[t]:>10.0f}" for t in args.threads))
    if nrsur:
        delta_t, f_min = 1.0 / 4096, 20.0
        # heavier than NRSUR, so that every draw starts below 20 Hz
        parameters = draws(dict(NRSUR, mass_1=50.0, mass_2=25.0), n_batch)
        calls = [lal_nrsur(row(parameters, i), delta_t, f_min) for i in range(n_batch)]
        lal_rate = n_batch / median_time(lambda: [c() for c in calls], max(2, repeats // 2))
        rates = {}
        for threads in args.threads:
            seconds = median_time(
                lambda: wf.batch_td_polarizations(
                    "NRSur7dq4", delta_t, f_min, parameters, threads=threads, strict=True
                ),
                max(2, repeats // 2),
            )
            rates[threads] = n_batch / seconds
        results["batch"]["NRSur7dq4 (time domain)"] = dict(lalsimulation=lal_rate, waveformrs=rates)
        print(f"{'NRSur7dq4 (TD)':<22}{lal_rate:>10.0f}" + "".join(f"{rates[t]:>10.0f}" for t in args.threads))

    if args.json:
        with open(args.json, "w") as file:
            json.dump(results, file, indent=2)


if __name__ == "__main__":
    main()

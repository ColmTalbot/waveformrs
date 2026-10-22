#!/usr/bin/env python
"""
Time lalsimulation's conversions between the time and frequency domains for
the configurations of ``examples/conversions.rs``: ``SimInspiralTD`` for
frequency-domain approximants, ``SimInspiralFD`` for NRSur7dq4, and LAL's
real FFTs (FFTW, ``FFTW_ESTIMATE`` plans as the conversions use them).

Usage
-----
    OMP_NUM_THREADS=1 LAL_DATA_PATH=/path/to/lal_data \\
        python scripts/time_conversions_lal.py [--output lal_timings.json]

then ``cargo run --release --example conversions -- lal_timings.json`` prints
both side by side.
"""
import argparse
import json
import time

import lal
import lalsimulation as lalsim
import numpy as np

ALIGNED = ([0.0, 0.0, 0.3], [0.0, 0.0, -0.4])
PRECESSING = ([0.4, 0.2, 0.3], [-0.2, 0.3, -0.1])
NO_MULTIBANDING = dict(PhenomXHMThresholdMband=0, PhenomXPHMThresholdMband=0)

# label, approximant, (m1, m2), spins, delta_t, f_min, LAL flags
TD_CONFIGURATIONS = [
    ("IMRPhenomD 36+29 Msun", "IMRPhenomD", (36.0, 29.0), ALIGNED, 1 / 4096, 20.0, {}),
    ("IMRPhenomXAS 36+29 Msun", "IMRPhenomXAS", (36.0, 29.0), ALIGNED, 1 / 4096, 20.0, {}),
    ("IMRPhenomXHM 36+29 Msun", "IMRPhenomXHM", (36.0, 29.0), ALIGNED, 1 / 4096, 20.0, NO_MULTIBANDING),
    ("IMRPhenomXP 36+29 Msun", "IMRPhenomXP", (36.0, 29.0), PRECESSING, 1 / 4096, 20.0, dict(PhenomXPrecVersion=320)),
    (
        "IMRPhenomXPHM 36+29 Msun",
        "IMRPhenomXPHM",
        (36.0, 29.0),
        PRECESSING,
        1 / 4096,
        20.0,
        dict(PhenomXPrecVersion=320, **NO_MULTIBANDING),
    ),
    ("IMRPhenomD 10+5 Msun, 10 Hz", "IMRPhenomD", (10.0, 5.0), ALIGNED, 1 / 4096, 10.0, {}),
    ("TaylorF2 1.4+1.3 Msun, 30 Hz", "TaylorF2", (1.4, 1.3), ([0, 0, 0], [0, 0, 0]), 1 / 4096, 30.0, {}),
]

# label, (m1, m2), spins, delta_f, f_min, f_max
FD_CONFIGURATIONS = [
    ("NRSur7dq4 40+20 Msun, df = 0, 2048 Hz", (40.0, 20.0), ([0.4, -0.3, 0.2], [0.1, 0.5, -0.3]), 0.0, 20.0, 2048.0),
    ("NRSur7dq4 40+20 Msun, df = 1/8, 2048 Hz", (40.0, 20.0), ([0.4, -0.3, 0.2], [0.1, 0.5, -0.3]), 0.125, 20.0, 2048.0),
]

FFT_LENGTHS = [2**16, 2**20]


def per_call(func, minimum=1.0):
    func()
    n = 0
    start = time.perf_counter()
    while time.perf_counter() - start < minimum:
        func()
        n += 1
    return (time.perf_counter() - start) / n


def lal_dict(flags):
    params = lal.CreateDict()
    for key, value in flags.items():
        getattr(lalsim, f"SimInspiralWaveformParamsInsert{key}")(params, value)
    return params


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", default=None)
    parser.add_argument("--minimum", type=float, default=1.0)
    args = parser.parse_args()
    results = dict()

    for label, approximant, (m1, m2), (s1, s2), delta_t, f_min, flags in TD_CONFIGURATIONS:
        params = lal_dict(flags)
        approx = lalsim.GetApproximantFromString(approximant)

        def call():
            return lalsim.SimInspiralTD(
                m1 * lal.MSUN_SI, m2 * lal.MSUN_SI, *s1, *s2, 500e6 * lal.PC_SI, 0.7, 1.2,
                0.0, 0.0, 0.0, delta_t, f_min, 0.0, params, approx,
            )

        hp, _ = call()
        results[label] = per_call(call, args.minimum)
        print(f"{label:<45} SimInspiralTD {results[label] * 1e3:9.3f} ms ({hp.data.length} samples)")

    for label, (m1, m2), (s1, s2), delta_f, f_min, f_max in FD_CONFIGURATIONS:

        def call():
            return lalsim.SimInspiralFD(
                m1 * lal.MSUN_SI, m2 * lal.MSUN_SI, *s1, *s2, 100e6 * lal.PC_SI, 0.7, 0.1,
                0.0, 0.0, 0.0, delta_f, f_min, f_max, 20.0, None, lalsim.NRSur7dq4,
            )

        hp, _ = call()
        results[label] = per_call(call, args.minimum)
        print(f"{label:<45} SimInspiralFD {results[label] * 1e3:9.3f} ms ({hp.data.length} samples)")

    rng = np.random.default_rng(0)
    for n in FFT_LENGTHS:
        series = lal.CreateREAL8TimeSeries("x", lal.LIGOTimeGPS(0), 0.0, 1 / 4096, lal.StrainUnit, n)
        series.data.data[:] = rng.standard_normal(n)
        spectrum = lal.CreateCOMPLEX16FrequencySeries("y", lal.LIGOTimeGPS(0), 0.0, 1.0, lal.StrainUnit, n // 2 + 1)

        def forward():
            plan = lal.CreateForwardREAL8FFTPlan(n, 0)
            lal.REAL8TimeFreqFFT(spectrum, series, plan)

        def reverse():
            plan = lal.CreateReverseREAL8FFTPlan(n, 0)
            lal.REAL8FreqTimeFFT(series, spectrum, plan)

        forward()
        for name, func in [("forward", forward), ("reverse", reverse)]:
            label = f"{name} real FFT, N = 2^{int(np.log2(n))}"
            results[label] = per_call(func, args.minimum)
            print(f"{label:<45} (plan + FFT) {results[label] * 1e3:9.3f} ms")

    if args.output:
        with open(args.output, "w") as ff:
            json.dump(dict(lalsimulation_version=lalsim.__version__, seconds=results), ff, indent=1)


if __name__ == "__main__":
    main()

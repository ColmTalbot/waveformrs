#!/usr/bin/env python
"""
Time ``SimInspiralSpinTaylorPNEvolveOrbit`` for the same IMRPhenomXPHM-like
calls as ``examples/spintaylor.rs``.

Usage
-----
    python scripts/time_spintaylor_lal.py
"""
import time

import lal
import lalsimulation as lalsim
import numpy as np

CASES = [
    # label, m1, m2, spin_1, spin_2, f_min, f_ref
    ("36+29 Msun, f_min = f_ref = 20 Hz", 36.0, 29.0, [0.4, 0.2, 0.3], [-0.2, 0.3, -0.1], 20.0, 20.0),
    ("60+15 Msun, f_min = 15, f_ref = 20 Hz", 60.0, 15.0, [0.6, -0.3, 0.5], [0.1, 0.1, 0.0], 15.0, 20.0),
    ("12+6 Msun, f_min = f_ref = 20 Hz", 12.0, 6.0, [0.2, 0.7, -0.1], [0.3, -0.3, 0.6], 20.0, 20.0),
]


def evolve(delta_t, m1, m2, f_start, f_end, s1, s2, lnhat=(0, 0, 1), e1=(1, 0, 0)):
    return lalsim.SimInspiralSpinTaylorPNEvolveOrbit(
        delta_t, m1 * lal.MSUN_SI, m2 * lal.MSUN_SI, f_start, f_end,
        *s1, *s2, *lnhat, *e1, 0.0, 0.0, 1.0, 1.0, 6, 12, 7, 0, lalsim.SpinTaylorT4,
    )


def xphm_like(m1, m2, s1, s2, f_min, f_ref):
    total_mass = m1 + m2
    f_cut = 0.2 / (total_mass * lal.MTSUN_SI)
    dt_coarse = 0.5 * 10 / f_cut
    n = 0
    if f_ref > f_min:
        n += evolve(dt_coarse, m1, m2, f_ref, f_min - 0.5, s1, s2)[0].data.length
        coarse = evolve(dt_coarse, m1, m2, f_ref, f_cut, s1, s2)
    else:
        coarse = evolve(dt_coarse, m1, m2, f_min, f_cut, s1, s2)
    low = coarse[0].data.length
    n += low
    index = low - 1 - min(9, low - 1)
    data = [ts.data.data[index] for ts in coarse]
    f_trans = data[0] ** 3 / (np.pi * total_mass * lal.MTSUN_SI)
    fine = evolve(0.5 / f_cut, m1, m2, f_trans, f_cut, data[2:5], data[5:8], data[8:11], data[11:14])
    return n + fine[0].data.length


def time_it(label, func):
    samples = func()
    start = time.perf_counter()
    repeats = 0
    while time.perf_counter() - start < 0.2:
        func()
        repeats += 1
    repeats = max(repeats * 5, 5)
    start = time.perf_counter()
    for _ in range(repeats):
        func()
    per_call = (time.perf_counter() - start) / repeats
    print(f"{label:<45} {per_call * 1e6:9.1f} us/call  ({samples} samples, {repeats} repeats)")


def main():
    print(f"lalsimulation {lalsim.__version__}")
    for label, m1, m2, s1, s2, f_min, f_ref in CASES:
        time_it(label, lambda: xphm_like(m1, m2, s1, s2, f_min, f_ref))
    f_cut = 0.2 / (65.0 * lal.MTSUN_SI)
    s1, s2 = [0.4, 0.2, 0.3], [-0.2, 0.3, -0.1]
    time_it("single evolve_orbit 36+29, 20 Hz -> fCut",
            lambda: evolve(5.0 / f_cut, 36.0, 29.0, 20.0, f_cut, s1, s2)[0].data.length)
    time_it("single evolve_orbit 36+29, 20 Hz -> fCut, dt/10",
            lambda: evolve(0.5 / f_cut, 36.0, 29.0, 20.0, f_cut, s1, s2)[0].data.length)


if __name__ == "__main__":
    main()

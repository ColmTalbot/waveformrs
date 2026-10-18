#!/usr/bin/env python
"""
Time lalsimulation's NRSur7dq4 for the configurations of
``examples/nrsur7dq4.rs``: ``SimInspiralChooseTDWaveform`` and
``SimInspiralChooseTDModes`` (q = 2, generic spins), and a sampler-like
sequence of random binaries.

Usage
-----
    OMP_NUM_THREADS=1 LAL_DATA_PATH=/path/to/lal_data \
        python scripts/time_nrsur7dq4_lal.py
"""
import time

import lal
import lalsimulation as lalsim
import numpy as np

# (total mass in Msun, sampling rate in Hz, f_min in Hz), as in the example
CONFIGURATIONS = [
    (60.0, 4096.0, 20.0),
    (60.0, 16384.0, 20.0),
    (20.0, 4096.0, 55.0),
    (20.0, 16384.0, 55.0),
    (100.0, 4096.0, 20.0),
    (100.0, 1024.0, 20.0),
    (200.0, 2048.0, 20.0),
    (300.0, 1024.0, 20.0),
]
SPIN_1 = [0.4, -0.3, 0.2]
SPIN_2 = [0.1, 0.5, -0.3]


def per_call(func, minimum=1.0):
    func()
    n = 0
    start = time.perf_counter()
    while time.perf_counter() - start < minimum:
        func()
        n += 1
    return (time.perf_counter() - start) / n


def polarizations(m1, m2, s1, s2, delta_t, f_min, f_ref, inclination=0.7, phi_ref=0.1):
    return lalsim.SimInspiralChooseTDWaveform(
        m1 * lal.MSUN_SI, m2 * lal.MSUN_SI, *s1, *s2, 1e8 * lal.PC_SI, inclination,
        phi_ref, 0.0, 0.0, 0.0, delta_t, f_min, f_ref, None, lalsim.NRSur7dq4,
    )


def modes(m1, m2, s1, s2, delta_t, f_min, f_ref):
    return lalsim.SimInspiralChooseTDModes(
        0.0, delta_t, m1 * lal.MSUN_SI, m2 * lal.MSUN_SI, *s1, *s2, f_min, f_ref,
        1e8 * lal.PC_SI, None, 4, lalsim.NRSur7dq4,
    )


def spin(rng):
    magnitude = rng.uniform(0, 0.99)
    cos_theta = rng.uniform(-1, 1)
    phi = rng.uniform(0, 2 * np.pi)
    sin_theta = np.sqrt(1 - cos_theta**2)
    return [magnitude * sin_theta * np.cos(phi), magnitude * sin_theta * np.sin(phi),
            magnitude * cos_theta]


def main():
    start = time.perf_counter()
    polarizations(40.0, 20.0, SPIN_1, SPIN_2, 1 / 4096, 20.0, 0.0)
    print(f"first call (loads the data): {1e3 * (time.perf_counter() - start):.1f} ms")
    print(f"{'M/Msun':>6} {'1/dt':>6} {'f_min':>6} {'samples':>7} | {'pols':>9} {'modes':>9}")
    for total_mass, rate, f_min in CONFIGURATIONS:
        args = (total_mass * 2 / 3, total_mass / 3, SPIN_1, SPIN_2, 1 / rate, f_min, 0.0)
        n = polarizations(*args)[0].data.length
        t_pol = per_call(lambda: polarizations(*args))
        t_modes = per_call(lambda: modes(*args))
        print(f"{total_mass:>6} {rate:>6} {f_min:>6} {n:>7} | {1e3 * t_pol:7.2f}ms {1e3 * t_modes:7.2f}ms")

    rng = np.random.default_rng(1)
    n, rejected, samples = 500, 0, 0
    binaries = list()
    for _ in range(n):
        total_mass = rng.uniform(40, 120)
        q = rng.uniform(1, 6)
        binaries.append((total_mass * q / (1 + q), total_mass / (1 + q), spin(rng), spin(rng),
                         rng.uniform(0, np.pi), rng.uniform(0, 2 * np.pi)))
    start = time.perf_counter()
    for m1, m2, s1, s2, inclination, phi_ref in binaries:
        try:
            hp, _ = polarizations(m1, m2, s1, s2, 1 / 4096, 20.0, 20.0, inclination, phi_ref)
            samples += hp.data.length
        except RuntimeError:
            rejected += 1
    elapsed = time.perf_counter() - start
    print(
        f"sampler-like sequence of {n} binaries (M = 40-120 Msun from 20 Hz, 4096 Hz): "
        f"{1e3 * elapsed / n:.2f} ms per waveform ({samples // (n - rejected)} samples on "
        f"average, {rejected} rejected)"
    )


if __name__ == "__main__":
    main()

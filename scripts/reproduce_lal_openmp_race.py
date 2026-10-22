"""
Reproduce the OpenMP data race in IMRPhenomXP with the MSA Euler angles.

Run the same IMRPhenomXP waveform several times in fresh processes with
different OMP_NUM_THREADS and compare against the single-threaded result.

    python scripts/reproduce_lal_openmp_race.py
"""
import os
import subprocess
import sys
import tempfile

import numpy as np

CHILD = r"""
import sys
import numpy as np
import lal, lalsimulation as ls

version = int(sys.argv[1])
params = lal.CreateDict()
ls.SimInspiralWaveformParamsInsertPhenomXPrecVersion(params, version)
m1, m2 = 30.0, 10.0
spins = (0.4, -0.3, 0.2, 0.1, 0.5, -0.3)
hp, hc = ls.SimInspiralChooseFDWaveform(
    m1 * lal.MSUN_SI, m2 * lal.MSUN_SI, *spins,
    400e6 * lal.PC_SI, 0.7, 0.1, 0.0, 0.0, 0.0,
    0.125, 20.0, 1024.0, 20.0, params, ls.IMRPhenomXP,
)
np.save(sys.argv[2], np.stack([hp.data.data, hc.data.data]))
"""


def run(version, threads, path):
    env = dict(os.environ, OMP_NUM_THREADS=str(threads))
    subprocess.run(
        [sys.executable, "-c", CHILD, str(version), path],
        env=env,
        check=True,
        stderr=subprocess.DEVNULL,
    )
    return np.load(path)


def compare(test, reference):
    diff = np.abs(test - reference).max(axis=0)
    scale = np.abs(reference).max()
    bad = diff > 1e-10 * scale
    return diff.max() / scale, int(bad.sum())


def main():
    tmp = tempfile.mkdtemp()
    path = os.path.join(tmp, "h.npy")
    for version in (223, 102, 320):
        reference = run(version, 1, path)
        repeat = run(version, 1, path)
        rows = [("1 (repeat)",) + compare(repeat, reference)]
        for threads in (2, 4, 8, 8, 8):
            test = run(version, threads, path)
            rows.append((str(threads),) + compare(test, reference))
        print(f"PhenomXPrecVersion {version}: {reference.shape[1]} frequency bins")
        for threads, err, nbad in rows:
            print(f"  OMP_NUM_THREADS={threads:<10} max |dh|/max|h| = {err:.2e}"
                  f"   bins differing by >1e-10: {nbad}")


if __name__ == "__main__":
    main()

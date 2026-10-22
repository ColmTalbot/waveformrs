"""The GIL is released while waveforms are computed."""

import os
import sys
import sysconfig
import threading
import time
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pytest

import waveformrs as wf

BINARY = dict(
    mass_1=36.0,
    mass_2=29.0,
    spin_1=(0.4, 0.2, 0.3),
    spin_2=(-0.2, 0.3, -0.1),
    luminosity_distance=500.0,
    inclination=0.7,
    phi_ref=1.2,
    f_ref=20.0,
    precession_version=223,
)
GRID = wf.FrequencyGrid.uniform(20.0, 1024.0, 0.125)


@pytest.mark.skipif(
    not sysconfig.get_config_var("Py_GIL_DISABLED"), reason="needs free-threaded Python"
)
def test_import_keeps_the_gil_disabled():
    """The extension declares free-threading support, so importing it does
    not switch the GIL back on (unless PYTHON_GIL=1 forces it)."""
    if os.environ.get("PYTHON_GIL") == "1":
        pytest.skip("PYTHON_GIL=1")
    assert not sys._is_gil_enabled()


def python_work(n=200_000):
    """Pure-Python work that needs the GIL."""
    total = 0
    for i in range(n):
        total += i
    return total


def test_python_runs_during_a_call():
    """Python code in another thread runs while a long call computes: with
    the GIL held it could only finish after the call."""
    parameters = dict(mass_1=np.full(48, 36.0), mass_2=29.0, spin_1=(0.4, 0.2, 0.3))
    # make sure the call takes a while (and warm up)
    start = time.perf_counter()
    wf.batch_polarizations("IMRPhenomXPHM", GRID, parameters, precession_version=223)
    duration = time.perf_counter() - start
    start = time.perf_counter()
    python_work()
    work = time.perf_counter() - start
    assert work < duration / 3, "the call is too short for the test"

    finished = {}

    def call():
        wf.batch_polarizations("IMRPhenomXPHM", GRID, parameters, precession_version=223)
        finished["call"] = time.perf_counter()

    thread = threading.Thread(target=call)
    thread.start()
    time.sleep(0.01)  # let the call start
    python_work()
    finished["work"] = time.perf_counter()
    thread.join()
    assert finished["work"] < finished["call"]


@pytest.mark.skipif((os.cpu_count() or 1) < 4, reason="needs four cores")
def test_python_threads_speed_up():
    """Calls from four Python threads run in parallel (indicative: the
    machine may be loaded, so the bound is loose)."""
    frequencies = GRID.frequencies
    n = 16

    def one(_):
        return wf.polarizations("IMRPhenomXPHM", frequencies, **BINARY)

    one(0)
    start = time.perf_counter()
    serial = [one(i) for i in range(n)]
    serial_time = time.perf_counter() - start
    with ThreadPoolExecutor(4) as pool:
        start = time.perf_counter()
        threaded = list(pool.map(one, range(n)))
        threaded_time = time.perf_counter() - start
    for (a, b), (c, d) in zip(serial, threaded):
        np.testing.assert_array_equal(a, c)
        np.testing.assert_array_equal(b, d)
    speedup = serial_time / threaded_time
    print(f"four Python threads: {speedup:.2f}x")
    assert speedup > 1.5, speedup

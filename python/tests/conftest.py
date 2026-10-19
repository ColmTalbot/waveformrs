"""Shared set-up of the tests.

* LALSimulation runs single-threaded (``OMP_NUM_THREADS=1``): its MSA and PNR
  angles race in the OpenMP loop of IMRPhenomXP.
* The NRSur7dq4 data are found through ``NRSUR7DQ4_DATA`` (the file) or
  ``LAL_DATA_PATH`` (its directory), else in ``lal_data/`` of the repository;
  ``LAL_DATA_PATH`` is set to the same directory for lalsimulation.
"""

import os
import platform
import sys
from pathlib import Path

# before lalsimulation (and its OpenMP runtime) is loaded
os.environ["OMP_NUM_THREADS"] = "1"

import numpy as np
import pytest

DATA_FILE = "NRSur7dq4_v1.0.h5"


def _nrsur_data():
    """The NRSur7dq4 data file, or None."""
    path = os.environ.get("NRSUR7DQ4_DATA")
    if path and Path(path).is_file():
        return Path(path)
    for directory in os.environ.get("LAL_DATA_PATH", "").split(os.pathsep):
        if directory and (Path(directory) / DATA_FILE).is_file():
            return Path(directory) / DATA_FILE
    candidate = Path(__file__).resolve().parents[2] / "lal_data" / DATA_FILE
    if candidate.is_file():
        return candidate
    return None


NRSUR_DATA = _nrsur_data()
if NRSUR_DATA is not None:
    os.environ["NRSUR7DQ4_DATA"] = str(NRSUR_DATA)
    paths = [p for p in os.environ.get("LAL_DATA_PATH", "").split(os.pathsep) if p]
    if str(NRSUR_DATA.parent) not in paths:
        os.environ["LAL_DATA_PATH"] = os.pathsep.join([str(NRSUR_DATA.parent), *paths])

#: Whether this is the build of lalsimulation that the crate mirrors
#: operation by operation (conda-forge, macOS arm64): elsewhere the
#: ill-conditioned parts of the higher-mode models amplify last-bit
#: differences (see docs/linux.md of the crate).
MIRRORED_BUILD = sys.platform == "darwin" and platform.machine() == "arm64"

requires_nrsur = pytest.mark.skipif(
    NRSUR_DATA is None,
    reason=f"{DATA_FILE} not found (set NRSUR7DQ4_DATA or LAL_DATA_PATH)",
)


@pytest.fixture(scope="session")
def lal():
    """The lal and lalsimulation modules (skips without them)."""
    lal = pytest.importorskip("lal")
    lalsim = pytest.importorskip("lalsimulation")
    return lal, lalsim


def compare(test, reference):
    """``(max |a - b| / max |b|, mismatch)`` with uniform weights."""
    test = np.asarray(test)
    reference = np.asarray(reference)
    assert test.shape == reference.shape
    scale = np.max(np.abs(reference))
    error = np.max(np.abs(test - reference)) / scale if scale > 0 else np.max(np.abs(test))
    aa = np.vdot(test, test).real
    bb = np.vdot(reference, reference).real
    mismatch = 0.0 if aa == 0 and bb == 0 else 1.0 - np.vdot(reference, test).real / np.sqrt(aa * bb)
    return error, mismatch

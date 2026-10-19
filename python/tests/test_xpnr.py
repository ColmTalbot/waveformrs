"""IMRPhenomXPNR with and without the parity fix of lalsuite !2593.

LALSimulation 6.2.0 and 6.2.1 refer the phase of the antisymmetric (2,2)
mode of IMRPhenomXPNR to the line of sight, which breaks the mirror
symmetry of the waveform; the fix (on lalsuite master, and the default here)
restores it. Agreement with 6.2.0 (``pnr_parity_fix=False``) is in
test_lalsimulation.py; the fixed model is compared with a patched
LALSimulation in the crate's tests (tests/xpnr.rs).
"""

import pickle

import numpy as np
import pytest

import waveformrs as wf

FREQUENCIES = np.geomspace(20.0, 1024.0, 400)
BINARIES = [
    dict(mass_1=36.0, mass_2=29.0, spin_1=(0.4, 0.2, 0.3), spin_2=(-0.2, 0.3, -0.1),
         luminosity_distance=500.0, inclination=0.7, phi_ref=1.2, f_ref=20.0),
    dict(mass_1=60.0, mass_2=15.0, spin_1=(0.6, -0.3, 0.5), spin_2=(0.1, 0.1, 0.0),
         luminosity_distance=500.0, inclination=2.3, phi_ref=4.0, f_ref=20.0),
    dict(mass_1=12.0, mass_2=6.0, spin_1=(0.2, 0.7, -0.1), spin_2=(0.3, -0.3, 0.6),
         luminosity_distance=500.0, inclination=0.3, phi_ref=5.5, f_ref=40.0),
]


def mirror(binary):
    """The mirror image in the orbital plane at f_ref (in-plane spins
    reversed), seen from the mirrored direction."""
    flip = lambda s: (-s[0], -s[1], s[2])  # noqa: E731
    return dict(binary, spin_1=flip(binary["spin_1"]), spin_2=flip(binary["spin_2"]),
                inclination=np.pi - binary["inclination"])


def mirror_violation(binary, **options):
    """max(|h+' - h+|, |hx' + hx|) / max |h| (zero for a symmetric model)."""
    hp, hc = wf.polarizations("IMRPhenomXPNR", FREQUENCIES, **binary, **options)
    hp_m, hc_m = wf.polarizations("IMRPhenomXPNR", FREQUENCIES, **mirror(binary), **options)
    scale = max(np.abs(hp).max(), np.abs(hc).max())
    return max(np.abs(hp_m - hp).max(), np.abs(hc_m + hc).max()) / scale


def test_option_defaults():
    assert wf.WaveformOptions().pnr_parity_fix
    assert not wf.WaveformOptions.lalsimulation_defaults().pnr_parity_fix
    options = wf.WaveformOptions(pnr_parity_fix=False)
    assert options.as_dict()["pnr_parity_fix"] is False
    assert pickle.loads(pickle.dumps(options)) == options
    assert options.replace(pnr_parity_fix=True) == wf.WaveformOptions()
    with pytest.raises(TypeError):
        wf.WaveformOptions(pnr_parity_fix="yes")


@pytest.mark.parametrize("index", range(len(BINARIES)))
def test_parity_fix_restores_the_mirror_symmetry(index):
    binary = BINARIES[index]
    assert mirror_violation(binary) < 1e-11
    assert mirror_violation(binary, pnr_parity_fix=False) > 1e-4


def test_fix_changes_only_xpnr():
    binary = BINARIES[0]
    fixed = wf.polarizations("IMRPhenomXPNR", FREQUENCIES, **binary)
    unfixed = wf.polarizations("IMRPhenomXPNR", FREQUENCIES, **binary, pnr_parity_fix=False)
    assert np.abs(fixed[0] - unfixed[0]).max() > 1e-4 * np.abs(unfixed[0]).max()
    for approximant in ["IMRPhenomXP", "IMRPhenomXPHM"]:
        a = wf.polarizations(approximant, FREQUENCIES, **binary)
        b = wf.polarizations(approximant, FREQUENCIES, **binary, pnr_parity_fix=False)
        assert all(np.array_equal(x, y) for x, y in zip(a, b))


@pytest.mark.parametrize("fix", [True, False])
def test_grid_and_batch(fix):
    """The grid path agrees with the sequence, the batch with single calls."""
    grid = wf.FrequencyGrid.uniform(20.0, 1024.0, 0.25)
    frequencies = np.asarray(grid.frequencies)
    hp, hc = wf.polarizations("IMRPhenomXPNR", frequencies, **BINARIES[0], pnr_parity_fix=fix)
    hp_grid, hc_grid = wf.polarizations_on_grid("IMRPhenomXPNR", grid, **BINARIES[0],
                                                pnr_parity_fix=fix)
    scale = np.abs(hp).max()
    assert np.abs(hp_grid - hp).max() < 1e-10 * scale
    assert np.abs(hc_grid - hc).max() < 1e-10 * scale

    parameters = {key: np.array([b[key] for b in BINARIES]) for key in BINARIES[0]}
    result = wf.batch_polarizations("IMRPhenomXPNR", grid, parameters, threads=3,
                                    pnr_parity_fix=fix)
    assert result.ok.all()
    for row, binary in enumerate(BINARIES):
        single = wf.polarizations_on_grid("IMRPhenomXPNR", grid, **binary, pnr_parity_fix=fix)
        assert np.array_equal(result.hplus[row], single[0])
        assert np.array_equal(result.hcross[row], single[1])

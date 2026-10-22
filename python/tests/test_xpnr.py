"""IMRPhenomXPNR with the parity fix of lalsuite !2593.

LALSimulation 6.2.0 and 6.2.1 refer the phase of the antisymmetric (2,2)
mode of IMRPhenomXPNR to the line of sight, which breaks the mirror
symmetry of the waveform; the fix (on lalsuite master, always applied here)
restores it. The model is compared with a patched LALSimulation in the
crate's tests (tests/xpnr.rs).
"""

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


def mirror_violation(binary):
    """max(|h+' - h+|, |hx' + hx|) / max |h| (zero for a symmetric model)."""
    hp, hc = wf.polarizations("IMRPhenomXPNR", FREQUENCIES, **binary)
    hp_m, hc_m = wf.polarizations("IMRPhenomXPNR", FREQUENCIES, **mirror(binary))
    scale = max(np.abs(hp).max(), np.abs(hc).max())
    return max(np.abs(hp_m - hp).max(), np.abs(hc_m + hc).max()) / scale


def test_the_option_is_gone():
    assert "pnr_parity_fix" not in wf.WaveformOptions().as_dict()
    with pytest.raises(TypeError):
        wf.WaveformOptions(pnr_parity_fix=False)


@pytest.mark.parametrize("index", range(len(BINARIES)))
def test_mirror_symmetry(index):
    assert mirror_violation(BINARIES[index]) < 1e-11


def test_grid_and_batch():
    """The grid path agrees with the sequence, the batch with single calls."""
    grid = wf.FrequencyGrid.uniform(20.0, 1024.0, 0.25)
    frequencies = np.asarray(grid.frequencies)
    hp, hc = wf.polarizations("IMRPhenomXPNR", frequencies, **BINARIES[0])
    hp_grid, hc_grid = wf.polarizations_on_grid("IMRPhenomXPNR", grid, **BINARIES[0])
    scale = np.abs(hp).max()
    assert np.abs(hp_grid - hp).max() < 1e-10 * scale
    assert np.abs(hc_grid - hc).max() < 1e-10 * scale

    parameters = {key: np.array([b[key] for b in BINARIES]) for key in BINARIES[0]}
    result = wf.batch_polarizations("IMRPhenomXPNR", grid, parameters, threads=3)
    assert result.ok.all()
    for row, binary in enumerate(BINARIES):
        single = wf.polarizations_on_grid("IMRPhenomXPNR", grid, **binary)
        assert np.array_equal(result.hplus[row], single[0])
        assert np.array_equal(result.hcross[row], single[1])

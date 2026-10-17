#!/usr/bin/env python
"""
Generate SpinTaylorT4 orbit evolutions from LALSimulation.

The output ``tests/data/spintaylor.json`` is consumed by the Rust integration
test ``tests/spintaylor.rs``.  Each case holds the inputs of
``SimInspiralSpinTaylorPNEvolveOrbit``, the exact length of the returned time
series and the 14 series (V, Phi, S1, S2, LNhat, E1) at a decimated set of
sample indices.

The cases mimic the calls IMRPhenomXPHM makes for its SpinTaylor Euler angles
(PrecVersion 320, coarse factor 10): a coarse integration with
``deltaT = 0.5 * 10 / fCut`` from f_min to ``fCut``, or from f_ref backwards to
``f_min - 0.5`` and forwards to ``fCut``, followed by a fine integration with
``deltaT = 0.5 / fCut`` from the transition sample, ten samples before the end
of the coarse integration.  ``fCut = (fRING + 8 fDAMP) / M`` is approximated
by ``0.2 / M``.

Usage
-----
    python scripts/generate_spintaylor_reference.py [--n-samples 20]

Requires ``python-lalsimulation`` (tested with lalsimulation 6.2.0).
"""
import argparse
import json
from pathlib import Path

import lal
import lalsimulation as lalsim
import numpy as np

DATA_DIR = Path(__file__).resolve().parent.parent / "tests" / "data"
SERIES = [
    "v", "phi", "s1_x", "s1_y", "s1_z", "s2_x", "s2_y", "s2_z",
    "lnhat_x", "lnhat_y", "lnhat_z", "e1_x", "e1_y", "e1_z",
]


def evolve(case):
    output = lalsim.SimInspiralSpinTaylorPNEvolveOrbit(
        case["delta_t"],
        case["m1_si"],
        case["m2_si"],
        case["f_start"],
        case["f_end"],
        *case["spin_1"],
        *case["spin_2"],
        *case["lnhat"],
        *case["e1"],
        *case["lambda"],
        *case["quadparam"],
        case["spin_order"],
        case["tidal_order"],
        case["phase_order"],
        int(case["lscorr"]),
        lalsim.SpinTaylorT4,
    )
    return np.array([series.data.data for series in output])


def make_case(label, m1, m2, spin_1, spin_2, delta_t, f_start, f_end, **kwargs):
    case = dict(
        label=label,
        m1_si=m1 * lal.MSUN_SI,
        m2_si=m2 * lal.MSUN_SI,
        spin_1=list(spin_1),
        spin_2=list(spin_2),
        lnhat=[0.0, 0.0, 1.0],
        e1=[1.0, 0.0, 0.0],
        delta_t=delta_t,
        f_start=f_start,
        f_end=f_end,
        **{"lambda": [0.0, 0.0]},
        quadparam=[1.0, 1.0],
        phase_order=7,
        spin_order=6,
        tidal_order=12,
        lscorr=False,
    )
    case.update(kwargs)
    return case


def f_cut(total_mass):
    return 0.2 / (total_mass * lal.MTSUN_SI)


def fine_case(coarse, data):
    """XPHM's high-frequency integration from the transition sample."""
    length = data.shape[1]
    index = length - 1 - min(9, length - 1)
    total_mass_si = coarse["m1_si"] + coarse["m2_si"]
    pi_gm = np.pi * total_mass_si * (lal.G_SI / lal.C_SI) / (lal.C_SI * lal.C_SI)
    f_trans = data[0, index] ** 3 / pi_gm
    return dict(
        coarse,
        label=coarse["label"] + ", fine",
        delta_t=0.5 / coarse["f_end"],
        f_start=f_trans,
        spin_1=data[2:5, index].tolist(),
        spin_2=data[5:8, index].tolist(),
        lnhat=data[8:11, index].tolist(),
        e1=data[11:14, index].tolist(),
    )


def cases():
    out = list()
    # (label, m1, m2, spin_1, spin_2, f_min, f_ref)
    configs = [
        ("36+29, f_ref = f_min", 36.0, 29.0, [0.4, 0.2, 0.3], [-0.2, 0.3, -0.1], 20.0, 20.0),
        ("60+15, f_ref > f_min", 60.0, 15.0, [0.6, -0.3, 0.5], [0.1, 0.1, 0.0], 15.0, 20.0),
        ("q = 1/20, maximal spins", 95.2380952380952, 4.76190476190476,
         [0.7, -0.6, 0.36], [-0.6, 0.6, -0.5], 10.0, 10.0),
        ("M = 5, q = 1/2", 3.33333333333333, 1.66666666666667,
         [0.0, 0.5, -0.6], [0.3, 0.0, 0.2], 20.0, 25.0),
        ("M = 300, anti-aligned", 200.0, 100.0, [0.5, 0.0, -0.8], [0.0, -0.3, -0.9], 5.0, 10.0),
    ]
    for label, m1, m2, spin_1, spin_2, f_min, f_ref in configs:
        fc = f_cut(m1 + m2)
        dt_coarse = 0.5 * 10 / fc
        if f_ref > f_min:
            out.append(make_case(label + ", backward", m1, m2, spin_1, spin_2, dt_coarse, f_ref, f_min - 0.5))
        coarse = make_case(label + ", forward", m1, m2, spin_1, spin_2, dt_coarse, f_ref, fc)
        out.append(coarse)
        out.append(fine_case(coarse, evolve(coarse)))
    # Orders as IMRPhenomXPHM passes them (-1), lscorr, matter, fEnd = 0
    out.append(make_case(
        "all orders, lscorr", 25.0, 12.0, [0.3, -0.6, 0.2], [0.5, 0.1, -0.4], 1 / 1024, 15.0, 600.0,
        phase_order=-1, spin_order=-1, tidal_order=-1, lscorr=True,
    ))
    out.append(make_case(
        "neutron stars, fEnd = 0", 1.6, 1.3, [0.05, 0.0, 0.02], [0.0, 0.03, 0.01], 1 / 4096, 60.0, 0.0,
        quadparam=[5.0, 7.0], **{"lambda": [400.0, 800.0]},
    ))
    return out


def decimate(length, n_samples):
    if length <= n_samples:
        return list(range(length))
    indices = set(range(3)) | set(range(length - 3, length))
    indices |= set(np.linspace(0, length - 1, n_samples - 6).round().astype(int).tolist())
    return sorted(indices)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-samples", type=int, default=20)
    args = parser.parse_args()
    output = list()
    for case in cases():
        data = evolve(case)
        length = data.shape[1]
        indices = decimate(length, args.n_samples)
        case = dict(case, length=length, indices=indices)
        for name, series in zip(SERIES, data):
            case[name] = series[indices].tolist()
        output.append(case)
        print(f"{case['label']:45s} {length:7d} samples")
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    path = DATA_DIR / "spintaylor.json"
    with open(path, "w") as ff:
        json.dump(
            dict(
                lalsimulation_version=lalsim.__version__,
                lal_version=lal.__version__,
                cases=output,
            ),
            ff,
        )
    print(f"Wrote {len(output)} cases to {path}")


if __name__ == "__main__":
    main()

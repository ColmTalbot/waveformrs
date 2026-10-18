#!/usr/bin/env python
"""
Compare two directories of reference data (scripts/generate_reference_data.py
and scripts/generate_spintaylor_reference.py output), e.g. the data of the
Linux lalsimulation with the committed osx-arm64 data in tests/data: how far
lalsimulation itself moves between platforms.

    compare_reference.py DIR REFERENCE_DIR [--label NAME --output FILE]

For each suite this prints (and with --output stores under NAME) the largest
error in the measures of the Rust tests: for frequency-domain waveforms the
relative error max|a - b| / max|b| and the mismatch with a flat weight
(tests/common/mod.rs), for time series max|a - b| / max|b|, for the angles
the absolute difference; and the fraction of bit-identical numbers.
"""
import argparse
import json
from pathlib import Path

import numpy as np


def weights(f):
    lo = np.concatenate([f[:1], f[:-1]])
    hi = np.concatenate([f[1:], f[-1:]])
    return 0.5 * (hi - lo)


def fd_compare(f, a, b):
    w = weights(np.asarray(f))
    ab = np.sum(w * (a * np.conj(b)).real)
    aa = np.sum(w * np.abs(a) ** 2)
    bb = np.sum(w * np.abs(b) ** 2)
    mismatch = 0.0 if aa == 0 and bb == 0 else float(1 - ab / np.sqrt(aa * bb))
    scale = np.max(np.abs(b))
    rel = float(np.max(np.abs(a - b)) / scale) if scale > 0 else float(np.max(np.abs(a - b)))
    return rel, mismatch


def relative(a, b):
    a, b = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
    if a.shape != b.shape:
        return float("inf")
    scale = np.max(np.abs(b)) if b.size else 0.0
    diff = np.max(np.abs(a - b)) if b.size else 0.0
    return float(diff / scale) if scale > 0 else float(diff)


class Tally:
    def __init__(self):
        self.rel = 0.0
        self.mismatch = 0.0
        self.identical = 0
        self.total = 0
        self.notes = []

    def bits(self, a, b):
        a, b = np.asarray(a, dtype=float).ravel(), np.asarray(b, dtype=float).ravel()
        if a.shape != b.shape:
            self.notes.append("length differs")
            return
        self.identical += int(np.sum(a.view(np.uint64) == b.view(np.uint64)))
        self.total += a.size

    def result(self):
        return dict(
            max_relative_error=self.rel,
            max_mismatch=self.mismatch,
            fraction_bit_identical=self.identical / self.total if self.total else None,
            notes=sorted(set(self.notes)),
        )


def complex_series(case, prefix):
    return np.asarray(case[f"{prefix}_real"]) + 1j * np.asarray(case[f"{prefix}_imag"])


def compare_suite(name, test, ref):
    tally = Tally()
    for a, b in zip(test["cases"], ref["cases"]):
        if "hplus_real" in b:
            for pol in ("hplus", "hcross"):
                x, y = complex_series(a, pol), complex_series(b, pol)
                rel, mm = fd_compare(b["frequencies"], x, y)
                tally.rel, tally.mismatch = max(tally.rel, rel), max(tally.mismatch, mm)
                tally.bits(np.concatenate([x.real, x.imag]), np.concatenate([y.real, y.imag]))
        elif name == "imrphenomxhm_modes":
            for ma, mb in zip(a["modes"], b["modes"]):
                x = np.asarray(ma["real"]) + 1j * np.asarray(ma["imag"])
                y = np.asarray(mb["real"]) + 1j * np.asarray(mb["imag"])
                rel, mm = fd_compare(b["frequencies"], x, y)
                tally.rel, tally.mismatch = max(tally.rel, rel), max(tally.mismatch, mm)
                tally.bits(np.concatenate([x.real, x.imag]), np.concatenate([y.real, y.imag]))
        elif name == "choose_fd_modes":
            for ma, mb in zip(a["modes"], b["modes"]):
                for half in ("positive", "negative"):
                    x = np.asarray(ma[f"{half}_real"]) + 1j * np.asarray(ma[f"{half}_imag"])
                    y = np.asarray(mb[f"{half}_real"]) + 1j * np.asarray(mb[f"{half}_imag"])
                    rel, mm = fd_compare(np.arange(len(y), dtype=float), x, y)
                    tally.rel, tally.mismatch = max(tally.rel, rel), max(tally.mismatch, mm)
                    tally.bits(np.concatenate([x.real, x.imag]),
                               np.concatenate([y.real, y.imag]))
        elif name == "imrphenomxp_spintaylor_angles":
            if a["length"] != b["length"]:
                tally.notes.append("length differs")
            for key in ("alpha", "cos_beta", "gamma"):
                tally.rel = max(tally.rel, float(np.max(np.abs(np.subtract(a[key], b[key])))))
                tally.bits(a[key], b[key])
        elif name == "spintaylor":
            if a["length"] != b["length"]:
                tally.notes.append("length differs")
                continue
            for key in ("v", "phi", "s1_x", "s1_y", "s1_z", "s2_x", "s2_y", "s2_z",
                        "lnhat_x", "lnhat_y", "lnhat_z", "e1_x", "e1_y", "e1_z"):
                tally.rel = max(tally.rel, relative(a[key], b[key]))
                tally.bits(a[key], b[key])
        elif name in ("nrsur7dq4", "nrsur7dq4_modes"):
            if a.get("lal_error") or b.get("lal_error"):
                if bool(a.get("lal_error")) != bool(b.get("lal_error")):
                    tally.notes.append("error differs")
                continue
            if a.get("epoch") != b.get("epoch"):
                tally.notes.append("epoch differs")
            series = (["hplus", "hcross"] if name == "nrsur7dq4" else [])
            for key in series:
                tally.rel = max(tally.rel, relative(a[key], b[key]))
                tally.bits(a[key], b[key])
            for ma, mb in zip(a.get("modes", []), b.get("modes", [])):
                for key in ("real", "imag"):
                    tally.rel = max(tally.rel, relative(ma[key], mb[key]))
                    tally.bits(ma[key], mb[key])
        elif name == "nrsur7dq4_dynamics":
            for key in ("q0", "qx", "qy", "qz", "orbphase", "chiAx", "chiAy", "chiAz",
                        "chiBx", "chiBy", "chiBz"):
                tally.rel = max(tally.rel, relative(a[key], b[key]))
                tally.bits(a[key], b[key])
        elif name == "spin_weighted_harmonics":
            # the angles are partly searched for on the generating platform:
            # compare the common ones (21 harmonics per angle, l <= 4)
            def by_angle(case):
                values = np.array([v[2:] for v in case["values"]])
                return {(t, p): values[21 * i:21 * (i + 1)]
                        for i, (t, p) in enumerate(zip(case["theta"], case["phi"]))}
            x, y = by_angle(a), by_angle(b)
            common = [key for key in y if key in x]
            tally.notes.append(f"{len(common)} of {len(y)} angles in both")
            if not common:
                continue
            x = np.array([x[key] for key in common])
            y = np.array([y[key] for key in common])
            tally.rel = max(tally.rel, float(np.max(np.abs(x - y))))
            tally.bits(x, y)
        else:
            tally.notes.append("not compared")
    return tally.result()


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("directory")
    parser.add_argument("reference")
    parser.add_argument("--label")
    parser.add_argument("--output")
    args = parser.parse_args()
    results = dict()
    for path in sorted(Path(args.reference).glob("*.json")):
        other = Path(args.directory) / path.name
        if not other.exists() or path.stem == "nrsur7dq4_hdf5_checksums":
            continue
        results[path.stem] = compare_suite(
            path.stem, json.loads(other.read_text()), json.loads(path.read_text()))
    print(f"{'suite':<36}{'max rel/abs':>13}{'mismatch':>11}{'identical':>11}")
    for name, r in results.items():
        frac = r["fraction_bit_identical"]
        print(f"{name:<36}{r['max_relative_error']:13.2e}{r['max_mismatch']:11.2e}"
              f"{'' if frac is None else f'{100 * frac:10.1f}%'} {' '.join(r['notes'])}")
    if args.output:
        output = Path(args.output)
        combined = json.loads(output.read_text()) if output.exists() else dict()
        combined[args.label] = results
        output.write_text(json.dumps(combined, indent=1) + "\n")


if __name__ == "__main__":
    main()

"""
Accuracy of the original waveformrs main (87f9d9b) against lalsimulation:
mismatch maximised over time and phase (aLIGO design PSD) and the maximum
relative amplitude error, over 100 random aligned-spin draws per model.

Usage
-----
    git worktree add /tmp/baseline 87f9d9b
    # save the program below as /tmp/baseline/examples/dump.rs
    (cd /tmp/baseline && cargo build --release --example dump)
    OMP_NUM_THREADS=1 python report/baseline_accuracy.py \
        /tmp/baseline/target/release/examples/dump report/data/baseline_main_accuracy.json

examples/dump.rs
----------------
    // Print h(f) = waveform(f, phi_c = 0) for the original implementation.
    // usage: dump <TaylorF2|IMRPhenomD> M q chi1 chi2 f_min f_max delta_f
    use waveformrs::{imrphenomd::IMRPhenomD, taylor::TaylorF2, waveform::Waveform};
    fn main() {
        let a: Vec<String> = std::env::args().collect();
        let p: Vec<f64> = a[2..].iter().map(|x| x.parse().unwrap()).collect();
        let (m, q, c1, c2, fmin, fmax, df) = (p[0], p[1], p[2], p[3], p[4], p[5], p[6]);
        let n = ((fmax - fmin) / df).round() as usize + 1;
        let f: Vec<f64> = (0..n).map(|i| fmin + i as f64 * df).collect();
        let h = match a[1].as_str() {
            "TaylorF2" => TaylorF2::new(m, q, c1, c2, 100.0).waveform(&f, 0.0),
            _ => IMRPhenomD::new(m, q, c1, c2, 100.0).waveform(&f, 0.0),
        };
        for x in h { println!("{:e} {:e}", x.re, x.im); }
    }
"""
import json
import subprocess
import sys

import lal
import lalsimulation as ls
import numpy as np
from scipy.optimize import minimize_scalar

DUMP = sys.argv[1]
OUT = sys.argv[2]
DF = 0.125
F_MIN = 20.0
MF_MAX = dict(TaylorF2=6**-1.5 / np.pi, IMRPhenomD=0.2)
Q_MIN = dict(TaylorF2=1 / 20, IMRPhenomD=1 / 18)


def old(model, m, q, c1, c2, f_max):
    out = subprocess.run(
        [DUMP, model, *map(repr, [m, q, c1, c2, F_MIN, f_max, DF])],
        capture_output=True, text=True, check=True,
    ).stdout.split()
    values = np.array(out, dtype=float).reshape(-1, 2)
    return values[:, 0] + 1j * values[:, 1]


def lal_h(model, m, q, c1, c2, freqs):
    m1 = m / (1 + q)
    m2 = m - m1
    fv = lal.CreateREAL8Vector(len(freqs))
    fv.data = freqs
    hp, _ = ls.SimInspiralChooseFDWaveformSequence(
        0.0, m1 * lal.MSUN_SI, m2 * lal.MSUN_SI, 0, 0, c1, 0, 0, c2, F_MIN,
        100e6 * lal.PC_SI, 0.0, None, ls.GetApproximantFromString(model), fv,
    )
    return hp.data.data.copy()


def maximised_match(a, b, freqs, weights):
    z = a * np.conj(b) * weights
    norm = np.sqrt(np.sum(np.abs(a) ** 2 * weights) * np.sum(np.abs(b) ** 2 * weights))
    n = 1 << int(np.ceil(np.log2(len(z) * 32)))
    coarse = np.abs(np.fft.ifft(z, n)) * n
    j = int(np.argmax(coarse))
    dt = 1 / (n * DF)
    t0 = j * dt

    def neg(t):
        return -np.abs(np.sum(z * np.exp(2j * np.pi * freqs * t)))

    res = minimize_scalar(neg, bounds=(t0 - dt, t0 + dt), method="bounded",
                          options=dict(xatol=1e-12))
    return max(-res.fun, coarse[j]) / norm


def main():
    rng = np.random.default_rng(1)
    results = dict()
    for model in ["TaylorF2", "IMRPhenomD"]:
        records = []
        while len(records) < 100:
            m = float(np.exp(rng.uniform(np.log(5), np.log(200))))
            q = float(rng.uniform(Q_MIN[model], 1))
            c1, c2 = (float(x) for x in rng.uniform(-0.9, 0.9, 2))
            f_max = min(MF_MAX[model] / (m * lal.MTSUN_SI), 1024.0)
            if f_max < 2 * F_MIN:
                continue
            f_max = float(F_MIN + np.floor((f_max - F_MIN) / DF) * DF)
            freqs = F_MIN + DF * np.arange(int(round((f_max - F_MIN) / DF)) + 1)
            h_old = old(model, m, q, c1, c2, f_max)
            h_lal = lal_h(model, m, q, c1, c2, freqs)
            weights = 1 / np.array([ls.SimNoisePSDaLIGOZeroDetHighPower(f) for f in freqs])
            match = maximised_match(h_old, h_lal, freqs, weights)
            amp_err = np.max(np.abs(np.abs(h_old) - np.abs(h_lal))) / np.max(np.abs(h_lal))
            records.append(dict(total_mass=m, mass_ratio=q, chi_1=c1, chi_2=c2,
                                mismatch=1 - match, amplitude_error=float(amp_err)))
        mm = np.array([r["mismatch"] for r in records])
        ae = np.array([r["amplitude_error"] for r in records])
        print(f"{model}: mismatch median {np.median(mm):.2e} max {mm.max():.2e} "
              f"(>1e-3: {np.sum(mm > 1e-3)}); amplitude error median {np.median(ae):.2e} max {ae.max():.2e}")
        results[model] = records
    json.dump(results, open(OUT, "w"))


if __name__ == "__main__":
    main()

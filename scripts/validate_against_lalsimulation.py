#!/usr/bin/env python
"""
Compare waveformrs against lalsimulation over random parameter draws.

For each approximant this draws random binaries, evaluates both codes on the
same uniform frequency grid, and records the noise-weighted mismatch (aLIGO
design sensitivity, no maximisation over time or phase), the maximum
relative error and the wall-clock time of both implementations.

Usage
-----
    OMP_NUM_THREADS=1 python scripts/validate_against_lalsimulation.py \
        --approximants TaylorF2 IMRPhenomD --n-draws 200 --outdir report/data

The Rust side is driven through ``cargo run --release --example lal_validation``.

``OMP_NUM_THREADS=1`` is required for reproducible MSA angles (IMRPhenomXP and
IMRPhenomXPHM with PrecVersion 223, and the SpinTaylor angles where they fall
back to MSA):
lalsimulation evaluates them inside an OpenMP loop that races on shared
state.

When the SpinTaylor integration fails through an ``XLAL_ERROR``, the C
waveform drivers fall back to the MSA angles and return success, but leave
``xlalErrno`` set, so the Python wrappers raise "Internal function call
failed". For such cases (identified by lalsimulation's fallback warning) the
reference is the PrecVersion 223 waveform, which the C-level fallback
reproduces bit for bit; they are counted as ``n_lal_errno_fallback``.

``IMRPhenomXPNR`` compares the model without the fix of lalsuite merge
request !2593 (``with_parity_fix(false)``) with lalsimulation 6.2.0.
``IMRPhenomXPNR_parityfix`` compares the model with the fix (its default)
with a liblalsimulation built with the fix
(``scripts/build_patched_lalsimulation.sh``), given with
``--lalsimulation-library`` and called through ctypes
(``scripts/lalsimulation_ctypes.py``).
"""
import argparse
import contextlib
import io
import json
import subprocess
import time
from pathlib import Path

import lal
import lalsimulation as lalsim
import numpy as np

ROOT = Path(__file__).resolve().parent.parent

PRECESSING = {
    "IMRPhenomXP",
    "IMRPhenomXP_NNLO",
    "IMRPhenomXP_MSA",
    "IMRPhenomXPHM",
    "IMRPhenomXPHM_NNLO",
    "IMRPhenomXPHM_MSA",
    "IMRPhenomXPNR",
    "IMRPhenomXPNR_parityfix",
}

# approximants that need a patched lalsimulation (--lalsimulation-library)
PATCHED = {"IMRPhenomXPNR_parityfix"}

# Approximant -> (dimensionless cutoff frequency, extra lal flags, max mass
# ratio, optionally the lalsimulation approximant if the name differs)
SETTINGS = dict(
    TaylorF2=dict(mf_max=6**-1.5 / np.pi, flags=dict(), q_min=1 / 20),
    IMRPhenomD=dict(mf_max=0.2, flags=dict(), q_min=1 / 18),
    IMRPhenomXAS=dict(mf_max=0.3, flags=dict(), q_min=1 / 20),
    IMRPhenomXHM=dict(
        mf_max=0.3, flags=dict(PhenomXHMThresholdMband=0.0), q_min=1 / 20
    ),
    IMRPhenomXP=dict(
        mf_max=0.3, flags=dict(PhenomXPrecVersion=320), q_min=1 / 20
    ),
    IMRPhenomXP_NNLO=dict(
        mf_max=0.3,
        flags=dict(PhenomXPrecVersion=102),
        q_min=1 / 20,
        lal_approximant="IMRPhenomXP",
    ),
    IMRPhenomXP_MSA=dict(
        mf_max=0.3,
        flags=dict(PhenomXPrecVersion=223),
        q_min=1 / 20,
        lal_approximant="IMRPhenomXP",
    ),
    IMRPhenomXPHM=dict(
        mf_max=0.3,
        flags=dict(
            PhenomXPrecVersion=320,
            PhenomXHMThresholdMband=0.0,
            PhenomXPHMThresholdMband=0.0,
        ),
        q_min=1 / 20,
    ),
    IMRPhenomXPHM_NNLO=dict(
        mf_max=0.3,
        flags=dict(
            PhenomXPrecVersion=102,
            PhenomXHMThresholdMband=0.0,
            PhenomXPHMThresholdMband=0.0,
        ),
        q_min=1 / 20,
        lal_approximant="IMRPhenomXPHM",
    ),
    IMRPhenomXPHM_MSA=dict(
        mf_max=0.3,
        flags=dict(
            PhenomXPrecVersion=223,
            PhenomXHMThresholdMband=0.0,
            PhenomXPHMThresholdMband=0.0,
        ),
        q_min=1 / 20,
        lal_approximant="IMRPhenomXPHM",
    ),
    # SimInspiralChooseFDWaveformSequence switches multibanding off
    IMRPhenomXPNR=dict(mf_max=0.3, flags=dict(), q_min=1 / 20),
    IMRPhenomXPNR_parityfix=dict(
        mf_max=0.3, flags=dict(), q_min=1 / 20, lal_approximant="IMRPhenomXPNR"
    ),
)


def lal_dict(flags):
    params = lal.CreateDict()
    for key, value in flags.items():
        getattr(lalsim, f"SimInspiralWaveformParamsInsert{key}")(params, value)
    return params


def draw(approximant, rng, delta_f, f_high):
    settings = SETTINGS[approximant]
    while True:
        total_mass = np.exp(rng.uniform(np.log(5.0), np.log(300.0)))
        mass_ratio = rng.uniform(settings["q_min"], 1.0)
        mass_1 = total_mass / (1 + mass_ratio)
        mass_2 = total_mass - mass_1
        f_min = rng.choice([10.0, 15.0, 20.0])
        f_max = min(settings["mf_max"] / (total_mass * lal.MTSUN_SI), f_high)
        if f_max > 2 * f_min:
            break
    if approximant in PRECESSING:
        spins = list()
        for _ in range(2):
            magnitude = rng.uniform(0, 0.99)
            cos_theta = rng.uniform(-1, 1)
            phi = rng.uniform(0, 2 * np.pi)
            sin_theta = np.sqrt(1 - cos_theta**2)
            spins.append(
                list(
                    magnitude
                    * np.array(
                        [sin_theta * np.cos(phi), sin_theta * np.sin(phi), cos_theta]
                    )
                )
            )
    else:
        spins = [[0.0, 0.0, rng.uniform(-0.99, 0.99)] for _ in range(2)]
    f_ref = f_min if rng.uniform() < 0.5 else rng.uniform(f_min, 3 * f_min)
    n_lo = int(np.ceil(f_min / delta_f))
    n_hi = int(np.floor(f_max / delta_f))
    frequencies = np.arange(n_lo, n_hi + 1) * delta_f
    return dict(
        approximant=approximant,
        mass_1=mass_1,
        mass_2=mass_2,
        spin_1=spins[0],
        spin_2=spins[1],
        luminosity_distance=float(np.exp(rng.uniform(np.log(100), np.log(3000)))),
        inclination=float(np.arccos(rng.uniform(-1, 1))),
        phi_ref=float(rng.uniform(0, 2 * np.pi)),
        f_ref=float(f_ref),
        f_min=float(frequencies[0]),
        frequencies=frequencies,
    )


def reports_msa_fallback(call):
    """Rerun a failed lalsimulation call capturing its XLAL warnings: whether
    it fell back from the SpinTaylor to the MSA angles."""
    level = lal.GetDebugLevel()
    lal.ClobberDebugLevel(level | lal.LALWARNING)
    redirect = lal.swig_redirect_standard_output_error(True)
    buffer = io.StringIO()
    try:
        with contextlib.redirect_stderr(buffer), contextlib.redirect_stdout(buffer):
            call()
    except RuntimeError:
        pass
    finally:
        lal.swig_redirect_standard_output_error(redirect)
        lal.ClobberDebugLevel(level)
    return "default to MSA angles" in buffer.getvalue()


def lal_waveform(case, flags=None):
    settings = SETTINGS[case["approximant"]]
    frequencies = lal.CreateREAL8Vector(len(case["frequencies"]))
    frequencies.data = case["frequencies"]
    params = lal_dict(settings["flags"] if flags is None else flags)
    approximant = lalsim.GetApproximantFromString(
        settings.get("lal_approximant", case["approximant"])
    )

    def call():
        return lalsim.SimInspiralChooseFDWaveformSequence(
            case["phi_ref"],
            case["mass_1"] * lal.MSUN_SI,
            case["mass_2"] * lal.MSUN_SI,
            *case["spin_1"],
            *case["spin_2"],
            case["f_ref"],
            case["luminosity_distance"] * 1e6 * lal.PC_SI,
            case["inclination"],
            params,
            approximant,
            frequencies,
        )

    try:
        hplus, hcross = call()
    except RuntimeError as error:
        version = settings["flags"].get("PhenomXPrecVersion")
        if flags is None and version == 320 and reports_msa_fallback(call):
            hp, hc, timing = lal_waveform(
                case, flags=dict(settings["flags"], PhenomXPrecVersion=223)
            )
            if hp is not None:
                case["lal_errno_fallback"] = True
            return hp, hc, timing
        return None, None, str(error)
    timings = list()
    for _ in range(5):
        start = time.perf_counter()
        call()
        timings.append(time.perf_counter() - start)
    return hplus.data.data.copy(), hcross.data.data.copy(), float(np.median(timings))


def rust_waveforms(cases):
    payload = json.dumps(
        dict(
            cases=[
                {**case, "frequencies": case["frequencies"].tolist()} for case in cases
            ]
        )
    )
    result = subprocess.run(
        ["cargo", "run", "--quiet", "--release", "--example", "lal_validation"],
        input=payload,
        capture_output=True,
        text=True,
        cwd=ROOT,
        check=True,
    )
    outputs = json.loads(result.stdout)["results"]
    return [
        (
            np.array(out["hplus_real"]) + 1j * np.array(out["hplus_imag"]),
            np.array(out["hcross_real"]) + 1j * np.array(out["hcross_imag"]),
            out["seconds"],
            out["error"],
        )
        for out in outputs
    ]


def psd(frequencies):
    values = np.array(
        [lalsim.SimNoisePSDaLIGOZeroDetHighPower(f) for f in frequencies]
    )
    return values


def mismatch(a, b, weights):
    inner = lambda x, y: np.sum((x * np.conj(y)).real * weights)
    aa, bb = inner(a, a), inner(b, b)
    if aa == 0 and bb == 0:
        return 0.0
    return float(1 - inner(a, b) / np.sqrt(aa * bb))


def max_relative_error(a, b):
    scale = np.max(np.abs(b))
    return float(np.max(np.abs(a - b)) / scale) if scale > 0 else 0.0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--approximants", nargs="+", default=["TaylorF2", "IMRPhenomD"])
    parser.add_argument("--n-draws", type=int, default=200)
    parser.add_argument("--delta-f", type=float, default=0.25)
    parser.add_argument("--f-high", type=float, default=2048.0)
    parser.add_argument("--seed", type=int, default=170817)
    parser.add_argument("--outdir", default=str(ROOT / "report" / "data"))
    parser.add_argument(
        "--lalsimulation-library",
        help="a liblalsimulation with the fix of lalsuite !2593 for IMRPhenomXPNR_parityfix",
    )
    args = parser.parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    global lalsim
    swig_lalsim = lalsim
    patched = None
    if args.lalsimulation_library is not None:
        from lalsimulation_ctypes import Library, Overlay

        library = Library(args.lalsimulation_library)
        patched = Overlay(swig_lalsim, library)
        patched_build = dict(library=Path(library.path).name, **library.vcs_info)

    for approximant in args.approximants:
        if approximant in PATCHED:
            if patched is None:
                raise SystemExit(f"{approximant} needs --lalsimulation-library")
            lalsim = patched
        else:
            lalsim = swig_lalsim
        rng = np.random.default_rng(args.seed)
        cases = [
            draw(approximant, rng, args.delta_f, args.f_high)
            for _ in range(args.n_draws)
        ]
        rust = rust_waveforms(cases)
        records = list()
        failures = list()
        for case, (hp_rs, hc_rs, t_rs, error) in zip(cases, rust):
            hp_lal, hc_lal, t_lal = lal_waveform(case)
            errno_fallback = case.pop("lal_errno_fallback", False)
            if error is not None or hp_lal is None:
                # Cases where lalsimulation silently falls back to angles that
                # are not implemented here (or either code raises an error).
                failures.append(
                    dict(
                        **{
                            key: value
                            for key, value in case.items()
                            if key not in ("frequencies", "approximant")
                        },
                        rust_error=error,
                        lal_error=t_lal if hp_lal is None else None,
                    )
                )
                continue
            weights = 1 / psd(case["frequencies"])
            records.append(
                dict(
                    **{
                        key: value
                        for key, value in case.items()
                        if key not in ("frequencies", "approximant")
                    },
                    n_frequencies=len(case["frequencies"]),
                    mismatch_plus=mismatch(hp_rs, hp_lal, weights),
                    mismatch_cross=mismatch(hc_rs, hc_lal, weights),
                    max_rel_err_plus=max_relative_error(hp_rs, hp_lal),
                    max_rel_err_cross=max_relative_error(hc_rs, hc_lal),
                    time_rust=t_rs,
                    time_lal=t_lal,
                    lal_errno_fallback=errno_fallback,
                )
            )
        mismatches = np.array(
            [max(r["mismatch_plus"], r["mismatch_cross"]) for r in records]
        )
        errors = np.array(
            [max(r["max_rel_err_plus"], r["max_rel_err_cross"]) for r in records]
        )
        speedup = np.array([r["time_lal"] / r["time_rust"] for r in records])
        summary = dict(
            approximant=approximant,
            lalsimulation_version=lalsim.__version__,
            n_draws=len(records) + len(failures),
            n_compared=len(records),
            n_lal_errno_fallback=sum(r["lal_errno_fallback"] for r in records),
            failures=failures,
            max_mismatch=float(mismatches.max()),
            median_mismatch=float(np.median(mismatches)),
            max_relative_error=float(errors.max()),
            median_relative_error=float(np.median(errors)),
            median_speedup=float(np.median(speedup)),
            records=records,
        )
        if approximant in PATCHED:
            summary["lalsimulation_build"] = patched_build
        with open(outdir / f"validation_{approximant}.json", "w") as ff:
            json.dump(summary, ff)
        print(
            f"{approximant:>14}: max mismatch {mismatches.max():.2e}, "
            f"median {np.median(mismatches):.2e}, max rel err {errors.max():.2e}, "
            f"median speed-up x{np.median(speedup):.2f}, "
            f"{len(failures)} not compared, "
            f"{summary['n_lal_errno_fallback']} via the 223 reference"
        )


if __name__ == "__main__":
    main()

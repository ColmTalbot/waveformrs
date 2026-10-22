//! Long, low-mass waveforms against LALSimulation.
//!
//! The reference data (`tests/data/long_waveforms.json`, from
//! `scripts/generate_long_reference.py`) cover every frequency-domain model
//! on grids of up to ~520k frequencies: `Delta f = 1/128 Hz` from 20 Hz to
//! 1-4 kHz for binaries of 3-6 solar masses (signals of 50-150 s), plus a
//! 2.6 solar-mass binary from 10 Hz with `Delta f = 1/1024 Hz` (~20 min) and
//! a 10 solar-mass one with `Delta f = 1/32 Hz`. Each case stores the
//! definition of its grid and LAL's polarizations at a subset of ~1250
//! indices (every k-th point, the band edges and log-spaced points at low
//! frequency). The waveforms are evaluated on the whole grid here, since the
//! SpinTaylor angles depend on it, and compared at those indices.
//!
//! # Tolerances
//!
//! The error of a case is the largest `|h - h_LAL|` at the stored indices,
//! relative to the largest `|h_LAL|`, as in `tests/lalsimulation.rs`. The
//! inspiral phase of these binaries reaches ~1e4 rad at `f_min` (4e4 rad
//! from 10 Hz), where one unit in the last place (ulp) is 2e-12 to 7e-12 rad,
//! and the waveform is largest at `f_min`, so the relative error is about
//! the error of the phase at `f_min` in radians.
//!
//! - The per-frequency paths of IMRPhenomXHM, XP and XPHM (`polarizations`,
//!   `polarizations_uniform`) follow LAL's operations, including the
//!   SpinTaylor integration over ~2 million output steps, but form the
//!   powers of `Mf` from `cbrt` rather than `pow(Mf, 1/6)`.
//! - TaylorF2, IMRPhenomD and IMRPhenomXAS evaluate their powers of `Mf`
//!   per frequency from `cbrt` and products (and the polynomials in a
//!   different order), and every grid path forms them as products of powers
//!   stored in the [`FrequencyGrid`]. Their phases are rounded differently
//!   from LAL's. Against an exact (50-digit) evaluation of the inspiral phase
//!   of IMRPhenomXAS, LAL's own value is off by up to ~9 ulp of the largest
//!   phase term (and by +0.8 ulp on average, because `pow(Mf, 1.0 / 6.0)`
//!   uses the rounded exponent), the per-frequency path by up to ~8 ulp and
//!   the grid path by up to ~11 ulp. They differ from LAL by up to ~20 ulp,
//!   3e-11 relative from 20 Hz and 1.4e-10 from 10 Hz. Storing correctly
//!   rounded powers in the grid (e.g. `f^(-5/3)`) halves the grid path's own
//!   error but not its difference from LAL, which is then set by LAL's
//!   rounding. These paths are held to [`PHASE_ULPS`] ulp of the
//!   leading-order phase at `f_min`, 3-5 times the largest difference seen.
//!
//! The tests are skipped in debug builds, which take minutes for these grids;
//! run them with `cargo test --release` (a few seconds).
mod common;

use std::collections::HashMap;
use std::sync::OnceLock;

use num_complex::Complex;
use serde::Deserialize;
use waveformrs::{
    grid::FrequencyGrid,
    imrphenomd::IMRPhenomD,
    imrphenomx::{
        IMRPhenomXAS, IMRPhenomXHM, IMRPhenomXP, IMRPhenomXPHM, PrecessingBinary, PrecessionVersion,
    },
    taylor::TaylorF2,
    waveform::Waveform,
};

type C64 = Complex<f64>;

/// Units in the last place of the leading-order phase at `f_min` allowed for
/// the paths that round the phase differently from LAL.
const PHASE_ULPS: f64 = 64.0;

/// The solar mass in seconds, `LAL_MTSUN_SI`.
const MTSUN_SI: f64 = 4.925490947641267e-6;

#[derive(Debug, Deserialize)]
struct LongFile {
    lalsimulation_version: String,
    cases: Vec<LongCase>,
    angles: Vec<AngleCase>,
}

#[derive(Debug, Deserialize)]
struct LongCase {
    name: String,
    approximant: String,
    /// `sequence` (ChooseFDWaveformSequence on `f_min + k delta_f`) or
    /// `uniform` (ChooseFDWaveform on `k delta_f`).
    generator: String,
    lal_flags: HashMap<String, f64>,
    mass_1: f64,
    mass_2: f64,
    spin_1: [f64; 3],
    spin_2: [f64; 3],
    luminosity_distance: f64,
    inclination: f64,
    phi_ref: f64,
    f_ref: f64,
    f_min: f64,
    f_max: f64,
    delta_f: f64,
    length: usize,
    indices: Vec<usize>,
    hplus_real: Vec<f64>,
    hplus_imag: Vec<f64>,
    hcross_real: Vec<f64>,
    hcross_imag: Vec<f64>,
}

/// Every `stride`-th value of `SimIMRPhenomXPSpinTaylorAngles`.
#[derive(Debug, Deserialize)]
struct AngleCase {
    name: String,
    mass_1: f64,
    mass_2: f64,
    spin_1: [f64; 3],
    spin_2: [f64; 3],
    f_min: f64,
    f_max: f64,
    delta_f: f64,
    f_ref: f64,
    phi_ref: f64,
    stride: usize,
    length: usize,
    alpha: Vec<f64>,
    cos_beta: Vec<f64>,
    gamma: Vec<f64>,
}

impl LongCase {
    /// The frequencies passed to `ChooseFDWaveformSequence`.
    fn sequence(&self) -> Vec<f64> {
        assert_eq!(self.generator, "sequence");
        (0..self.length)
            .map(|k| self.f_min + k as f64 * self.delta_f)
            .collect()
    }

    fn version(&self) -> PrecessionVersion {
        match self.lal_flags.get("PhenomXPrecVersion").copied() {
            Some(102.0) => PrecessionVersion::Nnlo,
            Some(223.0) => PrecessionVersion::Msa,
            Some(320.0) => PrecessionVersion::SpinTaylor,
            other => panic!("unexpected PhenomXPrecVersion {other:?}"),
        }
    }

    fn binary(&self) -> PrecessingBinary {
        PrecessingBinary {
            mass_1: self.mass_1,
            mass_2: self.mass_2,
            spin_1: self.spin_1,
            spin_2: self.spin_2,
            luminosity_distance: self.luminosity_distance,
            inclination: self.inclination,
            phi_ref: self.phi_ref,
            f_ref: self.f_ref,
        }
    }

    fn label(&self) -> String {
        let version = self
            .lal_flags
            .get("PhenomXPrecVersion")
            .map_or(String::new(), |v| format!(" ({v})"));
        format!(
            "{}{version} {} {} [{} frequencies]",
            self.approximant, self.name, self.generator, self.length
        )
    }

    /// The leading-order (Newtonian) SPA phase at `f_min`,
    /// `3 / (128 eta) (pi M f_min)^(-5/3)`: the size of the largest terms
    /// of the inspiral phase.
    fn phase_scale(&self) -> f64 {
        let total_mass = self.mass_1 + self.mass_2;
        let eta = self.mass_1 * self.mass_2 / (total_mass * total_mass);
        let v = (std::f64::consts::PI * total_mass * MTSUN_SI * self.f_min).cbrt();
        3.0 / (128.0 * eta) / v.powi(5)
    }

    /// Tolerance of the paths that round the phase differently from LAL.
    fn rounding_tolerance(&self) -> f64 {
        PHASE_ULPS * f64::EPSILON * self.phase_scale()
    }

    /// Largest difference from LAL at the stored indices, relative to the
    /// largest stored LAL value, for each polarization, and the frequency
    /// at which it occurs.
    fn errors(&self, hplus: &[C64], hcross: &[C64]) -> [(f64, f64); 2] {
        assert_eq!(hplus.len(), self.length, "{}", self.label());
        assert_eq!(hcross.len(), self.length, "{}", self.label());
        let frequency = |k: usize| match self.generator.as_str() {
            "sequence" => self.f_min + k as f64 * self.delta_f,
            _ => k as f64 * self.delta_f,
        };
        let mut out = [(0.0, 0.0); 2];
        for (slot, test, re, im) in [
            (0, hplus, &self.hplus_real, &self.hplus_imag),
            (1, hcross, &self.hcross_real, &self.hcross_imag),
        ] {
            let mut max_ref: f64 = 0.0;
            let mut max_diff: f64 = 0.0;
            let mut at = 0.0;
            for (j, &k) in self.indices.iter().enumerate() {
                let truth = C64::new(re[j], im[j]);
                max_ref = max_ref.max(truth.norm());
                let diff = (test[k] - truth).norm();
                if diff > max_diff {
                    max_diff = diff;
                    at = frequency(k);
                }
            }
            out[slot] = (max_diff / max_ref, at);
        }
        out
    }
}

fn reference() -> &'static LongFile {
    static REFERENCE: OnceLock<LongFile> = OnceLock::new();
    REFERENCE.get_or_init(|| {
        let path = format!(
            "{}/tests/data/long_waveforms.json",
            env!("CARGO_MANIFEST_DIR")
        );
        let contents = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{path}: {e}"));
        serde_json::from_str(&contents).unwrap()
    })
}

fn cases(approximant: &str) -> Vec<&'static LongCase> {
    let cases: Vec<_> = reference()
        .cases
        .iter()
        .filter(|case| case.approximant == approximant)
        .collect();
    assert!(!cases.is_empty(), "no cases for {approximant}");
    cases
}

/// The output buffers of a grid path.
fn buffers(n: usize) -> (Vec<C64>, Vec<C64>) {
    let zero = C64::new(0.0, 0.0);
    (vec![zero; n], vec![zero; n])
}

/// Failures of the comparisons of one test.
#[derive(Default)]
struct Failures(Vec<String>);

impl Failures {
    /// Compare one path of one case with LAL and record a failure if the
    /// error exceeds `tolerance`.
    fn check(
        &mut self,
        case: &LongCase,
        path: &str,
        (hplus, hcross): (Vec<C64>, Vec<C64>),
        tolerance: f64,
    ) {
        // at least the tolerance of the suites that amplify last-bit
        // differences (the SpinTaylor angles and the (3,2) mode;
        // tests/common)
        let floor = common::tolerances(&format!("long {}", case.approximant), (0.0, 0.0)).0;
        let tolerance = tolerance.max(floor);
        let [(plus, at_plus), (cross, at_cross)] = case.errors(&hplus, &hcross);
        println!(
            "{:<58} {path:<8} plus {plus:.2e} ({at_plus:7.2} Hz), cross {cross:.2e} ({at_cross:7.2} Hz), \
             tolerance {tolerance:.1e}",
            case.label()
        );
        if !(plus <= tolerance && cross <= tolerance) {
            self.0.push(format!(
                "{} {path}: {plus:.3e}, {cross:.3e} > {tolerance:.3e}",
                case.label()
            ));
        }
    }

    fn assert_none(self) {
        assert!(
            self.0.is_empty(),
            "against lalsimulation {}:\n{}",
            reference().lalsimulation_version,
            self.0.join("\n")
        );
    }
}

/// The per-frequency path (`waveform_modes`) and the grid path of an
/// aligned-spin model built by `model` (`phi_c` is LAL's orbital `phi_ref`).
fn check_aligned<W: Waveform>(approximant: &str, model: impl Fn(&LongCase) -> W) {
    let mut failures = Failures::default();
    for case in cases(approximant) {
        let frequencies = case.sequence();
        let waveform = model(case);
        let tolerance = case.rounding_tolerance();
        let polarizations = waveform.waveform_modes(&frequencies, case.phi_ref, case.inclination);
        failures.check(case, "per-freq", polarizations, tolerance);

        let grid = FrequencyGrid::new(&frequencies);
        let (mut hplus, mut hcross) = buffers(grid.len());
        waveform.waveform_modes_on_grid(
            &grid,
            case.phi_ref,
            case.inclination,
            &mut hplus,
            &mut hcross,
        );
        failures.check(case, "grid", (hplus, hcross), tolerance);
    }
    failures.assert_none();
}

#[test]
#[cfg_attr(debug_assertions, ignore = "long grids: run with --release")]
fn taylorf2_long() {
    check_aligned("TaylorF2", |case| {
        TaylorF2::new(
            case.mass_1 + case.mass_2,
            case.mass_2 / case.mass_1,
            case.spin_1[2],
            case.spin_2[2],
            case.luminosity_distance,
        )
        .with_reference_frequency(case.f_ref)
    });
}

#[test]
#[cfg_attr(debug_assertions, ignore = "long grids: run with --release")]
fn imrphenomd_long() {
    check_aligned("IMRPhenomD", |case| {
        IMRPhenomD::new(
            case.mass_1 + case.mass_2,
            case.mass_2 / case.mass_1,
            case.spin_1[2],
            case.spin_2[2],
            case.luminosity_distance,
        )
        .with_reference_frequency(case.f_ref)
    });
}

#[test]
#[cfg_attr(debug_assertions, ignore = "long grids: run with --release")]
fn imrphenomxas_long() {
    check_aligned("IMRPhenomXAS", |case| {
        IMRPhenomXAS::from_components(
            case.mass_1,
            case.mass_2,
            case.spin_1[2],
            case.spin_2[2],
            case.luminosity_distance,
        )
        .unwrap()
        .with_reference_frequency(case.f_ref)
    });
}

#[test]
#[cfg_attr(debug_assertions, ignore = "long grids: run with --release")]
fn imrphenomxhm_long() {
    let mut failures = Failures::default();
    for case in cases("IMRPhenomXHM") {
        let frequencies = case.sequence();
        let model = IMRPhenomXHM::from_components(
            case.mass_1,
            case.mass_2,
            case.spin_1[2],
            case.spin_2[2],
            case.luminosity_distance,
        )
        .unwrap()
        .with_reference_frequency(case.f_ref);
        let polarizations = model.polarizations(&frequencies, case.phi_ref, case.inclination);
        failures.check(case, "per-freq", polarizations, case.rounding_tolerance());

        let grid = FrequencyGrid::new(&frequencies);
        let (mut hplus, mut hcross) = buffers(grid.len());
        model.polarizations_on_grid(
            &grid,
            case.phi_ref,
            case.inclination,
            &mut hplus,
            &mut hcross,
        );
        failures.check(case, "grid", (hplus, hcross), case.rounding_tolerance());
    }
    failures.assert_none();
}

/// The per-frequency (`polarizations`) and grid paths on the sequence cases
/// and `polarizations_uniform` on the uniform ones, of a precessing model
/// whose three paths are given.
fn check_precessing(
    approximant: &str,
    per_frequency: impl Fn(&LongCase, &[f64]) -> (Vec<C64>, Vec<C64>),
    on_grid: impl Fn(&LongCase, &FrequencyGrid, &mut [C64], &mut [C64]),
    uniform: impl Fn(&LongCase) -> (Vec<C64>, Vec<C64>),
) {
    let mut failures = Failures::default();
    for case in cases(approximant) {
        if case.generator == "uniform" {
            failures.check(case, "uniform", uniform(case), case.rounding_tolerance());
            continue;
        }
        let frequencies = case.sequence();
        let polarizations = per_frequency(case, &frequencies);
        failures.check(case, "per-freq", polarizations, case.rounding_tolerance());

        let grid = FrequencyGrid::new(&frequencies);
        let (mut hplus, mut hcross) = buffers(grid.len());
        on_grid(case, &grid, &mut hplus, &mut hcross);
        failures.check(case, "grid", (hplus, hcross), case.rounding_tolerance());
    }
    failures.assert_none();
}

#[test]
#[cfg_attr(debug_assertions, ignore = "long grids: run with --release")]
fn imrphenomxp_long() {
    let model =
        |case: &LongCase| IMRPhenomXP::with_version(&case.binary(), case.version()).unwrap();
    check_precessing(
        "IMRPhenomXP",
        |case, frequencies| model(case).polarizations(frequencies).unwrap(),
        |case, grid, hplus, hcross| {
            model(case)
                .polarizations_on_grid(grid, hplus, hcross)
                .unwrap()
        },
        |case| {
            model(case)
                .polarizations_uniform(case.delta_f, case.f_min, case.f_max)
                .unwrap()
        },
    );
}

#[test]
#[cfg_attr(debug_assertions, ignore = "long grids: run with --release")]
fn imrphenomxphm_long() {
    let model =
        |case: &LongCase| IMRPhenomXPHM::with_version(&case.binary(), case.version()).unwrap();
    check_precessing(
        "IMRPhenomXPHM",
        |case, frequencies| model(case).polarizations(frequencies).unwrap(),
        |case, grid, hplus, hcross| {
            model(case)
                .polarizations_on_grid(grid, hplus, hcross)
                .unwrap()
        },
        |case| {
            model(case)
                .polarizations_uniform(case.delta_f, case.f_min, case.f_max)
                .unwrap()
        },
    );
}

/// The SpinTaylor Euler angles from a backward integration of the PN
/// equations over ~900 s (~2.3 million output samples of the coarse step)
/// against `SimIMRPhenomXPSpinTaylorAngles`: bit for bit on Apple silicon.
/// Elsewhere LAL's arithmetic depends on how its build contracts
/// multiply-adds, and the adaptive step control amplifies last-bit
/// differences, so only a loose bound is checked.
#[test]
#[cfg_attr(debug_assertions, ignore = "long grids: run with --release")]
fn spintaylor_angles_long() {
    // tests/common: 0 on Apple silicon, else 1e-5 (gamma, the integral of
    // alpha, differs by 1e-6 after ~900 s on Linux x86-64)
    let tolerance = common::tolerances("spintaylor_angles_long", (0.0, 0.0)).0;
    for case in &reference().angles {
        let binary = PrecessingBinary {
            mass_1: case.mass_1,
            mass_2: case.mass_2,
            spin_1: case.spin_1,
            spin_2: case.spin_2,
            luminosity_distance: 100.0,
            inclination: 0.0,
            phi_ref: case.phi_ref,
            f_ref: case.f_ref,
        };
        let angles =
            IMRPhenomXP::spin_taylor_angles(&binary, case.f_min, case.f_max, case.delta_f).unwrap();
        assert_eq!(angles.alpha.len(), case.length, "{}", case.name);
        for (name, test, truth) in [
            ("alpha", &angles.alpha, &case.alpha),
            ("cos(beta)", &angles.cos_beta, &case.cos_beta),
            ("gamma", &angles.gamma, &case.gamma),
        ] {
            let stored: Vec<f64> = test.iter().step_by(case.stride).copied().collect();
            assert_eq!(stored.len(), truth.len(), "{}", case.name);
            let max_diff = stored
                .iter()
                .zip(truth)
                .map(|(a, b)| (a - b).abs())
                .fold(0.0, f64::max);
            println!(
                "{} {name}: max |difference| {max_diff:.3e} over {} values",
                case.name,
                truth.len()
            );
            assert!(
                max_diff <= tolerance,
                "{} {name}: {max_diff:e} (lalsimulation {})",
                case.name,
                reference().lalsimulation_version
            );
        }
    }
}

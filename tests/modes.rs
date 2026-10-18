//! Spherical-harmonic modes in the frequency domain: `choose_fd_modes` and
//! `FrequencyDomainModel::modes` against LALSimulation's
//! `SimInspiralChooseFDModes` (reference data from
//! `scripts/generate_reference_data.py --suites choose_fd_modes`), and the
//! polarizations built from the modes against the models' polarizations.
mod common;

use num_complex::Complex;
use serde::Deserialize;
use waveformrs::{
    choose_fd_model, choose_fd_modes, imrphenomx::PrecessionVersion, Approximant, Binary, Error,
    FrequencyDomainModel, FrequencyModes, WaveformOptions,
};

type C64 = Complex<f64>;

/// `LAL_MTSUN_SI`
const MTSUN_SI: f64 = 4.925_490_947_641_267e-6;

#[derive(Deserialize)]
struct ModesFile {
    lalsimulation_version: String,
    cases: Vec<ModesCase>,
}

#[derive(Deserialize)]
struct ModesCase {
    approximant: String,
    mass_1: f64,
    mass_2: f64,
    spin_1: [f64; 3],
    spin_2: [f64; 3],
    luminosity_distance: f64,
    inclination: f64,
    phi_ref: f64,
    f_min: f64,
    f_ref: f64,
    delta_f: f64,
    f_max: f64,
    #[serde(default)]
    lal_flags: std::collections::HashMap<String, f64>,
    #[serde(default)]
    mode_array: Vec<(u32, i32)>,
    fd_modes_stride: usize,
    length: usize,
    modes: Vec<LalMode>,
}

#[derive(Deserialize)]
struct LalMode {
    ell: u32,
    emm: i32,
    positive_real: Vec<f64>,
    positive_imag: Vec<f64>,
    negative_real: Vec<f64>,
    negative_imag: Vec<f64>,
}

impl LalMode {
    /// The stored samples of `h_lm(k delta_f)` followed by those of
    /// `h_lm(-k delta_f)`.
    fn samples(&self) -> Vec<C64> {
        let zip = |re: &[f64], im: &[f64]| -> Vec<C64> {
            re.iter().zip(im).map(|(&r, &i)| C64::new(r, i)).collect()
        };
        let mut samples = zip(&self.positive_real, &self.positive_imag);
        samples.extend(zip(&self.negative_real, &self.negative_imag));
        samples
    }
}

impl ModesCase {
    fn binary(&self) -> Binary {
        Binary {
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

    fn approximant(&self) -> Approximant {
        self.approximant.parse().unwrap()
    }

    fn version(&self) -> PrecessionVersion {
        match self.lal_flags.get("PhenomXPrecVersion").copied() {
            Some(102.0) => PrecessionVersion::Nnlo,
            Some(223.0) => PrecessionVersion::Msa,
            _ => PrecessionVersion::SpinTaylor,
        }
    }

    fn options(&self) -> WaveformOptions {
        let mut options = WaveformOptions::default().with_precession_version(self.version());
        if !self.mode_array.is_empty() {
            options = options.with_mode_array(&self.mode_array);
        }
        options
    }

    fn label(&self) -> String {
        format!(
            "{} m1={} m2={} s1={:?} s2={:?} version={:?} f_ref={}",
            self.approximant,
            self.mass_1,
            self.mass_2,
            self.spin_1,
            self.spin_2,
            self.version(),
            self.f_ref
        )
    }

    /// The stored samples of a mode of `modes` (indices `k delta_f` of the
    /// full grid starting at `offset`).
    fn samples_of(&self, modes: &FrequencyModes, ell: u32, emm: i32, offset: usize) -> Vec<C64> {
        let mode = modes
            .mode(ell, emm)
            .unwrap_or_else(|| panic!("({ell},{emm}) missing: {}", self.label()));
        let zero = C64::new(0.0, 0.0);
        let at = |h: &[C64], k: usize| {
            k.checked_sub(offset)
                .and_then(|i| h.get(i))
                .copied()
                .unwrap_or(zero)
        };
        let stride = self.fd_modes_stride;
        let mut samples: Vec<C64> = (0..self.length)
            .step_by(stride)
            .map(|k| at(&mode.positive, k))
            .collect();
        samples.extend(
            (0..self.length)
                .step_by(stride)
                .map(|k| at(&mode.negative, k)),
        );
        samples
    }
}

fn load() -> ModesFile {
    let path = format!(
        "{}/tests/data/choose_fd_modes.json",
        env!("CARGO_MANIFEST_DIR")
    );
    serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap()
}

/// Compare every mode of `modes` with LAL's, returning the failures.
fn compare_modes(
    case: &ModesCase,
    modes: &FrequencyModes,
    offset: usize,
    what: &str,
    failures: &mut Vec<String>,
) {
    let mut lal_modes: Vec<(u32, i32)> = case.modes.iter().map(|m| (m.ell, m.emm)).collect();
    let mut ours: Vec<(u32, i32)> = modes.modes.iter().map(|m| (m.ell, m.m)).collect();
    lal_modes.sort_unstable();
    ours.sort_unstable();
    assert_eq!(ours, lal_modes, "{what}: mode list of {}", case.label());
    for lal in &case.modes {
        let truth = lal.samples();
        let test = case.samples_of(modes, lal.ell, lal.emm, offset);
        // the stored samples are equally spaced on each half
        let grid: Vec<f64> = (0..truth.len()).map(|i| i as f64).collect();
        let result = common::compare(&grid, &test, &truth);
        println!(
            "{what} ({},{}): max rel err {:.3e}, mismatch {:.3e} | {}",
            lal.ell,
            lal.emm,
            result.max_relative_error,
            result.mismatch,
            case.label()
        );
        if !(result.max_relative_error <= 1e-9 && result.mismatch <= 1e-12) {
            failures.push(format!(
                "{what} ({},{}) {result:?} {}",
                lal.ell,
                lal.emm,
                case.label()
            ));
        }
    }
}

/// `choose_fd_modes` against `SimInspiralChooseFDModes` (uniform grids).
#[test]
fn choose_fd_modes_matches_lalsimulation() {
    let reference = load();
    let mut failures = Vec::new();
    for case in &reference.cases {
        let modes = choose_fd_modes(
            case.approximant(),
            &case.binary(),
            case.delta_f,
            case.f_min,
            case.f_max,
            &case.options(),
        )
        .unwrap();
        assert_eq!(modes.frequencies.len(), case.length, "{}", case.label());
        assert_eq!(modes.frequencies[1], case.delta_f);
        compare_modes(case, &modes, 0, "choose_fd_modes", &mut failures);
    }
    assert!(
        failures.is_empty(),
        "modes differ from lalsimulation {}:\n{}",
        reference.lalsimulation_version,
        failures.join("\n")
    );
}

/// `FrequencyDomainModel::modes` on the frequencies `k delta_f` from `f_min`
/// to `min(f_max, f_cut)` gives the modes of `SimInspiralChooseFDModes`
/// where the frequency sequence and the uniform grid agree in LAL:
/// IMRPhenomXHM, and IMRPhenomXPHM with the MSA and NNLO angles (the
/// SpinTaylor angles of the uniform drivers start the PN integration
/// `3 delta_f` instead of 0.5 Hz below `f_min`).
#[test]
fn modes_on_frequency_sequences_match_lalsimulation() {
    let reference = load();
    let mut failures = Vec::new();
    for case in &reference.cases {
        if case.approximant() == Approximant::IMRPhenomXPHM
            && case.version() == PrecessionVersion::SpinTaylor
        {
            continue;
        }
        let total_mass = case.mass_1 + case.mass_2;
        let f_cut = 0.3 / (total_mass * MTSUN_SI);
        let f_max_prime = if case.f_max == 0.0 || case.f_max > f_cut {
            f_cut
        } else {
            case.f_max
        };
        let start = (case.f_min / case.delta_f) as usize;
        let stop = (f_max_prime / case.delta_f) as usize + 1;
        let band: Vec<f64> = (start..stop).map(|k| k as f64 * case.delta_f).collect();
        let mut binary = case.binary();
        if binary.f_ref == 0.0 {
            // the uniform drivers reference to f_min
            binary.f_ref = case.f_min;
        }
        let model = choose_fd_model(case.approximant(), &binary, &case.options()).unwrap();
        let modes = model.modes(&band).unwrap();
        compare_modes(case, &modes, start, "modes", &mut failures);
    }
    assert!(
        failures.is_empty(),
        "modes differ from lalsimulation {}:\n{}",
        reference.lalsimulation_version,
        failures.join("\n")
    );
}

fn aligned_binary() -> Binary {
    Binary {
        mass_1: 36.0,
        mass_2: 29.0,
        spin_1: [0.0, 0.0, 0.3],
        spin_2: [0.0, 0.0, -0.4],
        luminosity_distance: 500.0,
        inclination: 0.7,
        phi_ref: 1.2,
        f_ref: 25.0,
    }
}

fn precessing_binary() -> Binary {
    Binary {
        spin_1: [0.4, 0.2, 0.3],
        spin_2: [-0.2, 0.3, -0.1],
        ..aligned_binary()
    }
}

fn max_relative_difference(a: &[C64], b: &[C64]) -> f64 {
    let scale = b.iter().map(|x| x.norm()).fold(0.0, f64::max);
    a.iter()
        .zip(b)
        .map(|(x, y)| (x - y).norm())
        .fold(0.0, f64::max)
        / scale
}

/// The polarizations built from the modes (`FrequencyModes::polarizations`)
/// agree with the models' polarizations, below the cutoff `Mf = 0.3`.
#[test]
fn modes_reproduce_the_polarizations() {
    let frequencies: Vec<f64> = (0..400).map(|i| 20.0 * 1.008f64.powi(i)).collect();
    let f_cut = 0.3 / (65.0 * MTSUN_SI);
    assert!(frequencies[frequencies.len() - 1] < f_cut);
    let cases = [
        (Approximant::TaylorF2, PrecessionVersion::SpinTaylor, 1e-14),
        (
            Approximant::IMRPhenomD,
            PrecessionVersion::SpinTaylor,
            1e-14,
        ),
        (
            Approximant::IMRPhenomXAS,
            PrecessionVersion::SpinTaylor,
            1e-14,
        ),
        (
            Approximant::IMRPhenomXHM,
            PrecessionVersion::SpinTaylor,
            1e-13,
        ),
        (Approximant::IMRPhenomXP, PrecessionVersion::Nnlo, 1e-13),
        (Approximant::IMRPhenomXP, PrecessionVersion::Msa, 1e-13),
        (
            Approximant::IMRPhenomXP,
            PrecessionVersion::SpinTaylor,
            1e-13,
        ),
        (Approximant::IMRPhenomXPHM, PrecessionVersion::Nnlo, 1e-13),
        (Approximant::IMRPhenomXPHM, PrecessionVersion::Msa, 1e-13),
        // each J-frame mode has its own range of SpinTaylor angles, as in LAL
        (
            Approximant::IMRPhenomXPHM,
            PrecessionVersion::SpinTaylor,
            1e-5,
        ),
    ];
    for f_ref in [25.0, 0.0] {
        for (approximant, version, tolerance) in cases {
            let binary = if approximant.is_precessing() {
                precessing_binary()
            } else {
                aligned_binary()
            };
            let binary = Binary { f_ref, ..binary };
            let options = WaveformOptions::default().with_precession_version(version);
            let model = choose_fd_model(approximant, &binary, &options).unwrap();
            let (hp, hc) = model.polarizations(&frequencies).unwrap();
            let modes = model.modes(&frequencies).unwrap();
            assert_eq!(modes.frequencies, frequencies);
            let (hp_modes, hc_modes) = modes.polarizations();
            let error = max_relative_difference(&hp_modes, &hp)
                .max(max_relative_difference(&hc_modes, &hc));
            println!("{approximant} {version:?} f_ref={f_ref}: {error:.3e}");
            assert!(
                error < tolerance,
                "{approximant} {version:?} f_ref={f_ref}: {error:e}"
            );
        }
    }
}

/// The mode content of each model, and the halves on which aligned-spin
/// modes live.
#[test]
fn mode_content() {
    let frequencies: Vec<f64> = (0..50).map(|i| 20.0 + 4.0 * i as f64).collect();
    let expected = [
        (Approximant::TaylorF2, 2),
        (Approximant::IMRPhenomD, 2),
        (Approximant::IMRPhenomXAS, 2),
        (Approximant::IMRPhenomXHM, 10),
        (Approximant::IMRPhenomXP, 5),
        (Approximant::IMRPhenomXPHM, 21),
    ];
    for (approximant, count) in expected {
        let binary = if approximant.is_precessing() {
            precessing_binary()
        } else {
            aligned_binary()
        };
        let modes = choose_fd_model(approximant, &binary, &WaveformOptions::default())
            .unwrap()
            .modes(&frequencies)
            .unwrap();
        assert_eq!(modes.modes.len(), count, "{approximant}");
        for mode in &modes.modes {
            assert_eq!(mode.positive.len(), frequencies.len());
            assert_eq!(mode.negative.len(), frequencies.len());
            if !approximant.is_precessing() {
                // m < 0 at positive frequencies, m > 0 at negative ones
                let (live, zero) = if mode.m < 0 {
                    (&mode.positive, &mode.negative)
                } else {
                    (&mode.negative, &mode.positive)
                };
                assert!(live.iter().all(|h| h.norm() > 0.0), "{approximant}");
                assert!(zero.iter().all(|h| h.norm() == 0.0), "{approximant}");
            } else {
                assert!(mode.positive.iter().any(|h| h.norm() > 0.0));
                assert!(mode.negative.iter().any(|h| h.norm() > 0.0));
            }
        }
    }
    // the (2,2) modes of IMRPhenomXAS are those of IMRPhenomXHM
    let b = aligned_binary();
    let xas = choose_fd_model(Approximant::IMRPhenomXAS, &b, &WaveformOptions::default())
        .unwrap()
        .modes(&frequencies)
        .unwrap();
    let xhm = choose_fd_model(Approximant::IMRPhenomXHM, &b, &WaveformOptions::default())
        .unwrap()
        .modes(&frequencies)
        .unwrap();
    for m in [-2, 2] {
        let (a, b) = (xas.mode(2, m).unwrap(), xhm.mode(2, m).unwrap());
        for (x, y) in a
            .positive
            .iter()
            .zip(&b.positive)
            .chain(a.negative.iter().zip(&b.negative))
        {
            assert!((x - y).norm() <= 1e-12 * y.norm(), "{x} != {y}");
        }
    }
}

#[test]
fn choose_fd_modes_supports_what_lal_supports() {
    for approximant in [
        Approximant::TaylorF2,
        Approximant::IMRPhenomD,
        Approximant::IMRPhenomXAS,
        Approximant::IMRPhenomXP,
        Approximant::NRSur7dq4,
    ] {
        let binary = if approximant.is_precessing() {
            precessing_binary()
        } else {
            aligned_binary()
        };
        assert!(matches!(
            choose_fd_modes(
                approximant,
                &binary,
                1.0,
                20.0,
                512.0,
                &WaveformOptions::default()
            ),
            Err(Error::Unsupported { .. })
        ));
    }
    assert!(choose_fd_modes(
        Approximant::IMRPhenomXHM,
        &aligned_binary(),
        0.0,
        20.0,
        512.0,
        &WaveformOptions::default()
    )
    .is_err());
}

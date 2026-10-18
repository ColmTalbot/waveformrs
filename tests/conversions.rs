//! The time/frequency-domain conversions against LALSimulation:
//! `SimInspiralTD` for frequency-domain approximants, `SimInspiralFD` for
//! NRSur7dq4 (and for frequency-domain approximants), `SimInspiralTD` for
//! NRSur7dq4.
//!
//! Reference data from `scripts/generate_conversion_reference.py`. Epochs,
//! lengths and sampling intervals must agree exactly; the values (every
//! `stride`-th sample stored) to `1e-10` of the largest absolute value of
//! the series. The NRSur7dq4 suites need `NRSur7dq4_v1.0.h5` (see
//! `tests/nrsur.rs`) and are skipped without it.
#![cfg(feature = "fft")]

mod common;

use std::collections::HashMap;
use std::path::PathBuf;

use num_complex::Complex;
use serde::de::DeserializeOwned;
use serde::Deserialize;
use waveformrs::conversions::{
    nrsur7dq4_fd, nrsur7dq4_td, simulate_fd, simulate_td, FdWaveform, FrequencySeries, GpsTime,
    TimeSeries,
};
use waveformrs::imrphenomx::PrecessionVersion;
use waveformrs::nrsur::{NRSur7dq4, NRSur7dq4Data};
use waveformrs::{Approximant, Binary, WaveformOptions};

/// Tolerance on `max |rust - lal| / max |lal|`.
const TOLERANCE: f64 = 1e-10;

/// [`TOLERANCE`], or the cross-build tolerance of the higher-mode models on
/// builds the crate does not mirror (tests/common).
fn tolerance(approximant: &str) -> f64 {
    common::tolerances(&format!("conversions {approximant}"), (TOLERANCE, 0.0)).0
}

#[derive(Deserialize)]
struct File<C> {
    lalsimulation_version: String,
    cases: Vec<C>,
}

fn load<C: DeserializeOwned>(name: &str) -> File<C> {
    let path = format!("{}/tests/data/{name}.json", env!("CARGO_MANIFEST_DIR"));
    let contents = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{path}: {e}"));
    serde_json::from_str(&contents).unwrap()
}

#[derive(Deserialize, Clone, Debug)]
struct Parameters {
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
    #[serde(default)]
    lal_flags: HashMap<String, f64>,
}

impl Parameters {
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

    /// The `LALDict` of the case: the Euler angles (multibanding is always
    /// off in the reference data).
    fn options(&self) -> WaveformOptions {
        let version = match self.lal_flags.get("PhenomXPrecVersion").map(|&v| v as u32) {
            Some(102) => PrecessionVersion::Nnlo,
            Some(223) => PrecessionVersion::Msa,
            Some(320) | None => PrecessionVersion::SpinTaylor,
            Some(other) => panic!("PhenomXPrecVersion {other}"),
        };
        WaveformOptions::default().with_precession_version(version)
    }

    fn label(&self) -> String {
        format!(
            "{} m1={} m2={} s1={:?} s2={:?} f_min={} f_ref={}",
            self.approximant,
            self.mass_1,
            self.mass_2,
            self.spin_1,
            self.spin_2,
            self.f_min,
            self.f_ref
        )
    }
}

#[derive(Deserialize, Clone, Debug)]
struct TdCase {
    #[serde(flatten)]
    binary: Parameters,
    delta_t: f64,
    epoch_seconds: i32,
    epoch_nanoseconds: i64,
    delta_t_out: f64,
    length: usize,
    stride: usize,
    hplus: Vec<f64>,
    hcross: Vec<f64>,
    max_abs_hplus: f64,
    max_abs_hcross: f64,
}

#[derive(Deserialize, Clone, Debug)]
struct FdCase {
    #[serde(flatten)]
    binary: Parameters,
    delta_f: f64,
    f_max: f64,
    epoch_seconds: i32,
    epoch_nanoseconds: i64,
    delta_f_out: f64,
    length: usize,
    stride: usize,
    hplus_real: Vec<f64>,
    hplus_imag: Vec<f64>,
    hcross_real: Vec<f64>,
    hcross_imag: Vec<f64>,
    max_abs_hplus: f64,
    max_abs_hcross: f64,
}

/// `max |test[j stride] - reference[j]| / max_abs`.
fn relative_error<T: Copy, F: Fn(T, T) -> f64>(
    test: &[T],
    reference: &[T],
    stride: usize,
    max_abs: f64,
    distance: F,
) -> f64 {
    let worst = test
        .iter()
        .step_by(stride)
        .zip(reference)
        .map(|(&a, &b)| distance(a, b))
        .fold(0.0, f64::max);
    worst / max_abs
}

/// Compare a time-domain result with a reference case; the description of
/// any failure.
fn check_td(case: &TdCase, hplus: &TimeSeries, hcross: &TimeSeries) -> Vec<String> {
    let mut failures = Vec::new();
    let epoch = GpsTime::from_parts(case.epoch_seconds, case.epoch_nanoseconds);
    for (name, series, reference, max_abs) in [
        ("plus", hplus, &case.hplus, case.max_abs_hplus),
        ("cross", hcross, &case.hcross, case.max_abs_hcross),
    ] {
        if series.epoch != epoch {
            failures.push(format!(
                "{name}: epoch {:?} != {epoch:?} | {}",
                series.epoch,
                case.binary.label()
            ));
        }
        if series.data.len() != case.length || series.delta_t != case.delta_t_out {
            failures.push(format!(
                "{name}: length {} delta_t {} != {} {} | {}",
                series.data.len(),
                series.delta_t,
                case.length,
                case.delta_t_out,
                case.binary.label()
            ));
            continue;
        }
        let error = relative_error(&series.data, reference, case.stride, max_abs, |a, b| {
            (a - b).abs()
        });
        println!(
            "{name:>5}: max rel err {error:.2e} | n = {} delta_t = {} epoch = {} | {}",
            case.length,
            case.delta_t,
            epoch.seconds(),
            case.binary.label()
        );
        if error.is_nan() || error > tolerance(&case.binary.approximant) {
            failures.push(format!("{name}: error {error:e} | {}", case.binary.label()));
        }
    }
    failures
}

fn check_fd(case: &FdCase, hplus: &FrequencySeries, hcross: &FrequencySeries) -> Vec<String> {
    let mut failures = Vec::new();
    let epoch = GpsTime::from_parts(case.epoch_seconds, case.epoch_nanoseconds);
    let zip = |re: &[f64], im: &[f64]| -> Vec<Complex<f64>> {
        re.iter()
            .zip(im)
            .map(|(&r, &i)| Complex::new(r, i))
            .collect()
    };
    for (name, series, reference, max_abs) in [
        (
            "plus",
            hplus,
            zip(&case.hplus_real, &case.hplus_imag),
            case.max_abs_hplus,
        ),
        (
            "cross",
            hcross,
            zip(&case.hcross_real, &case.hcross_imag),
            case.max_abs_hcross,
        ),
    ] {
        if series.epoch != epoch {
            failures.push(format!(
                "{name}: epoch {:?} != {epoch:?} | {}",
                series.epoch,
                case.binary.label()
            ));
        }
        if series.data.len() != case.length || series.delta_f != case.delta_f_out {
            failures.push(format!(
                "{name}: length {} delta_f {} != {} {} | {}",
                series.data.len(),
                series.delta_f,
                case.length,
                case.delta_f_out,
                case.binary.label()
            ));
            continue;
        }
        let error = relative_error(&series.data, &reference, case.stride, max_abs, |a, b| {
            (a - b).norm()
        });
        println!(
            "{name:>5}: max rel err {error:.2e} | n = {} delta_f = {} f_max = {} epoch = {} | {}",
            case.length,
            case.delta_f,
            case.f_max,
            epoch.seconds(),
            case.binary.label()
        );
        if error.is_nan() || error > tolerance(&case.binary.approximant) {
            failures.push(format!("{name}: error {error:e} | {}", case.binary.label()));
        }
    }
    failures
}

fn report(name: &str, version: &str, failures: Vec<String>) {
    assert!(
        failures.is_empty(),
        "{name}: {} comparisons failed against lalsimulation {version}:\n{}",
        failures.len(),
        failures.join("\n")
    );
}

/// `SimInspiralTD` for a frequency-domain approximant.
fn check_td_suite(name: &str) {
    let reference: File<TdCase> = load(name);
    let mut failures = Vec::new();
    for case in &reference.cases {
        let b = &case.binary;
        match simulate_td(
            b.approximant(),
            &b.binary(),
            &b.options(),
            case.delta_t,
            b.f_min,
        ) {
            Ok(h) => failures.extend(check_td(case, &h.hplus, &h.hcross)),
            Err(err) => failures.push(format!("error {err} | {}", case.binary.label())),
        }
    }
    report(name, &reference.lalsimulation_version, failures);
}

#[test]
fn imrphenomd_td_matches_lalsimulation() {
    check_td_suite("conversions_td_imrphenomd");
}

#[test]
fn imrphenomxas_td_matches_lalsimulation() {
    check_td_suite("conversions_td_imrphenomxas");
}

#[test]
fn imrphenomxhm_td_matches_lalsimulation() {
    check_td_suite("conversions_td_imrphenomxhm");
}

#[test]
fn imrphenomxp_td_matches_lalsimulation() {
    check_td_suite("conversions_td_imrphenomxp");
}

#[test]
fn imrphenomxphm_td_matches_lalsimulation() {
    check_td_suite("conversions_td_imrphenomxphm");
}

#[test]
fn taylorf2_td_matches_lalsimulation() {
    check_td_suite("conversions_td_taylorf2");
}

/// `SimInspiralFD` for frequency-domain approximants (conditioning only).
#[test]
fn conditioned_fd_matches_lalsimulation() {
    let name = "conversions_fd";
    let reference: File<FdCase> = load(name);
    let mut failures = Vec::new();
    for case in &reference.cases {
        let b = &case.binary;
        let model = FdWaveform::new(b.approximant(), b.binary(), b.options()).unwrap();
        match model.conditioned(case.delta_f, b.f_min, case.f_max) {
            Ok(h) => failures.extend(check_fd(case, &h.hplus, &h.hcross)),
            Err(err) => failures.push(format!("error {err} | {}", case.binary.label())),
        }
    }
    report(name, &reference.lalsimulation_version, failures);
}

fn nrsur_data() -> Option<&'static NRSur7dq4Data> {
    if waveformrs::nrsur::find_data_file().is_none() {
        let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("lal_data")
            .join(waveformrs::nrsur::NRSUR7DQ4_DATAFILE);
        if path.is_file() {
            std::env::set_var("NRSUR7DQ4_DATA", path);
        } else {
            println!(
                "skipping: {} not found (set NRSUR7DQ4_DATA or LAL_DATA_PATH)",
                waveformrs::nrsur::NRSUR7DQ4_DATAFILE
            );
            return None;
        }
    }
    Some(NRSur7dq4Data::global().unwrap())
}

/// `SimInspiralFD` for NRSur7dq4.
#[test]
fn nrsur7dq4_fd_matches_lalsimulation() {
    let Some(data) = nrsur_data() else {
        return;
    };
    let name = "conversions_fd_nrsur7dq4";
    let reference: File<FdCase> = load(name);
    let mut failures = Vec::new();
    for case in &reference.cases {
        let model = NRSur7dq4::new(data, &case.binary.binary()).unwrap();
        match nrsur7dq4_fd(&model, case.delta_f, case.binary.f_min, case.f_max) {
            Ok(h) => {
                failures.extend(check_fd(case, &h.hplus, &h.hcross));
                // the same through the approximant (with the global data)
                let b = &case.binary;
                let again = simulate_fd(
                    b.approximant(),
                    &b.binary(),
                    &b.options(),
                    case.delta_f,
                    b.f_min,
                    case.f_max,
                )
                .unwrap();
                assert_eq!(again, h, "{}", b.label());
            }
            Err(err) => failures.push(format!("error {err} | {}", case.binary.label())),
        }
    }
    report(name, &reference.lalsimulation_version, failures);
}

/// `SimInspiralTD` for NRSur7dq4: the tapered waveform, bit for bit.
#[test]
fn nrsur7dq4_td_matches_lalsimulation() {
    let Some(data) = nrsur_data() else {
        return;
    };
    let name = "conversions_td_nrsur7dq4";
    let reference: File<TdCase> = load(name);
    let mut failures = Vec::new();
    for case in &reference.cases {
        let model = NRSur7dq4::new(data, &case.binary.binary()).unwrap();
        match nrsur7dq4_td(&model, case.delta_t, case.binary.f_min) {
            Ok(h) => {
                failures.extend(check_td(case, &h.hplus, &h.hcross));
                if cfg!(target_arch = "aarch64") && cfg!(target_vendor = "apple") {
                    for (series, reference) in [(&h.hplus, &case.hplus), (&h.hcross, &case.hcross)]
                    {
                        let same = series
                            .data
                            .iter()
                            .step_by(case.stride)
                            .zip(reference)
                            .all(|(a, b)| a.to_bits() == b.to_bits());
                        if !same {
                            failures.push(format!("not bit for bit | {}", case.binary.label()));
                        }
                    }
                }
            }
            Err(err) => failures.push(format!("error {err} | {}", case.binary.label())),
        }
    }
    report(name, &reference.lalsimulation_version, failures);
}

//! Evaluate waveforms for a batch of parameter sets read as JSON from stdin
//! and write the polarizations as JSON to stdout.
//!
//! Used by `scripts/validate_against_lalsimulation.py` to compare against
//! lalsimulation over many random draws.
//!
//! Input:  {"cases": [{"approximant": ..., "mass_1": ..., "mass_2": ...,
//!          "spin_1": [x, y, z], "spin_2": [x, y, z], "luminosity_distance": ...,
//!          "inclination": ..., "phi_ref": ..., "f_ref": ..., "f_min": ...,
//!          "frequencies": [...]}, ...]}
//! Output: {"results": [{"hplus_real": [...], "hplus_imag": [...],
//!          "hcross_real": [...], "hcross_imag": [...], "seconds": ...}, ...]}
use std::io::Read;
use std::time::Instant;

use num_complex::Complex;
use serde::{Deserialize, Serialize};

use waveformrs::{
    imrphenomd::IMRPhenomD, imrphenomx::IMRPhenomXAS, taylor::TaylorF2, waveform::Waveform,
};

#[derive(Deserialize)]
struct Input {
    cases: Vec<Case>,
}

#[derive(Deserialize)]
#[allow(dead_code)]
struct Case {
    approximant: String,
    mass_1: f64,
    mass_2: f64,
    spin_1: [f64; 3],
    spin_2: [f64; 3],
    luminosity_distance: f64,
    inclination: f64,
    phi_ref: f64,
    f_ref: f64,
    f_min: f64,
    frequencies: Vec<f64>,
}

#[derive(Serialize)]
struct Output {
    hplus_real: Vec<f64>,
    hplus_imag: Vec<f64>,
    hcross_real: Vec<f64>,
    hcross_imag: Vec<f64>,
    seconds: f64,
}

type Polarizations = (Vec<Complex<f64>>, Vec<Complex<f64>>);

fn generate(case: &Case) -> Polarizations {
    let total_mass = case.mass_1 + case.mass_2;
    let mass_ratio = case.mass_2 / case.mass_1;
    match case.approximant.as_str() {
        "TaylorF2" => TaylorF2::new(
            total_mass,
            mass_ratio,
            case.spin_1[2],
            case.spin_2[2],
            case.luminosity_distance,
        )
        .with_reference_frequency(case.f_ref)
        .waveform_modes(&case.frequencies, case.phi_ref, case.inclination),
        "IMRPhenomD" => IMRPhenomD::new(
            total_mass,
            mass_ratio,
            case.spin_1[2],
            case.spin_2[2],
            case.luminosity_distance,
        )
        .with_reference_frequency(case.f_ref)
        .waveform_modes(&case.frequencies, 2.0 * case.phi_ref, case.inclination),
        "IMRPhenomXAS" => IMRPhenomXAS::new(
            total_mass,
            mass_ratio,
            case.spin_1[2],
            case.spin_2[2],
            case.luminosity_distance,
        )
        .with_reference_frequency(case.f_ref)
        .waveform_modes(&case.frequencies, 2.0 * case.phi_ref, case.inclination),
        other => panic!("unsupported approximant {other}"),
    }
}

fn main() {
    let mut buffer = String::new();
    std::io::stdin().read_to_string(&mut buffer).unwrap();
    let input: Input = serde_json::from_str(&buffer).unwrap();
    let results: Vec<Output> = input
        .cases
        .iter()
        .map(|case| {
            // time the generation (median of a few repeats to reduce noise)
            let mut timings = Vec::new();
            let mut output = generate(case);
            for _ in 0..5 {
                let start = Instant::now();
                output = generate(case);
                timings.push(start.elapsed().as_secs_f64());
            }
            timings.sort_by(f64::total_cmp);
            let (hp, hc) = output;
            Output {
                hplus_real: hp.iter().map(|h| h.re).collect(),
                hplus_imag: hp.iter().map(|h| h.im).collect(),
                hcross_real: hc.iter().map(|h| h.re).collect(),
                hcross_imag: hc.iter().map(|h| h.im).collect(),
                seconds: timings[timings.len() / 2],
            }
        })
        .collect();
    println!(
        "{}",
        serde_json::to_string(&serde_json::json!({ "results": results })).unwrap()
    );
}

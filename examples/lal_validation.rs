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
//!          "hcross_real": [...], "hcross_imag": [...], "seconds": ...,
//!          "error": null | "..."}, ...]}
//!
//! Precessing approximants are `IMRPhenomXP` (SpinTaylor angles,
//! `PhenomXPrecVersion` 320), `IMRPhenomXP_NNLO` (102) and
//! `IMRPhenomXP_MSA` (223).
use std::io::Read;
use std::time::Instant;

use num_complex::Complex;
use serde::{Deserialize, Serialize};

use waveformrs::{
    imrphenomd::IMRPhenomD,
    imrphenomx::{IMRPhenomXAS, IMRPhenomXHM, IMRPhenomXP, PrecessingBinary, PrecessionVersion},
    taylor::TaylorF2,
    waveform::Waveform,
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
    error: Option<String>,
}

type Polarizations = (Vec<Complex<f64>>, Vec<Complex<f64>>);

fn precessing(case: &Case, version: PrecessionVersion) -> Result<Polarizations, String> {
    let binary = PrecessingBinary {
        mass_1: case.mass_1,
        mass_2: case.mass_2,
        spin_1: case.spin_1,
        spin_2: case.spin_2,
        luminosity_distance: case.luminosity_distance,
        inclination: case.inclination,
        phi_ref: case.phi_ref,
        f_ref: case.f_ref,
    };
    IMRPhenomXP::with_version(&binary, version)
        .and_then(|model| model.polarizations(&case.frequencies))
        .map_err(|err| err.to_string())
}

fn generate(case: &Case) -> Result<Polarizations, String> {
    let total_mass = case.mass_1 + case.mass_2;
    let mass_ratio = case.mass_2 / case.mass_1;
    Ok(match case.approximant.as_str() {
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
        "IMRPhenomXHM" => IMRPhenomXHM::from_components(
            case.mass_1,
            case.mass_2,
            case.spin_1[2],
            case.spin_2[2],
            case.luminosity_distance,
        )
        .map_err(|err| err.to_string())?
        .with_reference_frequency(case.f_ref)
        .polarizations(&case.frequencies, case.phi_ref, case.inclination),
        "IMRPhenomXP" => precessing(case, PrecessionVersion::SpinTaylor)?,
        "IMRPhenomXP_NNLO" => precessing(case, PrecessionVersion::Nnlo)?,
        "IMRPhenomXP_MSA" => precessing(case, PrecessionVersion::Msa)?,
        other => panic!("unsupported approximant {other}"),
    })
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
            let (error, (hp, hc)) = match output {
                Ok(polarizations) => (None, polarizations),
                Err(error) => (Some(error), (Vec::new(), Vec::new())),
            };
            Output {
                hplus_real: hp.iter().map(|h| h.re).collect(),
                hplus_imag: hp.iter().map(|h| h.im).collect(),
                hcross_real: hc.iter().map(|h| h.re).collect(),
                hcross_imag: hc.iter().map(|h| h.im).collect(),
                seconds: timings[timings.len() / 2],
                error,
            }
        })
        .collect();
    println!(
        "{}",
        serde_json::to_string(&serde_json::json!({ "results": results })).unwrap()
    );
}

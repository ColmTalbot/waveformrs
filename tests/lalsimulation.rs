//! Comparisons against LALSimulation.
//!
//! Reference data are produced by `scripts/generate_reference_data.py` using
//! `SimInspiralChooseFDWaveformSequence` and stored in `tests/data`.
mod common;

use serde::Deserialize;
use waveformrs::{
    imrphenomd::IMRPhenomD,
    imrphenomx::{IMRPhenomXAS, IMRPhenomXHM, IMRPhenomXP, PrecessingBinary, PrecessionVersion},
    taylor::TaylorF2,
    waveform::Waveform,
};

#[test]
fn taylorf2_matches_lalsimulation() {
    common::check_all("taylorf2", 1e-10, 1e-14, |case| {
        TaylorF2::new(
            case.total_mass(),
            case.mass_ratio(),
            case.spin_1[2],
            case.spin_2[2],
            case.luminosity_distance,
        )
        .with_reference_frequency(case.f_ref)
        .waveform_modes(&case.frequencies, case.phi_ref, case.inclination)
    });
}

#[test]
fn imrphenomd_matches_lalsimulation() {
    common::check_all("imrphenomd", 1e-10, 1e-14, |case| {
        IMRPhenomD::new(
            case.total_mass(),
            case.mass_ratio(),
            case.spin_1[2],
            case.spin_2[2],
            case.luminosity_distance,
        )
        .with_reference_frequency(case.f_ref)
        // IMRPhenomD's phi_c is the GW phase, twice the orbital phase.
        .waveform_modes(&case.frequencies, 2.0 * case.phi_ref, case.inclination)
    });
}

#[test]
fn imrphenomxas_matches_lalsimulation() {
    common::check_all("imrphenomxas", 1e-9, 1e-12, |case| {
        IMRPhenomXAS::new(
            case.total_mass(),
            case.mass_ratio(),
            case.spin_1[2],
            case.spin_2[2],
            case.luminosity_distance,
        )
        .with_reference_frequency(case.f_ref)
        // As for IMRPhenomD, phi_c is twice the orbital reference phase.
        .waveform_modes(&case.frequencies, 2.0 * case.phi_ref, case.inclination)
    });
}

fn precessing_binary(case: &common::Case) -> PrecessingBinary {
    PrecessingBinary {
        mass_1: case.mass_1,
        mass_2: case.mass_2,
        spin_1: case.spin_1,
        spin_2: case.spin_2,
        luminosity_distance: case.luminosity_distance,
        inclination: case.inclination,
        phi_ref: case.phi_ref,
        f_ref: case.f_ref,
    }
}

#[test]
fn imrphenomxp_nnlo_matches_lalsimulation() {
    common::check_all("imrphenomxp_nnlo", 1e-9, 1e-12, |case| {
        IMRPhenomXP::with_version(&precessing_binary(case), PrecessionVersion::Nnlo)
            .unwrap()
            .polarizations(&case.frequencies)
            .unwrap()
    });
}

#[test]
fn imrphenomxp_spintaylor_matches_lalsimulation() {
    common::check_all("imrphenomxp_spintaylor", 1e-9, 1e-12, |case| {
        IMRPhenomXP::new(&precessing_binary(case))
            .unwrap()
            .polarizations(&case.frequencies)
            .unwrap()
    });
}

#[test]
fn imrphenomxp_uniform_matches_lalsimulation() {
    // SimInspiralChooseFDWaveform on the grid k * delta_f
    common::check_all("imrphenomxp_uniform", 1e-9, 1e-12, |case| {
        let version = match case.lal_flags.get("PhenomXPrecVersion") {
            Some(&flag) if flag == 102.0 => PrecessionVersion::Nnlo,
            _ => PrecessionVersion::SpinTaylor,
        };
        let delta_f = case.frequencies[1];
        IMRPhenomXP::with_version(&precessing_binary(case), version)
            .unwrap()
            .polarizations_uniform(delta_f, case.f_min, case.f_max)
            .unwrap()
    });
}

#[derive(Debug, Deserialize)]
struct AngleFile {
    lalsimulation_version: String,
    cases: Vec<AngleCase>,
}

#[derive(Debug, Deserialize)]
struct AngleCase {
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

/// The SpinTaylor Euler angles against `SimIMRPhenomXPSpinTaylorAngles`.
/// On arm64 (where the port mirrors LAL's fused multiply-adds) they agree
/// bit for bit.
#[test]
fn imrphenomxp_spintaylor_angles_match_lalsimulation() {
    let path = format!(
        "{}/tests/data/imrphenomxp_spintaylor_angles.json",
        env!("CARGO_MANIFEST_DIR")
    );
    let reference: AngleFile =
        serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap();
    for case in &reference.cases {
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
        assert_eq!(angles.alpha.len(), case.length);
        let mut max_diff = 0.0f64;
        for (test, truth) in [
            (&angles.alpha, &case.alpha),
            (&angles.cos_beta, &case.cos_beta),
            (&angles.gamma, &case.gamma),
        ] {
            for (a, b) in test.iter().step_by(case.stride).zip(truth) {
                max_diff = max_diff.max((a - b).abs());
            }
        }
        println!(
            "max |angle difference| {max_diff:.3e} | m1={} m2={} (lalsimulation {})",
            case.mass_1, case.mass_2, reference.lalsimulation_version
        );
        assert!(max_diff <= 1e-11, "{max_diff}");
    }
}

#[test]
fn imrphenomxhm_matches_lalsimulation() {
    common::check_all("imrphenomxhm", 1e-9, 1e-12, |case| {
        // The component masses are passed exactly: the (3,2) mode of
        // IMRPhenomXHM amplifies rounding differences in its inputs.
        IMRPhenomXHM::from_components(
            case.mass_1,
            case.mass_2,
            case.spin_1[2],
            case.spin_2[2],
            case.luminosity_distance,
        )
        .unwrap()
        .with_reference_frequency(case.f_ref)
        .polarizations(&case.frequencies, case.phi_ref, case.inclination)
    });
}

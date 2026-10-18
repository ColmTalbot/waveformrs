//! Comparisons against LALSimulation.
//!
//! Reference data are produced by `scripts/generate_reference_data.py` using
//! `SimInspiralChooseFDWaveformSequence` and stored in `tests/data`.
mod common;

use serde::Deserialize;
use waveformrs::{
    imrphenomd::IMRPhenomD,
    imrphenomx::{
        IMRPhenomXAS, IMRPhenomXHM, IMRPhenomXP, IMRPhenomXPHM, PrecessingBinary, PrecessionVersion,
    },
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
        .waveform_modes(&case.frequencies, case.phi_ref, case.inclination)
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
        .waveform_modes(&case.frequencies, case.phi_ref, case.inclination)
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
        let delta_f = case.frequencies[1];
        IMRPhenomXP::with_version(&precessing_binary(case), precession_version(case))
            .unwrap()
            .polarizations_uniform(delta_f, case.f_min, case.f_max)
            .unwrap()
    });
}

/// Polarisations on the non-uniform grid of the case or, for the cases
/// generated with `SimInspiralChooseFDWaveform` (grid `k delta_f` from zero),
/// on the uniform grid.
fn xp_polarizations(
    case: &common::Case,
    version: PrecessionVersion,
) -> (Vec<common::C64>, Vec<common::C64>) {
    let model = IMRPhenomXP::with_version(&precessing_binary(case), version).unwrap();
    if case.frequencies[0] == 0.0 {
        model
            .polarizations_uniform(case.frequencies[1], case.f_min, case.f_max)
            .unwrap()
    } else {
        model.polarizations(&case.frequencies).unwrap()
    }
}

/// MSA angles (`PhenomXPrecVersion` 223), including a binary whose MSA
/// system fails to initialise (NNLO fallback). The reference data are
/// generated single-threaded (see `msa.rs` on LAL's OpenMP race).
#[test]
fn imrphenomxp_msa_matches_lalsimulation() {
    common::check_all("imrphenomxp_msa", 1e-9, 1e-12, |case| {
        xp_polarizations(case, PrecessionVersion::Msa)
    });
}

/// Heavy binaries for which the SpinTaylor integration does not reach
/// `f_min` and lalsimulation falls back to the MSA angles.
#[test]
fn imrphenomxp_spintaylor_fallback_matches_lalsimulation() {
    common::check_all("imrphenomxp_spintaylor_fallback", 1e-9, 1e-12, |case| {
        xp_polarizations(case, PrecessionVersion::SpinTaylor)
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

#[derive(Deserialize)]
struct XhmModeFile {
    lalsimulation_version: String,
    cases: Vec<XhmModeCase>,
}

#[derive(Deserialize)]
struct XhmModeCase {
    mass_1: f64,
    mass_2: f64,
    chi_1: f64,
    chi_2: f64,
    luminosity_distance: f64,
    phi_ref: f64,
    f_ref: f64,
    frequencies: Vec<f64>,
    modes: Vec<XhmMode>,
}

#[derive(Deserialize)]
struct XhmMode {
    ell: u32,
    emm: i32,
    real: Vec<f64>,
    imag: Vec<f64>,
}

/// Each mode of IMRPhenomXHM against `SimIMRPhenomXHMFrequencySequenceOneMode`.
#[test]
fn imrphenomxhm_modes_match_lalsimulation() {
    let path = format!(
        "{}/tests/data/imrphenomxhm_modes.json",
        env!("CARGO_MANIFEST_DIR")
    );
    let reference: XhmModeFile =
        serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap();
    let mut failures = Vec::new();
    for case in &reference.cases {
        let xhm = IMRPhenomXHM::from_components(
            case.mass_1,
            case.mass_2,
            case.chi_1,
            case.chi_2,
            case.luminosity_distance,
        )
        .unwrap()
        .with_reference_frequency(case.f_ref);
        for mode in &case.modes {
            let test = xhm.mode(mode.ell, mode.emm, &case.frequencies, case.phi_ref);
            let truth: Vec<common::C64> = mode
                .real
                .iter()
                .zip(&mode.imag)
                .map(|(&re, &im)| common::C64::new(re, im))
                .collect();
            let result = common::compare(&case.frequencies, &test, &truth);
            println!(
                "({},{}): max rel err {:.3e}, mismatch {:.3e} | m1={} m2={} chi=({}, {})",
                mode.ell,
                mode.emm,
                result.max_relative_error,
                result.mismatch,
                case.mass_1,
                case.mass_2,
                case.chi_1,
                case.chi_2
            );
            if !(result.max_relative_error <= 1e-9 && result.mismatch <= 1e-12) {
                failures.push(format!(
                    "({},{}) {:?} m1={} m2={}",
                    mode.ell, mode.emm, result, case.mass_1, case.mass_2
                ));
            }
        }
    }
    assert!(
        failures.is_empty(),
        "modes differ from lalsimulation {}:\n{}",
        reference.lalsimulation_version,
        failures.join("\n")
    );
}

fn precession_version(case: &common::Case) -> PrecessionVersion {
    if case.lal_flags.get("PhenomXPrecVersion") == Some(&102.0) {
        PrecessionVersion::Nnlo
    } else {
        PrecessionVersion::SpinTaylor
    }
}

#[test]
fn imrphenomxphm_spintaylor_matches_lalsimulation() {
    common::check_all("imrphenomxphm_spintaylor", 1e-9, 1e-12, |case| {
        IMRPhenomXPHM::new(&precessing_binary(case))
            .unwrap()
            .polarizations(&case.frequencies)
            .unwrap()
    });
}

#[test]
fn imrphenomxphm_nnlo_matches_lalsimulation() {
    common::check_all("imrphenomxphm_nnlo", 1e-9, 1e-12, |case| {
        IMRPhenomXPHM::with_version(&precessing_binary(case), PrecessionVersion::Nnlo)
            .unwrap()
            .polarizations(&case.frequencies)
            .unwrap()
    });
}

#[test]
fn imrphenomxphm_uniform_matches_lalsimulation() {
    // SimInspiralChooseFDWaveform on the grid k * delta_f
    common::check_all("imrphenomxphm_uniform", 1e-9, 1e-12, |case| {
        IMRPhenomXPHM::with_version(&precessing_binary(case), precession_version(case))
            .unwrap()
            .polarizations_uniform(case.frequencies[1], case.f_min, case.f_max)
            .unwrap()
    });
}

#[test]
fn imrphenomxphm_mode_arrays_match_lalsimulation() {
    common::check_all("imrphenomxphm_modes", 1e-9, 1e-12, |case| {
        let mut model =
            IMRPhenomXPHM::with_version(&precessing_binary(case), precession_version(case))
                .unwrap();
        if !case.mode_array.is_empty() {
            let modes: Vec<(u32, u32)> = case
                .mode_array
                .iter()
                .filter(|(_, m)| *m > 0)
                .map(|&(l, m)| (l, m as u32))
                .collect();
            model = model.with_mode_array(&modes).unwrap();
        }
        model.polarizations(&case.frequencies).unwrap()
    });
}

fn xphm_polarizations(
    case: &common::Case,
    version: PrecessionVersion,
) -> (Vec<common::C64>, Vec<common::C64>) {
    let model = IMRPhenomXPHM::with_version(&precessing_binary(case), version).unwrap();
    // uniform grids start at zero frequency
    if case.frequencies[0] == 0.0 {
        model
            .polarizations_uniform(case.frequencies[1], case.f_min, case.f_max)
            .unwrap()
    } else {
        model.polarizations(&case.frequencies).unwrap()
    }
}

#[test]
fn imrphenomxphm_msa_matches_lalsimulation() {
    common::check_all("imrphenomxphm_msa", 1e-9, 1e-12, |case| {
        xphm_polarizations(case, PrecessionVersion::Msa)
    });
}

/// Heavy binaries for which the SpinTaylor integration does not provide the
/// angles and lalsimulation falls back to the MSA angles.
#[test]
fn imrphenomxphm_spintaylor_fallback_matches_lalsimulation() {
    common::check_all("imrphenomxphm_spintaylor_fallback", 1e-9, 1e-12, |case| {
        let result = xphm_polarizations(case, PrecessionVersion::SpinTaylor);
        assert_eq!(
            result,
            xphm_polarizations(case, PrecessionVersion::Msa),
            "no fallback for {}",
            case.label()
        );
        result
    });
}

//! Comparisons against LALSimulation.
//!
//! Reference data are produced by `scripts/generate_reference_data.py` using
//! `SimInspiralChooseFDWaveformSequence` and stored in `tests/data`.
mod common;

use waveformrs::{imrphenomd::IMRPhenomD, taylor::TaylorF2, waveform::Waveform};

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

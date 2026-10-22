// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2024-2026 Colm Talbot. Parts of waveformrs are translated
// from LALSuite, GSL and fdlibm: see NOTICE.md for their authors.
//! IMRPhenomXPNR with LAL's moving average of the antisymmetric amplitude
//! ratio against lalsimulation 6.2.0 with the fix of lalsuite !2593, at the
//! tolerances of the rest of the model.
//!
//! The model averages the ratio as Ghosh et al. 2024 describe (Eq. 17 of
//! arXiv:2310.16980), not in place as LAL does, which changes the
//! waveforms by up to ~3e-3 (docs/lal-quirks.md). These tests switch to
//! LAL's average on their thread, so that everything else is compared with
//! LAL as closely as before; `tests/xpnr.rs` compares the model itself.
use crate::test_common as common;

use crate::imrphenomx::pnr::LAL_MOVING_AVERAGE;
use crate::imrphenomx::{IMRPhenomXPNR, PrecessingBinary};

fn model_of(case: &common::Case) -> IMRPhenomXPNR {
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
    let model = IMRPhenomXPNR::new(&binary).unwrap();
    if case.mode_array.is_empty() {
        return model;
    }
    let mut modes: Vec<(u32, u32)> = case
        .mode_array
        .iter()
        .map(|&(ell, emm)| (ell, emm.unsigned_abs()))
        .collect();
    modes.sort_unstable();
    modes.dedup();
    model.with_mode_array(&modes).unwrap()
}

#[test]
fn polarizations_match_patched_lalsimulation() {
    LAL_MOVING_AVERAGE.with(|lal| lal.set(true));
    common::check_all("imrphenomxpnr_parityfix", 1e-12, 1e-15, |case| {
        model_of(case).polarizations(&case.frequencies).unwrap()
    });
}

#[test]
fn uniform_grids_match_patched_lalsimulation() {
    LAL_MOVING_AVERAGE.with(|lal| lal.set(true));
    common::check_all("imrphenomxpnr_parityfix_uniform", 1e-12, 1e-15, |case| {
        model_of(case)
            .polarizations_uniform(case.frequencies[1], case.f_min, case.f_max)
            .unwrap()
    });
}

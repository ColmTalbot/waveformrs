// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2024-2026 Colm Talbot. Parts of waveformrs are translated
// from LALSuite, GSL and fdlibm: see NOTICE.md for their authors.
//! The cost of `NRSur7dq4::with_continuous`: median time of one waveform on
//! a detector-like grid, alternating LAL's and the continuous version.
//!
//! cargo run --release --example nrsur_continuous_timing

use std::path::PathBuf;
use std::time::Instant;

use waveformrs::nrsur::{NRSur7dq4, NRSur7dq4Data, NRSurWorkspace, TimeSeriesPolarizations};
use waveformrs::Binary;

fn main() {
    let path = std::env::var_os("NRSUR7DQ4_DATA")
        .map(PathBuf::from)
        .unwrap_or_else(|| {
            PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("lal_data/NRSur7dq4_v1.0.h5")
        });
    let data = NRSur7dq4Data::load(&path).unwrap();
    let mut ws = NRSurWorkspace::new();
    let mut outs = vec![TimeSeriesPolarizations::default(); 2];
    let mut times = [Vec::new(), Vec::new()];
    for i in 0..400 {
        let b = Binary {
            mass_1: 36.0 + 0.01 * i as f64,
            mass_2: 29.0,
            spin_1: [0.3, 0.1, 0.25],
            spin_2: [-0.1, 0.25, 0.1],
            luminosity_distance: 500.0,
            inclination: 0.6,
            phi_ref: 1.3,
            f_ref: 20.0,
        };
        for (k, continuous) in [false, true].into_iter().enumerate() {
            let model = NRSur7dq4::new(&data, &b)
                .unwrap()
                .with_continuous(continuous);
            let t = Instant::now();
            model
                .polarizations_on_grids_with(
                    &mut ws,
                    0.0,
                    1.0 / 2048.0,
                    &[-0.6, -0.59],
                    1434,
                    &mut outs,
                )
                .unwrap();
            times[k].push(t.elapsed().as_secs_f64() * 1e6);
        }
    }
    for (k, name) in ["LAL", "continuous"].iter().enumerate() {
        times[k].sort_by(f64::total_cmp);
        println!("{name}: median {:.0} us", times[k][times[k].len() / 2]);
    }
}

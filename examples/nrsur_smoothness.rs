// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2024-2026 Colm Talbot. Parts of waveformrs are translated
// from LALSuite, GSL and fdlibm: see NOTICE.md for their authors.
//! Smoothness of NRSur7dq4 in its parameters: `h_+` on a fixed time grid
//! around the peak along fine steps of one parameter, written as the norm
//! of the second differences relative to the norm of `h_+` (steady for a
//! smooth model, with spikes at jumps).
//!
//! cargo run --release --example nrsur_smoothness -- PARAMETER STEP [N] [OUTPUT]
//!
//! PARAMETER is one of total_mass, mass_ratio, chi1x, chi1z, chi2y, f_ref.

use std::path::PathBuf;

use waveformrs::nrsur::{NRSur7dq4, NRSur7dq4Data, NRSurWorkspace, TimeSeriesPolarizations};
use waveformrs::Binary;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let parameter = args.get(1).map_or("total_mass", |s| s.as_str()).to_string();
    let step: f64 = args.get(2).map_or(1e-4, |s| s.parse().unwrap());
    let n: usize = args.get(3).map_or(1000, |s| s.parse().unwrap());
    let output = args.get(4).cloned();
    let path = std::env::var_os("NRSUR7DQ4_DATA")
        .map(PathBuf::from)
        .unwrap_or_else(|| {
            PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("lal_data/NRSur7dq4_v1.0.h5")
        });
    let data = NRSur7dq4Data::load(&path).unwrap();
    let base = Binary {
        mass_1: 36.0,
        mass_2: 29.0,
        spin_1: [0.3, 0.1, 0.25],
        spin_2: [-0.1, 0.25, 0.1],
        luminosity_distance: 500.0,
        inclination: 0.6,
        phi_ref: 1.3,
        f_ref: 20.0,
    };
    let (delta_t, length) = (1.0 / 4096.0, 4096);
    let starts = [-0.8];
    let mut ws = NRSurWorkspace::new();
    let mut outs = vec![TimeSeriesPolarizations::default()];
    let waveforms: Vec<Vec<f64>> = (0..n)
        .map(|i| {
            let x = i as f64 * step;
            let mut b = base;
            match parameter.as_str() {
                "total_mass" => {
                    let s = 1.0 + x / 65.0;
                    b.mass_1 *= s;
                    b.mass_2 *= s;
                }
                "mass_ratio" => b.mass_2 += x,
                "chi1x" => b.spin_1[0] += x,
                "chi1z" => b.spin_1[2] += x,
                "chi2y" => b.spin_2[1] += x,
                "f_ref" => b.f_ref += x,
                other => panic!("unknown parameter {other}"),
            }
            let model = NRSur7dq4::new(&data, &b)
                .unwrap()
                .with_continuous(std::env::var("NRSUR_CONTINUOUS").is_ok());
            model
                .polarizations_on_grids_with(&mut ws, 0.0, delta_t, &starts, length, &mut outs)
                .unwrap();
            outs[0].hplus.clone()
        })
        .collect();
    let norm = waveforms[0].iter().map(|h| h * h).sum::<f64>().sqrt();
    let second: Vec<f64> = waveforms
        .windows(3)
        .map(|w| {
            let s: f64 = (0..length)
                .map(|k| (w[0][k] - 2. * w[1][k] + w[2][k]).powi(2))
                .sum();
            s.sqrt() / norm
        })
        .collect();
    let mut sorted = second.clone();
    sorted.sort_by(f64::total_cmp);
    let median = sorted[sorted.len() / 2];
    let spikes: Vec<usize> = (0..second.len())
        .filter(|&i| second[i] > 5.0 * median)
        .collect();
    println!(
        "{parameter}: step {step:e}, median {median:.3e}, max {:.3e}, {} spikes at {:?}",
        sorted[sorted.len() - 1],
        spikes.len(),
        spikes
            .iter()
            .take(10)
            .map(|&i| (i + 1) as f64 * step)
            .collect::<Vec<_>>()
    );
    if let Some(&worst) = spikes.get(
        std::env::var("SPIKE")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(0),
    ) {
        let w = &waveforms[worst..worst + 3];
        let mut per: Vec<(f64, usize)> = (0..length)
            .map(|k| ((w[0][k] - 2. * w[1][k] + w[2][k]).abs(), k))
            .collect();
        per.sort_by(|a, b| b.0.total_cmp(&a.0));
        for (v, k) in per.iter().take(6) {
            println!(
                "  sample {k} (t = {:.5} s) second difference {v:.3e}, |h| {:.3e}",
                starts[0] + *k as f64 * delta_t,
                w[1][*k].abs()
            );
        }
    }
    if let Some(output) = output {
        let mut text = String::new();
        for (i, v) in second.iter().enumerate() {
            text += &format!("{} {v:e}\n", (i + 1) as f64 * step);
        }
        std::fs::write(output, text).unwrap();
    }
}

// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2024-2026 Colm Talbot. Parts of waveformrs are translated
// from LALSuite, GSL and fdlibm: see NOTICE.md for their authors.
//! `PrecessionVersion::SpinTaylorContinuous`: IMRPhenomXPHM changes smoothly
//! with the masses, where the SpinTaylor angles of LAL (320) jump whenever
//! the PN solution gains or loses a sample.

use waveformrs::{
    choose_fd_model, imrphenomx::PrecessionVersion, Approximant, Binary, WaveformOptions,
};

/// The binary with chirp mass `chirp_mass` and mass ratio 0.8.
fn binary(chirp_mass: f64) -> Binary {
    let q: f64 = 0.8;
    let total = chirp_mass * (1. + q).powf(1.2) / q.powf(0.6);
    let (t1, t2, phi12): (f64, f64, f64) = (2.1, 1.65, 2.07);
    Binary {
        mass_1: total / (1. + q),
        mass_2: total * q / (1. + q),
        spin_1: [0.54 * t1.sin(), 0.0, 0.54 * t1.cos()],
        spin_2: [
            0.03 * t2.sin() * phi12.cos(),
            0.03 * t2.sin() * phi12.sin(),
            0.03 * t2.cos(),
        ],
        luminosity_distance: 500.0,
        inclination: 0.5,
        phi_ref: 2.6,
        f_ref: 20.0,
    }
}

/// The largest second difference of `h_+` along steps of `4e-5` solar
/// masses in the chirp mass, relative to the norm of `h_+`, and the median.
/// Over the 0.03 solar masses swept, the curvature of the waveform is
/// nearly constant, so a jump in the waveform of more than a few 1e-9 (in
/// the model's angles or anywhere else) stands out.
fn second_differences(approximant: Approximant, version: PrecessionVersion) -> (f64, f64) {
    let frequencies: Vec<f64> = (80..4096).map(|i| i as f64 * 0.25).collect();
    let options = WaveformOptions::default().with_precession_version(version);
    let waveforms: Vec<_> = (0..800)
        .map(|i| {
            let model =
                choose_fd_model(approximant, &binary(27.95 + 4e-5 * i as f64), &options).unwrap();
            model.polarizations(&frequencies).unwrap().0
        })
        .collect();
    let norm = waveforms[0]
        .iter()
        .map(|h| h.norm_sqr())
        .sum::<f64>()
        .sqrt();
    let mut values: Vec<f64> = waveforms
        .windows(3)
        .map(|w| {
            let sum: f64 = (0..frequencies.len())
                .map(|k| (w[0][k] - 2. * w[1][k] + w[2][k]).norm_sqr())
                .sum();
            sum.sqrt() / norm
        })
        .collect();
    values.sort_by(f64::total_cmp);
    (values[values.len() - 1], values[values.len() / 2])
}

#[test]
fn smooth_in_the_chirp_mass() {
    for approximant in [Approximant::IMRPhenomXP, Approximant::IMRPhenomXPHM] {
        let (lal_max, lal_median) = second_differences(approximant, PrecessionVersion::SpinTaylor);
        let (max, median) =
            second_differences(approximant, PrecessionVersion::SpinTaylorContinuous);
        println!("{approximant:?} SpinTaylor: max {lal_max:.2e} median {lal_median:.2e}");
        println!("{approximant:?} continuous: max {max:.2e} median {median:.2e}");
        // LAL's angles jump every ~0.01 solar masses in the chirp mass.
        assert!(lal_max > 1e3 * lal_median, "{lal_max:e} {lal_median:e}");
        // The continuous angles change smoothly: steady second differences
        // (the end of the PN solution, the start of the merger-ringdown
        // continuation and the integration of gamma across it all move
        // continuously).
        assert!(max < 1.05 * median, "{max:e} {median:e}");
    }
}

/// The continuous angles agree with LAL's to about the size of LAL's jumps.
#[test]
fn close_to_the_spintaylor_angles() {
    let frequencies: Vec<f64> = (80..4096).map(|i| i as f64 * 0.25).collect();
    let mut worst: f64 = 0.0;
    for i in 0..12 {
        let b = binary(28.0 + 1e-3 * i as f64);
        let h = |version| {
            let options = WaveformOptions::default().with_precession_version(version);
            choose_fd_model(Approximant::IMRPhenomXPHM, &b, &options)
                .unwrap()
                .polarizations(&frequencies)
                .unwrap()
                .0
        };
        let (lal, continuous) = (
            h(PrecessionVersion::SpinTaylor),
            h(PrecessionVersion::SpinTaylorContinuous),
        );
        let norm = lal.iter().map(|h| h.norm_sqr()).sum::<f64>().sqrt();
        let difference = lal
            .iter()
            .zip(&continuous)
            .map(|(a, b)| (a - b).norm_sqr())
            .sum::<f64>()
            .sqrt();
        worst = worst.max(difference / norm);
    }
    println!("largest relative difference {worst:.2e}");
    assert!(worst < 1e-2, "{worst:e}");
}

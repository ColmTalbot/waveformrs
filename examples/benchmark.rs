//! Timing of waveform generation for the models in the crate.
//!
//! Each model is constructed and evaluated on a uniform 0.25 Hz grid between
//! 20 and 1500 Hz (5921 frequencies), the same set-up as `examples/tf2.rs`.
//!
//!     cargo run --release --example benchmark
use std::hint::black_box;
use std::time::Instant;

use waveformrs::{
    imrphenomd::IMRPhenomD, imrphenomx::IMRPhenomXAS, taylor::TaylorF2, waveform::Waveform,
};

fn frequencies() -> Vec<f64> {
    let delta_f = 0.25;
    (80..=6000).map(|i| i as f64 * delta_f).collect()
}

fn time<F: FnMut()>(name: &str, n_frequencies: usize, mut func: F) {
    // warm up
    for _ in 0..20 {
        func();
    }
    let mut best = f64::INFINITY;
    let mut total = 0.0;
    let repeats = 10;
    let iterations = 200;
    for _ in 0..repeats {
        let start = Instant::now();
        for _ in 0..iterations {
            func();
        }
        let elapsed = start.elapsed().as_secs_f64() / iterations as f64;
        best = best.min(elapsed);
        total += elapsed;
    }
    println!(
        "{name:<32} best {:>9.2} us  mean {:>9.2} us  ({:.2} ns / frequency)",
        best * 1e6,
        total / repeats as f64 * 1e6,
        best * 1e9 / n_frequencies as f64
    );
}

fn main() {
    let freqs = frequencies();
    let n = freqs.len();
    let (total_mass, mass_ratio, chi_1, chi_2, distance) = (60.0, 0.5, 0.3, -0.2, 400.0);

    time("TaylorF2 new", 1, || {
        black_box(TaylorF2::new(total_mass, mass_ratio, chi_1, chi_2, distance));
    });
    time("TaylorF2 new + waveform", n, || {
        let wf = TaylorF2::new(total_mass, mass_ratio, chi_1, chi_2, distance);
        black_box(wf.waveform(black_box(&freqs), 0.1));
    });
    time("TaylorF2 new + waveform_modes", n, || {
        let wf = TaylorF2::new(total_mass, mass_ratio, chi_1, chi_2, distance);
        black_box(wf.waveform_modes(black_box(&freqs), 0.1, 0.4));
    });
    time("IMRPhenomD new", 1, || {
        black_box(IMRPhenomD::new(total_mass, mass_ratio, chi_1, chi_2, distance));
    });
    time("IMRPhenomD new + waveform", n, || {
        let wf = IMRPhenomD::new(total_mass, mass_ratio, chi_1, chi_2, distance);
        black_box(wf.waveform(black_box(&freqs), 0.1));
    });
    time("IMRPhenomD new + waveform_modes", n, || {
        let wf = IMRPhenomD::new(total_mass, mass_ratio, chi_1, chi_2, distance);
        black_box(wf.waveform_modes(black_box(&freqs), 0.1, 0.4));
    });
    time("IMRPhenomXAS new", 1, || {
        black_box(IMRPhenomXAS::new(total_mass, mass_ratio, chi_1, chi_2, distance));
    });
    time("IMRPhenomXAS new + waveform_modes", n, || {
        let wf = IMRPhenomXAS::new(total_mass, mass_ratio, chi_1, chi_2, distance)
            .with_reference_frequency(20.0);
        black_box(wf.waveform_modes(black_box(&freqs), 0.1, 0.4));
    });
}

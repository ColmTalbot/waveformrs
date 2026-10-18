//! Per-frequency evaluation against evaluation on a prepared
//! [`FrequencyGrid`] (as when generating many waveforms on the same
//! frequencies): timings interleaved so that machine load affects both alike,
//! and the largest difference between the two.
//!
//!     cargo run --release --example grid_benchmark
use num_complex::Complex;
use std::hint::black_box;
use std::time::Instant;
use waveformrs::{
    grid::FrequencyGrid,
    imrphenomd::IMRPhenomD,
    imrphenomx::{IMRPhenomXAS, IMRPhenomXHM},
    taylor::TaylorF2,
    waveform::Waveform,
};

fn ab<W: Waveform>(name: &str, make: impl Fn(usize) -> W, grid: &FrequencyGrid, phi_c: f64) {
    let n = grid.len();
    let (mut hp, mut hc) = (
        vec![Complex::new(0.0, 0.0); n],
        vec![Complex::new(0.0, 0.0); n],
    );
    let (mut t_old, mut t_new) = (f64::INFINITY, f64::INFINITY);
    for rep in 0..30 {
        let start = Instant::now();
        for k in 0..50 {
            black_box(make(k).waveform_modes(black_box(grid.frequencies()), phi_c, 0.4));
        }
        t_old = t_old.min(start.elapsed().as_secs_f64() / 50.0);
        let start = Instant::now();
        for k in 0..50 {
            make(k).waveform_modes_on_grid(black_box(grid), phi_c, 0.4, &mut hp, &mut hc);
            black_box(&hp);
        }
        t_new = t_new.min(start.elapsed().as_secs_f64() / 50.0);
        let _ = rep;
    }
    let (p_old, _) = make(3).waveform_modes(grid.frequencies(), phi_c, 0.4);
    make(3).waveform_modes_on_grid(grid, phi_c, 0.4, &mut hp, &mut hc);
    let scale = p_old.iter().map(|x| x.norm()).fold(0.0, f64::max);
    let diff = p_old
        .iter()
        .zip(&hp)
        .map(|(a, b)| (a - b).norm())
        .fold(0.0, f64::max)
        / scale;
    println!(
        "{name:<14} per-frequency {:7.1} us   grid {:7.1} us   x{:.2}   max rel diff {diff:.1e}",
        t_old * 1e6,
        t_new * 1e6,
        t_old / t_new
    );
}

fn main() {
    let grid = FrequencyGrid::uniform(20.0, 1024.0, 0.125);
    let m = |k: usize| 30.0 + 2.0 * k as f64;
    ab(
        "TaylorF2",
        |k| TaylorF2::new(m(k), 0.5, 0.3, -0.2, 400.0).with_reference_frequency(20.0),
        &grid,
        0.1,
    );
    ab(
        "IMRPhenomD",
        |k| IMRPhenomD::new(m(k), 0.5, 0.3, -0.2, 400.0).with_reference_frequency(20.0),
        &grid,
        0.2,
    );
    ab(
        "IMRPhenomXAS",
        |k| IMRPhenomXAS::new(m(k), 0.5, 0.3, -0.2, 400.0).with_reference_frequency(20.0),
        &grid,
        0.2,
    );
    // IMRPhenomXHM
    let n = grid.len();
    let (mut hp, mut hc) = (
        vec![Complex::new(0.0, 0.0); n],
        vec![Complex::new(0.0, 0.0); n],
    );
    let make =
        |k: usize| IMRPhenomXHM::new(m(k), 0.5, 0.3, -0.2, 400.0).with_reference_frequency(20.0);
    let (mut t_old, mut t_new) = (f64::INFINITY, f64::INFINITY);
    for _ in 0..10 {
        let start = Instant::now();
        for k in 0..20 {
            black_box(make(k).polarizations(black_box(grid.frequencies()), 0.1, 0.4));
        }
        t_old = t_old.min(start.elapsed().as_secs_f64() / 20.0);
        let start = Instant::now();
        for k in 0..20 {
            make(k).polarizations_on_grid(black_box(&grid), 0.1, 0.4, &mut hp, &mut hc);
            black_box(&hp);
        }
        t_new = t_new.min(start.elapsed().as_secs_f64() / 20.0);
    }
    let (p_old, _) = make(3).polarizations(grid.frequencies(), 0.1, 0.4);
    make(3).polarizations_on_grid(&grid, 0.1, 0.4, &mut hp, &mut hc);
    let scale = p_old.iter().map(|x| x.norm()).fold(0.0, f64::max);
    let diff = p_old
        .iter()
        .zip(&hp)
        .map(|(a, b)| (a - b).norm())
        .fold(0.0, f64::max)
        / scale;
    println!(
        "{:<14} per-frequency {:7.1} us   grid {:7.1} us   x{:.2}   max rel diff {diff:.1e}",
        "IMRPhenomXHM",
        t_old * 1e6,
        t_new * 1e6,
        t_old / t_new
    );
}

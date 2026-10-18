//! Per-frequency evaluation against evaluation on a prepared
//! [`FrequencyGrid`] (as when generating many waveforms on the same
//! frequencies): timings interleaved so that machine load affects both alike,
//! and the largest difference between the two.
//!
//! The precessing models use the binary of `examples/benchmark.rs`
//! (m1 = 40, m2 = 20 Msun, chi1 = (0.4, -0.3, 0.2), chi2 = (0.1, 0.5, -0.3)).
//!
//!     cargo run --release --example grid_benchmark
use num_complex::Complex;
use std::hint::black_box;
use std::time::Instant;
use waveformrs::{
    grid::FrequencyGrid,
    imrphenomd::IMRPhenomD,
    imrphenomx::{
        IMRPhenomXAS, IMRPhenomXHM, IMRPhenomXP, IMRPhenomXPHM, PrecessingBinary, PrecessionVersion,
    },
    taylor::TaylorF2,
    waveform::Waveform,
};

type Polarizations = (Vec<Complex<f64>>, Vec<Complex<f64>>);

/// Time `per_frequency(k)` against `on_grid(k, hplus, hcross)` for waveforms
/// `k = 0, 1, ...` (best of `repeats` batches of `batch` waveforms each) and
/// print the timings and the largest difference relative to the largest
/// per-frequency value.
fn ab(
    name: &str,
    grid: &FrequencyGrid,
    (repeats, batch): (usize, usize),
    per_frequency: impl Fn(usize) -> Polarizations,
    on_grid: impl Fn(usize, &mut [Complex<f64>], &mut [Complex<f64>]),
) {
    let n = grid.len();
    let (mut hp, mut hc) = (
        vec![Complex::new(0.0, 0.0); n],
        vec![Complex::new(0.0, 0.0); n],
    );
    let (mut t_old, mut t_new) = (f64::INFINITY, f64::INFINITY);
    // alternate the order: running one right after the other can slow it
    // down (seen for IMRPhenomD on the grid after the per-frequency path)
    for repeat in 0..repeats {
        for path in [repeat % 2, 1 - repeat % 2] {
            let start = Instant::now();
            for k in 0..batch {
                if path == 0 {
                    black_box(per_frequency(k));
                } else {
                    on_grid(k, &mut hp, &mut hc);
                    black_box((&hp, &hc));
                }
            }
            let elapsed = start.elapsed().as_secs_f64() / batch as f64;
            if path == 0 {
                t_old = t_old.min(elapsed);
            } else {
                t_new = t_new.min(elapsed);
            }
        }
    }
    let (p_old, c_old) = per_frequency(3);
    on_grid(3, &mut hp, &mut hc);
    let scale = p_old
        .iter()
        .chain(&c_old)
        .map(|x| x.norm())
        .fold(0.0, f64::max);
    let diff = p_old
        .iter()
        .zip(&hp)
        .chain(c_old.iter().zip(&hc))
        .map(|(a, b)| (a - b).norm())
        .fold(0.0, f64::max)
        / scale;
    println!(
        "{name:<19} per-frequency {:7.1} us   grid {:7.1} us   x{:.2}   max rel diff {diff:.1e}",
        t_old * 1e6,
        t_new * 1e6,
        t_old / t_new
    );
}

/// [`ab`] for a model with the [`Waveform`] interface.
fn ab_waveform<W: Waveform>(
    name: &str,
    make: impl Fn(usize) -> W,
    grid: &FrequencyGrid,
    phi_c: f64,
) {
    ab(
        name,
        grid,
        (30, 50),
        |k| make(k).waveform_modes(black_box(grid.frequencies()), phi_c, 0.4),
        |k, hp, hc| make(k).waveform_modes_on_grid(black_box(grid), phi_c, 0.4, hp, hc),
    );
}

fn main() {
    let grid = FrequencyGrid::uniform(20.0, 1024.0, 0.125);
    let m = |k: usize| 30.0 + 2.0 * k as f64;
    ab_waveform(
        "TaylorF2",
        |k| TaylorF2::new(m(k), 0.5, 0.3, -0.2, 400.0).with_reference_frequency(20.0),
        &grid,
        0.1,
    );
    ab_waveform(
        "IMRPhenomD",
        |k| IMRPhenomD::new(m(k), 0.5, 0.3, -0.2, 400.0).with_reference_frequency(20.0),
        &grid,
        0.2,
    );
    ab_waveform(
        "IMRPhenomXAS",
        |k| IMRPhenomXAS::new(m(k), 0.5, 0.3, -0.2, 400.0).with_reference_frequency(20.0),
        &grid,
        0.2,
    );
    let xhm =
        |k: usize| IMRPhenomXHM::new(m(k), 0.5, 0.3, -0.2, 400.0).with_reference_frequency(20.0);
    ab(
        "IMRPhenomXHM",
        &grid,
        (10, 20),
        |k| xhm(k).polarizations(black_box(grid.frequencies()), 0.1, 0.4),
        |k, hp, hc| xhm(k).polarizations_on_grid(black_box(&grid), 0.1, 0.4, hp, hc),
    );

    let binary = PrecessingBinary {
        mass_1: 40.0,
        mass_2: 20.0,
        spin_1: [0.4, -0.3, 0.2],
        spin_2: [0.1, 0.5, -0.3],
        luminosity_distance: 400.0,
        inclination: 0.4,
        phi_ref: 0.1,
        f_ref: 20.0,
    };
    for (label, version) in [
        ("102", PrecessionVersion::Nnlo),
        ("223", PrecessionVersion::Msa),
        ("320", PrecessionVersion::SpinTaylor),
    ] {
        let xp = |_| IMRPhenomXP::with_version(&binary, version).unwrap();
        ab(
            &format!("IMRPhenomXP ({label})"),
            &grid,
            (10, 20),
            |k| xp(k).polarizations(black_box(grid.frequencies())).unwrap(),
            |k, hp, hc| {
                xp(k)
                    .polarizations_on_grid(black_box(&grid), hp, hc)
                    .unwrap()
            },
        );
    }
    for (label, version) in [
        ("102", PrecessionVersion::Nnlo),
        ("223", PrecessionVersion::Msa),
        ("320", PrecessionVersion::SpinTaylor),
    ] {
        let xphm = |_| IMRPhenomXPHM::with_version(&binary, version).unwrap();
        ab(
            &format!("IMRPhenomXPHM ({label})"),
            &grid,
            (10, 10),
            |k| {
                xphm(k)
                    .polarizations(black_box(grid.frequencies()))
                    .unwrap()
            },
            |k, hp, hc| {
                xphm(k)
                    .polarizations_on_grid(black_box(&grid), hp, hc)
                    .unwrap()
            },
        );
    }
}

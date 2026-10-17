//! Time the SpinTaylorT4 orbit evolution for IMRPhenomXPHM-like calls.
//!
//! Run with `cargo run --release --example spintaylor`.  The equivalent
//! lalsimulation timing (`SimInspiralSpinTaylorPNEvolveOrbit` through the
//! SWIG Python bindings) is printed by `scripts/time_spintaylor_lal.py`.
use std::hint::black_box;
use std::time::Instant;

use waveformrs::spintaylor::SpinTaylorT4;

const MTSUN_SI: f64 = 4.925490947641267e-6;

struct Case {
    label: &'static str,
    m1: f64,
    m2: f64,
    spin_1: [f64; 3],
    spin_2: [f64; 3],
    f_min: f64,
    f_ref: f64,
}

/// The calls IMRPhenomXPHM (PrecVersion 320, default coarse factor 10)
/// makes: a coarse integration from f_ref (backwards to f_min - 0.5 and
/// forwards to fCut, or only forwards if f_ref = f_min) followed by a fine
/// integration from the transition sample to fCut.
fn xphm_like(case: &Case) -> usize {
    let binary = SpinTaylorT4::from_solar_masses(case.m1, case.m2, case.spin_1, case.spin_2);
    let total_mass = case.m1 + case.m2;
    // fCut = (fRING + 8 fDAMP) / M, roughly 0.2 / M
    let f_cut = 0.2 / (total_mass * MTSUN_SI);
    let delta_t_coarse = 0.5 * 10.0 / f_cut;
    let mut n = 0;
    let coarse = if case.f_ref > case.f_min {
        let backward = binary
            .evolve_orbit(delta_t_coarse, case.f_ref, case.f_min - 0.5)
            .unwrap();
        let forward = binary
            .evolve_orbit(delta_t_coarse, case.f_ref, f_cut)
            .unwrap();
        n += backward.len();
        forward
    } else {
        binary
            .evolve_orbit(delta_t_coarse, case.f_min, f_cut)
            .unwrap()
    };
    n += coarse.len();
    let low = coarse.len();
    let index = low - 1 - 9.min(low - 1);
    let v = coarse.v[index];
    let f_trans = v * v * v / (std::f64::consts::PI * total_mass * MTSUN_SI);
    let fine = SpinTaylorT4::from_solar_masses(
        case.m1,
        case.m2,
        [coarse.s1_x[index], coarse.s1_y[index], coarse.s1_z[index]],
        [coarse.s2_x[index], coarse.s2_y[index], coarse.s2_z[index]],
    )
    .with_orientation(
        [
            coarse.lnhat_x[index],
            coarse.lnhat_y[index],
            coarse.lnhat_z[index],
        ],
        [coarse.e1_x[index], coarse.e1_y[index], coarse.e1_z[index]],
    )
    .evolve_orbit(0.5 / f_cut, f_trans, f_cut)
    .unwrap();
    n + fine.len()
}

fn time<F: FnMut() -> usize>(label: &str, mut func: F) {
    let samples = func();
    // aim for roughly a second of work
    let start = Instant::now();
    let mut repeats = 0;
    while start.elapsed().as_secs_f64() < 0.2 {
        black_box(func());
        repeats += 1;
    }
    let repeats = (repeats * 5).max(5);
    let start = Instant::now();
    for _ in 0..repeats {
        black_box(func());
    }
    let per_call = start.elapsed().as_secs_f64() / repeats as f64;
    println!(
        "{label:<45} {:>9.1} us/call  ({samples} samples, {repeats} repeats)",
        per_call * 1e6
    );
}

fn main() {
    let cases = [
        Case {
            label: "36+29 Msun, f_min = f_ref = 20 Hz",
            m1: 36.0,
            m2: 29.0,
            spin_1: [0.4, 0.2, 0.3],
            spin_2: [-0.2, 0.3, -0.1],
            f_min: 20.0,
            f_ref: 20.0,
        },
        Case {
            label: "60+15 Msun, f_min = 15, f_ref = 20 Hz",
            m1: 60.0,
            m2: 15.0,
            spin_1: [0.6, -0.3, 0.5],
            spin_2: [0.1, 0.1, 0.0],
            f_min: 15.0,
            f_ref: 20.0,
        },
        Case {
            label: "12+6 Msun, f_min = f_ref = 20 Hz",
            m1: 12.0,
            m2: 6.0,
            spin_1: [0.2, 0.7, -0.1],
            spin_2: [0.3, -0.3, 0.6],
            f_min: 20.0,
            f_ref: 20.0,
        },
    ];
    for case in &cases {
        time(case.label, || xphm_like(case));
    }

    // A single call, matching scripts/time_spintaylor_lal.py
    let binary = SpinTaylorT4::from_solar_masses(36.0, 29.0, [0.4, 0.2, 0.3], [-0.2, 0.3, -0.1]);
    let f_cut = 0.2 / (65.0 * MTSUN_SI);
    time("single evolve_orbit 36+29, 20 Hz -> fCut", || {
        binary.evolve_orbit(5.0 / f_cut, 20.0, f_cut).unwrap().len()
    });
    time("single evolve_orbit 36+29, 20 Hz -> fCut, dt/10", || {
        binary.evolve_orbit(0.5 / f_cut, 20.0, f_cut).unwrap().len()
    });
}

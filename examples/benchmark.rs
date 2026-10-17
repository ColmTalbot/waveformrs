//! Timing of waveform generation for the models in the crate.
//!
//! Each model is constructed and evaluated on a uniform 0.25 Hz grid between
//! 20 and 1500 Hz (5921 frequencies). The aligned-spin models use
//! M = 60 Msun, q = 0.5, chi = (0.3, -0.2); the precessing models use
//! m1 = 40, m2 = 20 Msun with chi1 = (0.4, -0.3, 0.2), chi2 = (0.1, 0.5, -0.3).
//! `report/lal_benchmark.py` times lalsimulation on the same grid.
//!
//!     cargo run --release --example benchmark [-- --json]
use std::hint::black_box;
use std::time::Instant;

use waveformrs::{
    imrphenomd::IMRPhenomD,
    imrphenomx::{
        IMRPhenomXAS, IMRPhenomXHM, IMRPhenomXP, PrecessingBinary, PrecessionVersion,
    },
    taylor::TaylorF2,
    waveform::Waveform,
};

fn frequencies() -> Vec<f64> {
    let delta_f = 0.25;
    (80..=6000).map(|i| i as f64 * delta_f).collect()
}

/// Best and mean time per call over 10 repeats of a batch sized to take
/// roughly 0.2 s.
fn time<F: FnMut()>(mut func: F) -> (f64, f64) {
    let start = Instant::now();
    let mut warmup = 0;
    while warmup < 5 || start.elapsed().as_secs_f64() < 0.05 {
        func();
        warmup += 1;
    }
    let per_call = start.elapsed().as_secs_f64() / warmup as f64;
    let iterations = ((0.2 / per_call) as usize).clamp(5, 2000);
    let mut best = f64::INFINITY;
    let mut total = 0.0;
    let repeats = 10;
    for _ in 0..repeats {
        let start = Instant::now();
        for _ in 0..iterations {
            func();
        }
        let elapsed = start.elapsed().as_secs_f64() / iterations as f64;
        best = best.min(elapsed);
        total += elapsed;
    }
    (best, total / repeats as f64)
}

fn main() {
    let json = std::env::args().any(|arg| arg == "--json");
    let freqs = frequencies();
    let n = freqs.len();
    let (total_mass, mass_ratio, chi_1, chi_2, distance) = (60.0, 0.5, 0.3, -0.2, 400.0);
    let (f_ref, phi_ref, inclination) = (20.0, 0.1, 0.4);
    let binary = PrecessingBinary {
        mass_1: 40.0,
        mass_2: 20.0,
        spin_1: [0.4, -0.3, 0.2],
        spin_2: [0.1, 0.5, -0.3],
        luminosity_distance: distance,
        inclination,
        phi_ref,
        f_ref,
    };

    let mut results: Vec<(String, (f64, f64))> = Vec::new();
    let mut run = |name: &str, func: &mut dyn FnMut()| {
        let timing = time(func);
        if !json {
            println!(
                "{name:<32} best {:>9.2} us  mean {:>9.2} us  ({:.2} ns / frequency)",
                timing.0 * 1e6,
                timing.1 * 1e6,
                timing.0 * 1e9 / n as f64
            );
        }
        results.push((name.to_string(), timing));
    };

    run("TaylorF2", &mut || {
        let wf = TaylorF2::new(total_mass, mass_ratio, chi_1, chi_2, distance)
            .with_reference_frequency(f_ref);
        black_box(wf.waveform_modes(black_box(&freqs), phi_ref, inclination));
    });
    run("IMRPhenomD", &mut || {
        let wf = IMRPhenomD::new(total_mass, mass_ratio, chi_1, chi_2, distance)
            .with_reference_frequency(f_ref);
        black_box(wf.waveform_modes(black_box(&freqs), 2.0 * phi_ref, inclination));
    });
    run("IMRPhenomXAS", &mut || {
        let wf = IMRPhenomXAS::new(total_mass, mass_ratio, chi_1, chi_2, distance)
            .with_reference_frequency(f_ref);
        black_box(wf.waveform_modes(black_box(&freqs), 2.0 * phi_ref, inclination));
    });
    run("IMRPhenomXHM", &mut || {
        let wf = IMRPhenomXHM::new(total_mass, mass_ratio, chi_1, chi_2, distance)
            .with_reference_frequency(f_ref);
        black_box(wf.polarizations(black_box(&freqs), phi_ref, inclination));
    });
    for (name, version) in [
        ("IMRPhenomXP (102)", PrecessionVersion::Nnlo),
        ("IMRPhenomXP (223)", PrecessionVersion::Msa),
        ("IMRPhenomXP (320)", PrecessionVersion::SpinTaylor),
    ] {
        run(name, &mut || {
            let wf = IMRPhenomXP::with_version(&binary, version).unwrap();
            black_box(wf.polarizations(black_box(&freqs)).unwrap());
        });
    }

    if json {
        let entries: Vec<String> = results
            .iter()
            .map(|(name, (best, mean))| {
                format!("\"{name}\": {{\"best\": {best:e}, \"mean\": {mean:e}}}")
            })
            .collect();
        println!("{{{}}}", entries.join(", "));
    }
}

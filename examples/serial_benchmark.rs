//! Throughput of generating many waveforms one after another on a fixed
//! frequency grid, as a sampler does.
//!
//! Reads `{"frequencies": [...], "draws": [{"mass_1", "mass_2", "spin_1",
//! "spin_2", "luminosity_distance", "inclination", "phi_ref", "f_ref"}, ...]}`
//! from the JSON file given as the first argument (written by
//! `report/serial_benchmark.py`) and prints the mean time per waveform for
//! each model as JSON.
//!
//! Models with a grid path are also timed on a prepared
//! [`waveformrs::grid::FrequencyGrid`] with output buffers reused across
//! draws (`grid_seconds`).
//!
//!     cargo run --release --example serial_benchmark -- draws.json
use std::hint::black_box;
use std::time::Instant;

use num_complex::Complex;
use serde::Deserialize;

use waveformrs::{
    grid::FrequencyGrid,
    imrphenomd::IMRPhenomD,
    imrphenomx::{
        IMRPhenomXAS, IMRPhenomXError, IMRPhenomXHM, IMRPhenomXP, IMRPhenomXPHM, IMRPhenomXPNR,
        PrecessingBinary, PrecessionVersion,
    },
    nrsur::{NRSur7dq4, NRSur7dq4Data, NRSurWorkspace, TimeSeriesPolarizations},
    taylor::TaylorF2,
    waveform::Waveform,
};

#[derive(Deserialize)]
struct Input {
    frequencies: Vec<f64>,
    draws: Vec<Draw>,
}

#[derive(Deserialize, Clone, Copy)]
struct Draw {
    mass_1: f64,
    mass_2: f64,
    spin_1: [f64; 3],
    spin_2: [f64; 3],
    luminosity_distance: f64,
    inclination: f64,
    phi_ref: f64,
    f_ref: f64,
}

impl Draw {
    fn binary(&self) -> PrecessingBinary {
        PrecessingBinary {
            mass_1: self.mass_1,
            mass_2: self.mass_2,
            spin_1: self.spin_1,
            spin_2: self.spin_2,
            luminosity_distance: self.luminosity_distance,
            inclination: self.inclination,
            phi_ref: self.phi_ref,
            f_ref: self.f_ref,
        }
    }
    fn total_mass(&self) -> f64 {
        self.mass_1 + self.mass_2
    }
    fn mass_ratio(&self) -> f64 {
        self.mass_2 / self.mass_1
    }
}

type Polarizations = (Vec<Complex<f64>>, Vec<Complex<f64>>);

fn generate(name: &str, draw: &Draw, freqs: &[f64]) -> Option<Polarizations> {
    let (m, q) = (draw.total_mass(), draw.mass_ratio());
    let (c1, c2) = (draw.spin_1[2], draw.spin_2[2]);
    let d = draw.luminosity_distance;
    let version = |name: &str| match name {
        n if n.ends_with("(102)") => PrecessionVersion::Nnlo,
        n if n.ends_with("(223)") => PrecessionVersion::Msa,
        _ => PrecessionVersion::SpinTaylor,
    };
    match name {
        "TaylorF2" => Some(
            TaylorF2::new(m, q, c1, c2, d)
                .with_reference_frequency(draw.f_ref)
                .waveform_modes(freqs, draw.phi_ref, draw.inclination),
        ),
        "IMRPhenomD" => Some(
            IMRPhenomD::new(m, q, c1, c2, d)
                .with_reference_frequency(draw.f_ref)
                .waveform_modes(freqs, draw.phi_ref, draw.inclination),
        ),
        "IMRPhenomXAS" => IMRPhenomXAS::from_components(draw.mass_1, draw.mass_2, c1, c2, d)
            .ok()
            .map(|wf| {
                wf.with_reference_frequency(draw.f_ref).waveform_modes(
                    freqs,
                    draw.phi_ref,
                    draw.inclination,
                )
            }),
        "IMRPhenomXHM" => IMRPhenomXHM::from_components(draw.mass_1, draw.mass_2, c1, c2, d)
            .ok()
            .map(|wf| {
                wf.with_reference_frequency(draw.f_ref).polarizations(
                    freqs,
                    draw.phi_ref,
                    draw.inclination,
                )
            }),
        "IMRPhenomXPNR" => xpnr(draw).and_then(|wf| wf.polarizations(freqs)).ok(),
        n if n.starts_with("IMRPhenomXPHM") => {
            IMRPhenomXPHM::with_version(&draw.binary(), version(n))
                .and_then(|wf| wf.polarizations(freqs))
                .ok()
        }
        n if n.starts_with("IMRPhenomXP") => IMRPhenomXP::with_version(&draw.binary(), version(n))
            .and_then(|wf| wf.polarizations(freqs))
            .ok(),
        other => panic!("unknown model {other}"),
    }
}

/// The grid path (`*_on_grid`) with buffers reused across draws; `None` if
/// the model has no grid path.
fn generate_on_grid(
    name: &str,
    draw: &Draw,
    grid: &FrequencyGrid,
    hp: &mut [Complex<f64>],
    hc: &mut [Complex<f64>],
) -> Option<bool> {
    let (m, q) = (draw.total_mass(), draw.mass_ratio());
    let (c1, c2) = (draw.spin_1[2], draw.spin_2[2]);
    let d = draw.luminosity_distance;
    match name {
        "TaylorF2" => {
            TaylorF2::new(m, q, c1, c2, d)
                .with_reference_frequency(draw.f_ref)
                .waveform_modes_on_grid(grid, draw.phi_ref, draw.inclination, hp, hc);
            Some(true)
        }
        "IMRPhenomD" => {
            IMRPhenomD::new(m, q, c1, c2, d)
                .with_reference_frequency(draw.f_ref)
                .waveform_modes_on_grid(grid, draw.phi_ref, draw.inclination, hp, hc);
            Some(true)
        }
        "IMRPhenomXAS" => Some(
            IMRPhenomXAS::from_components(draw.mass_1, draw.mass_2, c1, c2, d)
                .map(|wf| {
                    wf.with_reference_frequency(draw.f_ref)
                        .waveform_modes_on_grid(grid, draw.phi_ref, draw.inclination, hp, hc)
                })
                .is_ok(),
        ),
        "IMRPhenomXHM" => Some(
            IMRPhenomXHM::from_components(draw.mass_1, draw.mass_2, c1, c2, d)
                .map(|wf| {
                    wf.with_reference_frequency(draw.f_ref)
                        .polarizations_on_grid(grid, draw.phi_ref, draw.inclination, hp, hc)
                })
                .is_ok(),
        ),
        "IMRPhenomXPNR" => Some(
            xpnr(draw)
                .and_then(|wf| wf.polarizations_on_grid(grid, hp, hc))
                .is_ok(),
        ),
        n if n.starts_with("IMRPhenomXPHM") => Some(
            IMRPhenomXPHM::with_version(&draw.binary(), precession_version(n))
                .and_then(|wf| wf.polarizations_on_grid(grid, hp, hc))
                .is_ok(),
        ),
        n if n.starts_with("IMRPhenomXP") => Some(
            IMRPhenomXP::with_version(&draw.binary(), precession_version(n))
                .and_then(|wf| wf.polarizations_on_grid(grid, hp, hc))
                .is_ok(),
        ),
        _ => None,
    }
}

/// IMRPhenomXPNR as lalsimulation 6.2.0 computes it (without the parity
/// fix of lalsuite !2593, which costs the same), to compare with LAL.
fn xpnr(draw: &Draw) -> Result<IMRPhenomXPNR, IMRPhenomXError> {
    Ok(IMRPhenomXPNR::new(&draw.binary())?.with_parity_fix(false))
}

fn precession_version(name: &str) -> PrecessionVersion {
    match name {
        n if n.ends_with("(102)") => PrecessionVersion::Nnlo,
        n if n.ends_with("(223)") => PrecessionVersion::Msa,
        _ => PrecessionVersion::SpinTaylor,
    }
}

const MODELS: [&str; 11] = [
    "TaylorF2",
    "IMRPhenomD",
    "IMRPhenomXAS",
    "IMRPhenomXHM",
    "IMRPhenomXP (102)",
    "IMRPhenomXP (223)",
    "IMRPhenomXP (320)",
    "IMRPhenomXPHM (102)",
    "IMRPhenomXPHM (223)",
    "IMRPhenomXPHM (320)",
    "IMRPhenomXPNR",
];

fn main() {
    let path = std::env::args()
        .nth(1)
        .expect("usage: serial_benchmark draws.json");
    let input: Input = serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap();
    let freqs = &input.frequencies;
    let mut entries = Vec::new();
    for name in MODELS {
        // warm up (and fault in the allocator) on the first few draws
        for draw in input.draws.iter().take(5) {
            black_box(generate(name, draw, freqs));
        }
        let mut failures = 0;
        let start = Instant::now();
        for draw in &input.draws {
            if black_box(generate(name, draw, black_box(freqs))).is_none() {
                failures += 1;
            }
        }
        let per_waveform = start.elapsed().as_secs_f64() / input.draws.len() as f64;

        // grid path, buffers reused across draws
        let grid = FrequencyGrid::new(freqs);
        let zero = Complex::new(0.0, 0.0);
        let (mut hp, mut hc) = (vec![zero; grid.len()], vec![zero; grid.len()]);
        let grid_per_waveform = generate_on_grid(name, &input.draws[0], &grid, &mut hp, &mut hc)
            .map(|_| {
                for draw in input.draws.iter().take(5) {
                    generate_on_grid(name, draw, &grid, &mut hp, &mut hc);
                }
                let start = Instant::now();
                for draw in &input.draws {
                    black_box(generate_on_grid(
                        name,
                        draw,
                        black_box(&grid),
                        &mut hp,
                        &mut hc,
                    ));
                }
                start.elapsed().as_secs_f64() / input.draws.len() as f64
            });
        eprintln!(
            "{name:<22} {:>9.1} us per waveform ({failures} errors)",
            per_waveform * 1e6
        );
        if let Some(t) = grid_per_waveform {
            eprintln!(
                "{:<22} {:>9.1} us per waveform on a FrequencyGrid",
                "",
                t * 1e6
            );
        }
        let grid_entry = grid_per_waveform.map_or("null".to_string(), |t| format!("{t:e}"));
        entries.push(format!(
            "\"{name}\": {{\"seconds\": {per_waveform:e}, \"failures\": {failures}, \
             \"grid_seconds\": {grid_entry}}}"
        ));
    }
    if let Ok(data) = NRSur7dq4Data::global() {
        entries.push(nrsur_serial(data, &input.draws));
    }
    println!("{{{}}}", entries.join(", "));
}

/// NRSur7dq4 in the time domain at 4096 Hz from the start of the surrogate
/// (`f_min = f_ref = 0`): a fresh workspace per call, and one reused
/// workspace and output (`grid_seconds`).
fn nrsur_serial(data: &NRSur7dq4Data, draws: &[Draw]) -> String {
    let delta_t = 1.0 / 4096.0;
    let binary = |draw: &Draw| PrecessingBinary {
        f_ref: 0.0,
        ..draw.binary()
    };
    for draw in draws.iter().take(5) {
        let _ = NRSur7dq4::new(data, &binary(draw)).and_then(|wf| wf.polarizations(delta_t, 0.0));
    }
    let mut failures = 0;
    let start = Instant::now();
    for draw in draws {
        let result =
            NRSur7dq4::new(data, &binary(draw)).and_then(|wf| wf.polarizations(delta_t, 0.0));
        if black_box(result).is_err() {
            failures += 1;
        }
    }
    let per_call = start.elapsed().as_secs_f64() / draws.len() as f64;
    let mut ws = NRSurWorkspace::new();
    let mut out = TimeSeriesPolarizations::default();
    let start = Instant::now();
    for draw in draws {
        let _ = NRSur7dq4::new(data, &binary(draw))
            .and_then(|wf| wf.polarizations_with(&mut ws, delta_t, 0.0, &mut out));
        black_box(&out);
    }
    let reused = start.elapsed().as_secs_f64() / draws.len() as f64;
    eprintln!(
        "{:<22} {:>9.1} us per waveform ({failures} errors)\n{:<22} {:>9.1} us with a reused workspace",
        "NRSur7dq4",
        per_call * 1e6,
        "",
        reused * 1e6
    );
    format!(
        "\"NRSur7dq4\": {{\"seconds\": {per_call:e}, \"failures\": {failures}, \
         \"grid_seconds\": {reused:e}}}"
    )
}

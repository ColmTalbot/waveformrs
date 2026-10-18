//! Throughput of generating many waveforms one after another on a fixed
//! frequency grid, as a sampler does.
//!
//! Reads `{"frequencies": [...], "draws": [{"mass_1", "mass_2", "spin_1",
//! "spin_2", "luminosity_distance", "inclination", "phi_ref", "f_ref"}, ...]}`
//! from the JSON file given as the first argument (written by
//! `report/serial_benchmark.py`) and prints the mean time per waveform for
//! each model as JSON.
//!
//!     cargo run --release --example serial_benchmark -- draws.json
use std::hint::black_box;
use std::time::Instant;

use num_complex::Complex;
use serde::Deserialize;

use waveformrs::{
    imrphenomd::IMRPhenomD,
    imrphenomx::{
        IMRPhenomXAS, IMRPhenomXHM, IMRPhenomXP, IMRPhenomXPHM, PrecessingBinary, PrecessionVersion,
    },
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
                .waveform_modes(freqs, 2.0 * draw.phi_ref, draw.inclination),
        ),
        "IMRPhenomXAS" => IMRPhenomXAS::from_components(draw.mass_1, draw.mass_2, c1, c2, d)
            .ok()
            .map(|wf| {
                wf.with_reference_frequency(draw.f_ref).waveform_modes(
                    freqs,
                    2.0 * draw.phi_ref,
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

const MODELS: [&str; 10] = [
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
        eprintln!(
            "{name:<22} {:>9.1} us per waveform ({failures} errors)",
            per_waveform * 1e6
        );
        entries.push(format!(
            "\"{name}\": {{\"seconds\": {per_waveform:e}, \"failures\": {failures}}}"
        ));
    }
    println!("{{{}}}", entries.join(", "));
}

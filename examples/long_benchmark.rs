//! Timings on long frequency grids (up to ~520k frequencies), single
//! thread: the per-frequency methods (construction plus polarizations,
//! output allocated per call), the grid paths (`*_on_grid`, grid prepared
//! once and buffers reused) and, for IMRPhenomXP/XPHM, `polarizations_uniform`.
//!
//! Reads `{"repeats": n, "cases": [{"label", "model", "version", "mass_1",
//! "mass_2", "spin_1", "spin_2", "luminosity_distance", "inclination",
//! "phi_ref", "f_ref", "f_min", "f_max", "delta_f"}, ...]}` from the JSON
//! file given as the first argument (written by `report/long_benchmark.py`)
//! and prints the best of `repeats` timings of each path as JSON (each path
//! is first run once to warm up; with `repeats = 0` nothing is timed, which
//! is used to measure the peak memory of one evaluation of each path). The
//! paths of a case alternate, so that load from other processes affects
//! them alike.
//!
//!     cargo run --release --example long_benchmark -- cases.json
use std::hint::black_box;
use std::time::Instant;

use num_complex::Complex;
use serde::Deserialize;

use waveformrs::{
    grid::FrequencyGrid,
    imrphenomd::IMRPhenomD,
    imrphenomx::{
        IMRPhenomXAS, IMRPhenomXHM, IMRPhenomXP, IMRPhenomXPHM, PrecessingBinary, PrecessionVersion,
    },
    taylor::TaylorF2,
    waveform::Waveform,
};

#[derive(Deserialize)]
struct Input {
    repeats: usize,
    cases: Vec<Case>,
}

#[derive(Deserialize)]
struct Case {
    label: String,
    model: String,
    #[serde(default)]
    version: Option<u32>,
    mass_1: f64,
    mass_2: f64,
    spin_1: [f64; 3],
    spin_2: [f64; 3],
    luminosity_distance: f64,
    inclination: f64,
    phi_ref: f64,
    f_ref: f64,
    f_min: f64,
    f_max: f64,
    delta_f: f64,
}

type C64 = Complex<f64>;

impl Case {
    fn frequencies(&self) -> Vec<f64> {
        let n = ((self.f_max - self.f_min) / self.delta_f).round() as usize + 1;
        (0..n)
            .map(|k| self.f_min + k as f64 * self.delta_f)
            .collect()
    }

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

    fn version(&self) -> PrecessionVersion {
        match self.version {
            Some(102) => PrecessionVersion::Nnlo,
            Some(223) => PrecessionVersion::Msa,
            Some(320) => PrecessionVersion::SpinTaylor,
            other => panic!("{}: unexpected version {other:?}", self.label),
        }
    }

    fn total_mass(&self) -> f64 {
        self.mass_1 + self.mass_2
    }

    fn mass_ratio(&self) -> f64 {
        self.mass_2 / self.mass_1
    }

    /// Construction plus polarizations at `frequencies`.
    fn per_call(&self, frequencies: &[f64]) -> (Vec<C64>, Vec<C64>) {
        let (m, q, d) = (
            self.total_mass(),
            self.mass_ratio(),
            self.luminosity_distance,
        );
        let (c1, c2) = (self.spin_1[2], self.spin_2[2]);
        match self.model.as_str() {
            "TaylorF2" => TaylorF2::new(m, q, c1, c2, d)
                .with_reference_frequency(self.f_ref)
                .waveform_modes(frequencies, self.phi_ref, self.inclination),
            "IMRPhenomD" => IMRPhenomD::new(m, q, c1, c2, d)
                .with_reference_frequency(self.f_ref)
                .waveform_modes(frequencies, self.phi_ref, self.inclination),
            "IMRPhenomXAS" => IMRPhenomXAS::from_components(self.mass_1, self.mass_2, c1, c2, d)
                .unwrap()
                .with_reference_frequency(self.f_ref)
                .waveform_modes(frequencies, self.phi_ref, self.inclination),
            "IMRPhenomXHM" => IMRPhenomXHM::from_components(self.mass_1, self.mass_2, c1, c2, d)
                .unwrap()
                .with_reference_frequency(self.f_ref)
                .polarizations(frequencies, self.phi_ref, self.inclination),
            "IMRPhenomXP" => IMRPhenomXP::with_version(&self.binary(), self.version())
                .and_then(|model| model.polarizations(frequencies))
                .unwrap(),
            "IMRPhenomXPHM" => IMRPhenomXPHM::with_version(&self.binary(), self.version())
                .and_then(|model| model.polarizations(frequencies))
                .unwrap(),
            other => panic!("unknown model {other}"),
        }
    }

    /// Construction plus polarizations on a prepared grid into `hp`, `hc`.
    fn on_grid(&self, grid: &FrequencyGrid, hp: &mut [C64], hc: &mut [C64]) {
        let (m, q, d) = (
            self.total_mass(),
            self.mass_ratio(),
            self.luminosity_distance,
        );
        let (c1, c2) = (self.spin_1[2], self.spin_2[2]);
        match self.model.as_str() {
            "TaylorF2" => TaylorF2::new(m, q, c1, c2, d)
                .with_reference_frequency(self.f_ref)
                .waveform_modes_on_grid(grid, self.phi_ref, self.inclination, hp, hc),
            "IMRPhenomD" => IMRPhenomD::new(m, q, c1, c2, d)
                .with_reference_frequency(self.f_ref)
                .waveform_modes_on_grid(grid, self.phi_ref, self.inclination, hp, hc),
            "IMRPhenomXAS" => IMRPhenomXAS::from_components(self.mass_1, self.mass_2, c1, c2, d)
                .unwrap()
                .with_reference_frequency(self.f_ref)
                .waveform_modes_on_grid(grid, self.phi_ref, self.inclination, hp, hc),
            "IMRPhenomXHM" => IMRPhenomXHM::from_components(self.mass_1, self.mass_2, c1, c2, d)
                .unwrap()
                .with_reference_frequency(self.f_ref)
                .polarizations_on_grid(grid, self.phi_ref, self.inclination, hp, hc),
            "IMRPhenomXP" => IMRPhenomXP::with_version(&self.binary(), self.version())
                .and_then(|model| model.polarizations_on_grid(grid, hp, hc))
                .unwrap(),
            "IMRPhenomXPHM" => IMRPhenomXPHM::with_version(&self.binary(), self.version())
                .and_then(|model| model.polarizations_on_grid(grid, hp, hc))
                .unwrap(),
            other => panic!("unknown model {other}"),
        }
    }

    /// `polarizations_uniform` (as `SimInspiralChooseFDWaveform` without
    /// multibanding), for the models that have it.
    fn uniform(&self) -> Option<(Vec<C64>, Vec<C64>)> {
        let (delta_f, f_min, f_max) = (self.delta_f, self.f_min, self.f_max);
        match self.model.as_str() {
            "IMRPhenomXP" => Some(
                IMRPhenomXP::with_version(&self.binary(), self.version())
                    .and_then(|model| model.polarizations_uniform(delta_f, f_min, f_max))
                    .unwrap(),
            ),
            "IMRPhenomXPHM" => Some(
                IMRPhenomXPHM::with_version(&self.binary(), self.version())
                    .and_then(|model| model.polarizations_uniform(delta_f, f_min, f_max))
                    .unwrap(),
            ),
            _ => None,
        }
    }
}

/// Seconds taken by `f`.
fn time<T>(f: impl FnOnce() -> T) -> f64 {
    let start = Instant::now();
    black_box(f());
    start.elapsed().as_secs_f64()
}

fn main() {
    let path = std::env::args()
        .nth(1)
        .expect("usage: long_benchmark cases.json");
    let input: Input = serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap();
    let mut entries = Vec::new();
    for case in &input.cases {
        let frequencies = case.frequencies();
        let grid_setup = time(|| FrequencyGrid::new(&frequencies));
        let grid = FrequencyGrid::new(&frequencies);
        let zero = C64::new(0.0, 0.0);
        let (mut hp, mut hc) = (vec![zero; grid.len()], vec![zero; grid.len()]);
        // warm up each path once
        black_box(case.per_call(&frequencies));
        case.on_grid(&grid, &mut hp, &mut hc);
        black_box(case.uniform());
        // best of `repeats` (None for no repeats, or no such path)
        let best = |current: Option<f64>, t: f64| Some(current.map_or(t, |c: f64| c.min(t)));
        let (mut per_call, mut on_grid, mut uniform) = (None, None, None);
        for _ in 0..input.repeats {
            per_call = best(per_call, time(|| case.per_call(black_box(&frequencies))));
            on_grid = best(
                on_grid,
                time(|| case.on_grid(black_box(&grid), &mut hp, &mut hc)),
            );
            if case.model.starts_with("IMRPhenomXP") {
                uniform = best(uniform, time(|| case.uniform()));
            }
        }
        let ms = |t: Option<f64>| t.map_or("-".to_string(), |t| format!("{:.2} ms", t * 1e3));
        eprintln!(
            "{:<28} {:>8} frequencies: per call {:>11}, grid {:>11}, uniform {:>11}",
            case.label,
            frequencies.len(),
            ms(per_call),
            ms(on_grid),
            ms(uniform)
        );
        let json = |t: Option<f64>| t.map_or("null".to_string(), |t| format!("{t:e}"));
        entries.push(format!(
            "\"{}\": {{\"n_frequencies\": {}, \"per_call\": {}, \"grid\": {}, \"uniform\": {}, \
             \"grid_setup\": {grid_setup:e}}}",
            case.label,
            frequencies.len(),
            json(per_call),
            json(on_grid),
            json(uniform)
        ));
    }
    println!("{{{}}}", entries.join(", "));
}

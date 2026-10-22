//! Evaluate NRSur7dq4 polarizations for a batch of parameter sets read as
//! JSON from stdin and write them as JSON to stdout.
//!
//! Used by `scripts/validate_nrsur7dq4.py` to compare against
//! `SimInspiralChooseTDWaveform` over random draws. The data file is found
//! through `NRSUR7DQ4_DATA` or `LAL_DATA_PATH`.
//!
//! Input:  {"cases": [{"mass_1": ..., "mass_2": ..., "spin_1": [x, y, z],
//!          "spin_2": [x, y, z], "luminosity_distance": ..., "inclination": ...,
//!          "phi_ref": ..., "f_ref": ..., "f_min": ..., "delta_t": ...}, ...]}
//! Output: {"results": [{"epoch": ..., "hplus": [...], "hcross": [...],
//!          "seconds": ..., "seconds_fresh": ..., "error": null | "..."}, ...]}
//!
//! `seconds` is the median time of repeated calls for the same binary
//! reusing an [`NRSurWorkspace`] and the output buffers, which reuse the
//! cached inertial modes (the cost when only extrinsic parameters change);
//! `seconds_fresh` that of [`NRSur7dq4::new`] plus
//! [`NRSur7dq4::polarizations`].
use std::io::Read;
use std::time::Instant;

use serde::{Deserialize, Serialize};

use waveformrs::nrsur::{
    NRSur7dq4, NRSur7dq4Data, NRSurWorkspace, PrecessingBinary, TimeSeriesPolarizations,
};

#[derive(Deserialize)]
struct Input {
    cases: Vec<Case>,
}

#[derive(Deserialize)]
struct Case {
    mass_1: f64,
    mass_2: f64,
    spin_1: [f64; 3],
    spin_2: [f64; 3],
    luminosity_distance: f64,
    inclination: f64,
    phi_ref: f64,
    f_ref: f64,
    f_min: f64,
    delta_t: f64,
}

#[derive(Serialize)]
struct Output {
    epoch: f64,
    hplus: Vec<f64>,
    hcross: Vec<f64>,
    seconds: f64,
    seconds_fresh: f64,
    error: Option<String>,
}

fn median(mut times: Vec<f64>) -> f64 {
    times.sort_by(f64::total_cmp);
    times[times.len() / 2]
}

fn main() {
    let data = NRSur7dq4Data::global().unwrap();
    let mut buffer = String::new();
    std::io::stdin().read_to_string(&mut buffer).unwrap();
    let input: Input = serde_json::from_str(&buffer).unwrap();
    let mut workspace = NRSurWorkspace::new();
    let mut polarizations = TimeSeriesPolarizations::default();
    let repeats = 5;
    let results: Vec<Output> = input
        .cases
        .iter()
        .map(|case| {
            let binary = PrecessingBinary {
                mass_1: case.mass_1,
                mass_2: case.mass_2,
                spin_1: case.spin_1,
                spin_2: case.spin_2,
                luminosity_distance: case.luminosity_distance,
                inclination: case.inclination,
                phi_ref: case.phi_ref,
                f_ref: case.f_ref,
            };
            let run = |ws: &mut NRSurWorkspace, out: &mut TimeSeriesPolarizations| {
                NRSur7dq4::new(data, &binary)
                    .and_then(|model| model.polarizations_with(ws, case.delta_t, case.f_min, out))
                    .map_err(|err| err.to_string())
            };
            let result = run(&mut workspace, &mut polarizations);
            let mut reused = Vec::new();
            let mut fresh = Vec::new();
            for _ in 0..repeats {
                let start = Instant::now();
                let _ = std::hint::black_box(run(&mut workspace, &mut polarizations));
                reused.push(start.elapsed().as_secs_f64());
                let start = Instant::now();
                let _ = std::hint::black_box(
                    NRSur7dq4::new(data, &binary)
                        .and_then(|model| model.polarizations(case.delta_t, case.f_min)),
                );
                fresh.push(start.elapsed().as_secs_f64());
            }
            let (error, epoch, hplus, hcross) = match result {
                Ok(()) => (
                    None,
                    polarizations.epoch,
                    polarizations.hplus.clone(),
                    polarizations.hcross.clone(),
                ),
                Err(error) => (Some(error), 0.0, Vec::new(), Vec::new()),
            };
            Output {
                epoch,
                hplus,
                hcross,
                seconds: median(reused),
                seconds_fresh: median(fresh),
                error,
            }
        })
        .collect();
    println!(
        "{}",
        serde_json::to_string(&serde_json::json!({ "results": results })).unwrap()
    );
}

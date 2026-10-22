//! Timing of NRSur7dq4: loading the data, the polarizations and modes for
//! several total masses and sampling rates (q = 2, generic spins), and
//! sampler-like sequences of waveforms.
//!
//! Run with the data file in `LAL_DATA_PATH` (or `NRSUR7DQ4_DATA`):
//!
//! ```sh
//! LAL_DATA_PATH=/path/to/lal_data cargo run --release --example nrsur7dq4
//! ```
//!
//! `scripts/time_nrsur7dq4_lal.py` times `SimInspiralChooseTDWaveform` and
//! `SimInspiralChooseTDModes` (single thread) for the same configurations.
use std::time::Instant;

use waveformrs::nrsur::{
    NRSur7dq4, NRSur7dq4Data, NRSurWorkspace, PrecessingBinary, TimeSeriesModes,
    TimeSeriesPolarizations,
};

/// (total mass in Msun, sampling rate in Hz, f_min in Hz)
pub const CONFIGURATIONS: [(f64, f64, f64); 8] = [
    (60.0, 4096.0, 20.0),
    (60.0, 16384.0, 20.0),
    (20.0, 4096.0, 55.0),
    (20.0, 16384.0, 55.0),
    (100.0, 4096.0, 20.0),
    (100.0, 1024.0, 20.0),
    (200.0, 2048.0, 20.0),
    (300.0, 1024.0, 20.0),
];

fn binary(total_mass: f64) -> PrecessingBinary {
    PrecessingBinary {
        mass_1: total_mass * 2.0 / 3.0,
        mass_2: total_mass / 3.0,
        spin_1: [0.4, -0.3, 0.2],
        spin_2: [0.1, 0.5, -0.3],
        luminosity_distance: 100.0,
        inclination: 0.7,
        phi_ref: 0.1,
        f_ref: 0.0,
    }
}

fn time<F: FnMut()>(iterations: u32, mut f: F) -> std::time::Duration {
    f();
    let start = Instant::now();
    for _ in 0..iterations {
        f();
    }
    start.elapsed() / iterations
}

/// A small deterministic generator for the sampler-like sequences.
struct Rng(u64);

impl Rng {
    fn uniform(&mut self, lo: f64, hi: f64) -> f64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        lo + (hi - lo) * (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }

    fn spin(&mut self) -> [f64; 3] {
        let magnitude = self.uniform(0.0, 0.99);
        let cos_theta = self.uniform(-1.0, 1.0);
        let phi = self.uniform(0.0, 2.0 * std::f64::consts::PI);
        let sin_theta = (1.0 - cos_theta * cos_theta).sqrt();
        [
            magnitude * sin_theta * phi.cos(),
            magnitude * sin_theta * phi.sin(),
            magnitude * cos_theta,
        ]
    }
}

fn main() {
    let path = waveformrs::nrsur::find_data_file()
        .expect("set NRSUR7DQ4_DATA or LAL_DATA_PATH to find NRSur7dq4_v1.0.h5");
    let json = std::env::args().any(|arg| arg == "--json");
    let mut rows = Vec::new();
    let start = Instant::now();
    let data = NRSur7dq4Data::load(&path).unwrap();
    let first = start.elapsed();
    let reload = time(20, || {
        std::hint::black_box(NRSur7dq4Data::load(&path).unwrap());
    });
    eprintln!(
        "load {} ({} fit coefficients): {first:.1?} (first), {reload:.1?} (cached file)",
        path.display(),
        data.n_coefficients(),
    );

    let mut workspace = NRSurWorkspace::new();
    let mut polarizations = TimeSeriesPolarizations::default();
    let mut modes = TimeSeriesModes::default();
    eprintln!(
        "{:>6} {:>6} {:>6} {:>7} | {:>10} {:>10} {:>10} {:>10}",
        "M/Msun", "1/dt", "f_min", "samples", "pols", "pols (ws)", "extrinsic", "modes (ws)"
    );
    for (total_mass, sampling_rate, f_min) in CONFIGURATIONS {
        let binary = binary(total_mass);
        let delta_t = 1.0 / sampling_rate;
        let model = NRSur7dq4::new(&data, &binary).unwrap();
        let iterations = 200;
        // a new model and workspace for every waveform
        let fresh = time(iterations, || {
            let model = NRSur7dq4::new(&data, &binary).unwrap();
            std::hint::black_box(model.polarizations(delta_t, f_min).unwrap());
        });
        // reusing the workspace and output buffers, alternating between two
        // binaries so that nothing is cached
        let other = NRSur7dq4::new(
            &data,
            &PrecessingBinary {
                spin_1: [0.3, 0.1, -0.2],
                ..binary
            },
        )
        .unwrap();
        let mut flip = false;
        let reused = time(iterations, || {
            flip = !flip;
            let model = if flip { &model } else { &other };
            model
                .polarizations_with(&mut workspace, delta_t, f_min, &mut polarizations)
                .unwrap();
            std::hint::black_box(&polarizations);
        });
        // only the extrinsic parameters change: the inertial modes are reused
        let mut inclination = 0.0;
        let extrinsic = time(iterations, || {
            inclination += 0.01;
            let model = NRSur7dq4::new(
                &data,
                &PrecessingBinary {
                    inclination,
                    ..binary
                },
            )
            .unwrap();
            model
                .polarizations_with(&mut workspace, delta_t, f_min, &mut polarizations)
                .unwrap();
            std::hint::black_box(&polarizations);
        });
        let mut flip = false;
        let modes_time = time(iterations / 2, || {
            flip = !flip;
            let model = if flip { &model } else { &other };
            model
                .modes_with(&mut workspace, delta_t, f_min, &mut modes)
                .unwrap();
            std::hint::black_box(&modes);
        });
        eprintln!(
            "{total_mass:>6} {sampling_rate:>6} {f_min:>6} {:>7} | {:>10.1?} {:>10.1?} {:>10.1?} {:>10.1?}",
            polarizations.hplus.len(),
            fresh,
            reused,
            extrinsic,
            modes_time
        );
        rows.push(format!(
            "{{\"total_mass\": {total_mass}, \"sampling_rate\": {sampling_rate}, \
             \"f_min\": {f_min}, \"samples\": {}, \"fresh\": {:e}, \"reused\": {:e}, \
             \"extrinsic\": {:e}, \"modes\": {:e}}}",
            polarizations.hplus.len(),
            fresh.as_secs_f64(),
            reused.as_secs_f64(),
            extrinsic.as_secs_f64(),
            modes_time.as_secs_f64()
        ));
    }

    // A sampler-like sequence: new intrinsic and extrinsic parameters for
    // every waveform (q in [1, 6], |chi| <= 0.99, M in [40, 120] Msun, from
    // 20 Hz at 4096 Hz), reusing one workspace.
    let mut rng = Rng(0x9E3779B97F4A7C15);
    let n = 2000;
    let binaries: Vec<PrecessingBinary> = (0..n)
        .map(|_| {
            let total_mass = rng.uniform(40.0, 120.0);
            let q = rng.uniform(1.0, 6.0);
            PrecessingBinary {
                mass_1: total_mass * q / (1.0 + q),
                mass_2: total_mass / (1.0 + q),
                spin_1: rng.spin(),
                spin_2: rng.spin(),
                luminosity_distance: rng.uniform(100.0, 2000.0),
                inclination: rng.uniform(0.0, std::f64::consts::PI),
                phi_ref: rng.uniform(0.0, 2.0 * std::f64::consts::PI),
                f_ref: 20.0,
            }
        })
        .collect();
    let start = Instant::now();
    let mut samples = 0;
    let mut rejected = 0;
    for binary in &binaries {
        let model = NRSur7dq4::new(&data, binary).unwrap();
        // binaries whose surrogate starts above 20 Hz are rejected, as by LAL
        match model.polarizations_with(&mut workspace, 1.0 / 4096.0, 20.0, &mut polarizations) {
            Ok(()) => samples += polarizations.hplus.len(),
            Err(_) => rejected += 1,
        }
    }
    let per_waveform = start.elapsed() / n as u32;
    eprintln!(
        "sampler-like sequence of {n} binaries (M = 40-120 Msun from 20 Hz, 4096 Hz): \
         {per_waveform:.1?} per waveform ({} samples on average, {rejected} rejected)",
        samples / (n - rejected).max(1)
    );
    if json {
        println!(
            "{{\"load_first\": {:e}, \"load_cached_file\": {:e}, \"sampler\": {:e}, \
             \"configurations\": [{}]}}",
            first.as_secs_f64(),
            reload.as_secs_f64(),
            per_waveform.as_secs_f64(),
            rows.join(", ")
        );
    }
}

//! Timing of the conversions between the time and frequency domains:
//! `SimInspiralTD` for the frequency-domain models
//! ([`FdWaveform::time_domain`]), `SimInspiralFD` for NRSur7dq4
//! ([`nrsur7dq4_fd`]) and the real FFTs, for the configurations of
//! `scripts/time_conversions_lal.py`.
//!
//! ```sh
//! OMP_NUM_THREADS=1 LAL_DATA_PATH=/path/to/lal_data \
//!     python scripts/time_conversions_lal.py --output lal_timings.json
//! LAL_DATA_PATH=/path/to/lal_data \
//!     cargo run --release --example conversions -- lal_timings.json [--json out.json]
//! ```
//!
//! Without the LAL timings only the Rust times are printed; without the
//! NRSur7dq4 data the NRSur7dq4 rows are skipped.
use std::collections::HashMap;
use std::time::{Duration, Instant};

use num_complex::Complex;
use waveformrs::conversions::{nrsur7dq4_fd, FdWaveform};
use waveformrs::fft::RealFft;
use waveformrs::nrsur::{NRSur7dq4, NRSur7dq4Data};
use waveformrs::{Approximant, Binary, WaveformOptions};

const ALIGNED: ([f64; 3], [f64; 3]) = ([0.0, 0.0, 0.3], [0.0, 0.0, -0.4]);
const PRECESSING: ([f64; 3], [f64; 3]) = ([0.4, 0.2, 0.3], [-0.2, 0.3, -0.1]);

/// Mean time per call over at least `minimum`, after one warm-up call.
fn per_call<F: FnMut()>(minimum: Duration, mut f: F) -> Duration {
    f();
    let start = Instant::now();
    let mut n = 0u32;
    while start.elapsed() < minimum {
        f();
        n += 1;
    }
    start.elapsed() / n
}

fn binary(m1: f64, m2: f64, spins: ([f64; 3], [f64; 3])) -> Binary {
    Binary {
        mass_1: m1,
        mass_2: m2,
        spin_1: spins.0,
        spin_2: spins.1,
        luminosity_distance: 500.0,
        inclination: 0.7,
        phi_ref: 1.2,
        f_ref: 0.0,
    }
}

struct Table {
    lal: HashMap<String, f64>,
    rows: std::cell::RefCell<Vec<String>>,
}

impl Table {
    fn row(&self, label: &str, what: &str, rust: Duration, samples: usize) {
        let rust_ms = rust.as_secs_f64() * 1e3;
        let lal_seconds = self
            .lal
            .get(label)
            .map_or("null".to_string(), |t| format!("{t:e}"));
        self.rows.borrow_mut().push(format!(
            "{{\"label\": {label:?}, \"what\": {what:?}, \"rust\": {:e}, \"lal\": {lal_seconds}, \
             \"samples\": {samples}}}",
            rust.as_secs_f64()
        ));
        match self.lal.get(label) {
            Some(&lal) => println!(
                "{label:<42} {what:<14} {rust_ms:9.3} ms  LAL {:9.3} ms  x{:6.2}  ({samples} samples)",
                lal * 1e3,
                lal * 1e3 / rust_ms
            ),
            None => println!("{label:<42} {what:<14} {rust_ms:9.3} ms  ({samples} samples)"),
        }
    }
}

fn main() {
    // arguments: [lal_timings.json] [--json output.json]
    let args: Vec<String> = std::env::args().skip(1).collect();
    let json_output = args
        .iter()
        .position(|a| a == "--json")
        .and_then(|i| args.get(i + 1).cloned());
    let lal = args
        .iter()
        .take_while(|a| *a != "--json")
        .next()
        .cloned()
        .map(|path| {
            let contents = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{path}: {e}"));
            let value: serde_json::Value = serde_json::from_str(&contents).unwrap();
            value["seconds"]
                .as_object()
                .unwrap()
                .iter()
                .map(|(k, v)| (k.clone(), v.as_f64().unwrap()))
                .collect()
        })
        .unwrap_or_default();
    let table = Table {
        lal,
        rows: Default::default(),
    };
    let minimum = Duration::from_secs(1);

    let td_configurations = [
        (
            "IMRPhenomD 36+29 Msun",
            Approximant::IMRPhenomD,
            (36.0, 29.0),
            ALIGNED,
            20.0,
        ),
        (
            "IMRPhenomXAS 36+29 Msun",
            Approximant::IMRPhenomXAS,
            (36.0, 29.0),
            ALIGNED,
            20.0,
        ),
        (
            "IMRPhenomXHM 36+29 Msun",
            Approximant::IMRPhenomXHM,
            (36.0, 29.0),
            ALIGNED,
            20.0,
        ),
        (
            "IMRPhenomXP 36+29 Msun",
            Approximant::IMRPhenomXP,
            (36.0, 29.0),
            PRECESSING,
            20.0,
        ),
        (
            "IMRPhenomXPHM 36+29 Msun",
            Approximant::IMRPhenomXPHM,
            (36.0, 29.0),
            PRECESSING,
            20.0,
        ),
        (
            "IMRPhenomD 10+5 Msun, 10 Hz",
            Approximant::IMRPhenomD,
            (10.0, 5.0),
            ALIGNED,
            10.0,
        ),
        (
            "TaylorF2 1.4+1.3 Msun, 30 Hz",
            Approximant::TaylorF2,
            (1.4, 1.3),
            ([0.0; 3], [0.0; 3]),
            30.0,
        ),
    ];
    for (label, approximant, (m1, m2), spins, f_min) in td_configurations {
        // PhenomXPrecVersion 320 (the default here), no multibanding
        let model = FdWaveform::new(
            approximant,
            binary(m1, m2, spins),
            WaveformOptions::default(),
        )
        .unwrap();
        let mut samples = 0;
        let time = per_call(minimum, || {
            let h = model.time_domain(1.0 / 4096.0, f_min).unwrap();
            samples = h.hplus.data.len();
        });
        table.row(label, "SimInspiralTD", time, samples);
    }

    match NRSur7dq4Data::global() {
        Ok(data) => {
            for (label, delta_f) in [
                ("NRSur7dq4 40+20 Msun, df = 0, 2048 Hz", 0.0),
                ("NRSur7dq4 40+20 Msun, df = 1/8, 2048 Hz", 0.125),
            ] {
                let binary = Binary {
                    mass_1: 40.0,
                    mass_2: 20.0,
                    spin_1: [0.4, -0.3, 0.2],
                    spin_2: [0.1, 0.5, -0.3],
                    luminosity_distance: 100.0,
                    inclination: 0.7,
                    phi_ref: 0.1,
                    f_ref: 20.0,
                };
                let mut samples = 0;
                let time = per_call(minimum, || {
                    let model = NRSur7dq4::new(data, &binary).unwrap();
                    let h = nrsur7dq4_fd(&model, delta_f, 20.0, 2048.0).unwrap();
                    samples = h.hplus.data.len();
                });
                table.row(label, "SimInspiralFD", time, samples);
            }
        }
        Err(err) => println!("skipping NRSur7dq4: {err}"),
    }

    for log2n in [16u32, 20] {
        let n = 1usize << log2n;
        let input: Vec<f64> = (0..n).map(|j| ((j * j) as f64 * 1e-3).sin()).collect();
        let mut spectrum = vec![Complex::new(0.0, 0.0); n / 2 + 1];
        let mut output = vec![0.0; n];
        let forward = per_call(minimum, || {
            RealFft::new(n).forward(&input, &mut spectrum);
        });
        table.row(
            &format!("forward real FFT, N = 2^{log2n}"),
            "(plan + FFT)",
            forward,
            n,
        );
        let reverse = per_call(minimum, || {
            RealFft::new(n).inverse(&spectrum, &mut output);
        });
        table.row(
            &format!("reverse real FFT, N = 2^{log2n}"),
            "(plan + FFT)",
            reverse,
            n,
        );
        let plan = RealFft::new(n);
        let cached = per_call(minimum, || plan.forward(&input, &mut spectrum));
        table.row(
            &format!("forward real FFT, N = 2^{log2n}"),
            "(cached plan)",
            cached,
            n,
        );
    }
    if let Some(path) = json_output {
        let rows = table.rows.borrow().join(",\n  ");
        std::fs::write(&path, format!("{{\"rows\": [\n  {rows}\n]}}\n"))
            .unwrap_or_else(|e| panic!("{path}: {e}"));
    }
}

//! Comparison of the SpinTaylorT4 orbit evolution against LALSimulation.
//!
//! Reference data are produced by `scripts/generate_spintaylor_reference.py`
//! with `SimInspiralSpinTaylorPNEvolveOrbit` and stored in
//! `tests/data/spintaylor.json`.
use serde::Deserialize;
use waveformrs::spintaylor::{SpinTaylorError, SpinTaylorOrbit, SpinTaylorT4, Termination};

/// The reference data were generated with the osx-arm64 build of
/// lalsimulation, which this crate reproduces bit for bit on the same
/// platform.  Elsewhere LAL itself is compiled differently (no fused
/// multiply-adds on x86-64, a different libm), and the adaptive step size
/// control amplifies last-bit differences.
const EXACT: bool = cfg!(all(target_arch = "aarch64", target_os = "macos"));

#[derive(Debug, Deserialize)]
struct ReferenceFile {
    lalsimulation_version: String,
    cases: Vec<Case>,
}

#[derive(Debug, Deserialize)]
struct Case {
    label: String,
    m1_si: f64,
    m2_si: f64,
    spin_1: [f64; 3],
    spin_2: [f64; 3],
    lnhat: [f64; 3],
    e1: [f64; 3],
    delta_t: f64,
    f_start: f64,
    f_end: f64,
    lambda: [f64; 2],
    quadparam: [f64; 2],
    phase_order: i32,
    spin_order: i32,
    tidal_order: i32,
    lscorr: bool,
    length: usize,
    indices: Vec<usize>,
    v: Vec<f64>,
    phi: Vec<f64>,
    s1_x: Vec<f64>,
    s1_y: Vec<f64>,
    s1_z: Vec<f64>,
    s2_x: Vec<f64>,
    s2_y: Vec<f64>,
    s2_z: Vec<f64>,
    lnhat_x: Vec<f64>,
    lnhat_y: Vec<f64>,
    lnhat_z: Vec<f64>,
    e1_x: Vec<f64>,
    e1_y: Vec<f64>,
    e1_z: Vec<f64>,
}

impl Case {
    fn binary(&self) -> SpinTaylorT4 {
        SpinTaylorT4::new(self.m1_si, self.m2_si, self.spin_1, self.spin_2)
            .with_orientation(self.lnhat, self.e1)
            .with_matter(self.lambda, self.quadparam)
            .with_pn_orders(self.phase_order, self.spin_order, self.tidal_order)
            .with_lscorr(self.lscorr)
    }

    fn reference(&self) -> [(&'static str, &[f64]); 14] {
        [
            ("v", &self.v),
            ("phi", &self.phi),
            ("s1_x", &self.s1_x),
            ("s1_y", &self.s1_y),
            ("s1_z", &self.s1_z),
            ("s2_x", &self.s2_x),
            ("s2_y", &self.s2_y),
            ("s2_z", &self.s2_z),
            ("lnhat_x", &self.lnhat_x),
            ("lnhat_y", &self.lnhat_y),
            ("lnhat_z", &self.lnhat_z),
            ("e1_x", &self.e1_x),
            ("e1_y", &self.e1_y),
            ("e1_z", &self.e1_z),
        ]
    }
}

fn series(orbit: &SpinTaylorOrbit) -> [&[f64]; 14] {
    [
        &orbit.v,
        &orbit.phi,
        &orbit.s1_x,
        &orbit.s1_y,
        &orbit.s1_z,
        &orbit.s2_x,
        &orbit.s2_y,
        &orbit.s2_z,
        &orbit.lnhat_x,
        &orbit.lnhat_y,
        &orbit.lnhat_z,
        &orbit.e1_x,
        &orbit.e1_y,
        &orbit.e1_z,
    ]
}

fn load() -> ReferenceFile {
    let path = format!("{}/tests/data/spintaylor.json", env!("CARGO_MANIFEST_DIR"));
    let contents = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{path}: {e}"));
    serde_json::from_str(&contents).unwrap()
}

#[test]
fn spintaylor_t4_matches_lalsimulation() {
    let reference = load();
    // max |rust - lal| / max |lal| for each series
    let tolerance = if EXACT { 1e-12 } else { 1e-6 };
    let mut failures = Vec::new();
    for case in &reference.cases {
        let orbit = case
            .binary()
            .evolve_orbit(case.delta_t, case.f_start, case.f_end)
            .unwrap();
        if orbit.len() != case.length {
            failures.push(format!(
                "{}: length {} != {}",
                case.label,
                orbit.len(),
                case.length
            ));
            continue;
        }
        let mut worst: f64 = 0.0;
        for ((name, truth), test) in case.reference().iter().zip(series(&orbit)) {
            assert_eq!(test.len(), case.length);
            let scale = truth.iter().fold(0.0f64, |acc, x| acc.max(x.abs()));
            let max_diff = case
                .indices
                .iter()
                .zip(truth.iter())
                .fold(0.0f64, |acc, (&i, t)| acc.max((test[i] - t).abs()));
            let error = if scale > 0.0 {
                max_diff / scale
            } else {
                max_diff
            };
            worst = worst.max(error);
            if error.is_nan() || error > tolerance {
                failures.push(format!("{}: {name} relative error {error:e}", case.label));
            }
        }
        println!(
            "{:45} {:7} samples, {:?}, max relative error {worst:.3e}",
            case.label,
            orbit.len(),
            orbit.termination()
        );
    }
    assert!(
        failures.is_empty(),
        "comparison with lalsimulation {} failed:\n{}",
        reference.lalsimulation_version,
        failures.join("\n")
    );
}

#[test]
fn backward_integration_ends_at_start_frequency() {
    let binary = SpinTaylorT4::from_solar_masses(30.0, 20.0, [0.5, 0.3, 0.2], [-0.3, 0.4, 0.1]);
    let delta_t = 1.0 / 512.0;
    let backward = binary.evolve_orbit(delta_t, 25.0, 14.5).unwrap();
    let forward = binary.evolve_orbit(delta_t, 25.0, 400.0).unwrap();
    assert_eq!(backward.termination(), Termination::FrequencyBound);
    assert!(forward.termination().is_normal());
    // Both start from the same state at f_start: the last backward sample
    // and the first forward sample.
    let last = backward.len() - 1;
    assert_eq!(backward.v[last], forward.v[0]);
    assert_eq!(backward.phi[last], 0.0);
    assert_eq!(forward.phi[0], 0.0);
    assert_eq!(backward.s1_x[last], forward.s1_x[0]);
    assert_eq!(backward.lnhat_z[last], 1.0);
    // ordered forward in time: v increases, the first sample is below f_end
    assert!(backward.v.windows(2).all(|w| w[0] < w[1]));
    assert!(forward.v.windows(2).all(|w| w[0] < w[1]));
    let mtotal = 50.0 * 4.925490947641267e-6;
    let f = |v: f64| v * v * v / (std::f64::consts::PI * mtotal);
    assert!(f(backward.v[0]) < 14.5 && f(backward.v[1]) >= 14.5);
    assert_eq!(backward.epoch, -((backward.len() - 1) as f64) * delta_t);
    // spins keep their magnitude, LNhat and E1 stay orthonormal (up to the
    // integration and interpolation errors)
    let mut drift: f64 = 0.0;
    for orbit in [&backward, &forward] {
        for i in 0..orbit.len() {
            let s1 = orbit.s1_x[i].powi(2) + orbit.s1_y[i].powi(2) + orbit.s1_z[i].powi(2);
            let ln = [orbit.lnhat_x[i], orbit.lnhat_y[i], orbit.lnhat_z[i]];
            let e1 = [orbit.e1_x[i], orbit.e1_y[i], orbit.e1_z[i]];
            drift = drift
                .max((s1 - 0.38).abs())
                .max((ln.iter().map(|x| x * x).sum::<f64>() - 1.0).abs())
                .max(ln.iter().zip(&e1).map(|(a, b)| a * b).sum::<f64>().abs());
        }
    }
    println!("max drift of conserved quantities {drift:e}");
    assert!(drift < 1e-7, "{drift:e}");
}

#[test]
fn invalid_inputs_are_rejected() {
    let binary = SpinTaylorT4::from_solar_masses(30.0, 20.0, [0.1, 0.0, 0.0], [0.0, 0.1, 0.0]);
    assert_eq!(
        binary.evolve_orbit(1e-3, 0.0, 100.0).unwrap_err(),
        SpinTaylorError::InvalidStartFrequency(0.0)
    );
    assert_eq!(
        binary.evolve_orbit(1e-3, 20.0, -1.0).unwrap_err(),
        SpinTaylorError::InvalidEndFrequency(-1.0)
    );
    assert_eq!(
        binary
            .with_pn_orders(7, 7, 12)
            .evolve_orbit(1e-3, 20.0, 100.0)
            .unwrap_err(),
        SpinTaylorError::InvalidSpinOrder(7)
    );
}

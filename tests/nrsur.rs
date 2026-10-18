//! Tests of the HDF5 reader and of NRSur7dq4 against LALSimulation.
//!
//! These need the surrogate data file `NRSur7dq4_v1.0.h5`, looked up from
//! `NRSUR7DQ4_DATA` (the file), then `LAL_DATA_PATH` (directories), then
//! `<CARGO_MANIFEST_DIR>/lal_data/`. Without it the tests are skipped with
//! a message.
use std::collections::HashMap;
use std::path::PathBuf;

use serde::Deserialize;
use waveformrs::hdf5::{self, Datatype, ObjectKind};
use waveformrs::nrsur::{
    NRSur7dq4, NRSur7dq4Data, NRSurWorkspace, PrecessingBinary, TimeSeriesModes,
    TimeSeriesPolarizations,
};

fn data_file() -> Option<PathBuf> {
    let found = waveformrs::nrsur::find_data_file().or_else(|| {
        let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("lal_data")
            .join(waveformrs::nrsur::NRSUR7DQ4_DATAFILE);
        path.is_file().then_some(path)
    });
    if found.is_none() {
        println!(
            "skipping: {} not found (set NRSUR7DQ4_DATA or LAL_DATA_PATH)",
            waveformrs::nrsur::NRSUR7DQ4_DATAFILE
        );
    }
    found
}

#[derive(Deserialize)]
struct ChecksumFile {
    md5: String,
    datasets: HashMap<String, DatasetSummary>,
}

#[derive(Deserialize)]
struct DatasetSummary {
    shape: Vec<usize>,
    dtype: String,
    checksum: String,
}

/// The checksum of `scripts/hdf5_checksums.py`.
fn checksum(bits: impl Iterator<Item = u64>) -> u64 {
    const K1: u64 = 0x9E37_79B9_7F4A_7C15;
    const K2: u64 = 0xBF58_476D_1CE4_E5B9;
    bits.enumerate().fold(0u64, |sum, (i, b)| {
        sum.wrapping_add((b ^ (i as u64).wrapping_mul(K1)).wrapping_mul(K2))
    })
}

fn walk<'f>(group: &hdf5::Group<'f>, prefix: &str, out: &mut Vec<(String, hdf5::Dataset<'f>)>) {
    for name in group.member_names() {
        let path = if prefix.is_empty() {
            name.to_string()
        } else {
            format!("{prefix}/{name}")
        };
        match group.kind(name).unwrap() {
            ObjectKind::Group => walk(&group.group(name).unwrap(), &path, out),
            ObjectKind::Dataset => out.push((path, group.dataset(name).unwrap())),
        }
    }
}

#[test]
fn hdf5_reader_matches_h5py() {
    let Some(path) = data_file() else { return };
    let reference: ChecksumFile = serde_json::from_str(
        &std::fs::read_to_string(concat!(
            env!("CARGO_MANIFEST_DIR"),
            "/tests/data/nrsur7dq4_hdf5_checksums.json"
        ))
        .unwrap(),
    )
    .unwrap();
    let file = hdf5::File::open(&path).unwrap();
    let mut datasets = Vec::new();
    walk(&file.root().unwrap(), "", &mut datasets);
    assert_eq!(datasets.len(), reference.datasets.len());
    for (name, dataset) in &datasets {
        let expected = reference
            .datasets
            .get(name)
            .unwrap_or_else(|| panic!("{name} not in the h5py dump"));
        assert_eq!(dataset.shape(), expected.shape.as_slice(), "{name}");
        let bits: Vec<u64> = match (dataset.dtype(), expected.dtype.as_str()) {
            (
                Datatype::Float {
                    size: 8,
                    big_endian: false,
                },
                "<f8",
            ) => dataset
                .read_f64()
                .unwrap()
                .iter()
                .map(|v| v.to_bits())
                .collect(),
            (
                Datatype::Integer {
                    size: 8,
                    signed: true,
                    big_endian: false,
                },
                "<i8",
            ) => dataset
                .read_i64()
                .unwrap()
                .iter()
                .map(|&v| v as u64)
                .collect(),
            (dtype, expected) => panic!("{name}: {dtype:?} vs {expected}"),
        };
        assert_eq!(
            format!("{:016x}", checksum(bits.into_iter())),
            expected.checksum,
            "{name}"
        );
    }
    println!(
        "{} datasets match h5py (file md5 {})",
        datasets.len(),
        reference.md5
    );
    // path lookups and errors
    let t_ds = file.dataset("/t_ds").unwrap().read_f64().unwrap();
    assert_eq!(t_ds.len(), 233);
    let coefs = file.dataset("ds_node_0/omega_coefs").unwrap();
    assert!(coefs.read_i64().is_err());
    assert!(matches!(
        file.dataset("ds_node_0/missing"),
        Err(hdf5::Hdf5Error::NotFound(_))
    ));
    assert!(matches!(
        file.dataset("ds_node_0"),
        Err(hdf5::Hdf5Error::Type(_))
    ));
}

fn load_data() -> Option<&'static NRSur7dq4Data> {
    static DATA: std::sync::OnceLock<Option<NRSur7dq4Data>> = std::sync::OnceLock::new();
    DATA.get_or_init(|| data_file().map(|path| NRSur7dq4Data::load(path).unwrap()))
        .as_ref()
}

fn reference<T: serde::de::DeserializeOwned>(name: &str) -> T {
    let path = format!("{}/tests/data/{name}.json", env!("CARGO_MANIFEST_DIR"));
    serde_json::from_str(&std::fs::read_to_string(&path).unwrap()).unwrap()
}

#[derive(Deserialize)]
struct ReferenceFile<T> {
    lalsimulation_version: String,
    cases: Vec<T>,
}

#[derive(Deserialize, Debug)]
struct TdCase {
    mass_1: f64,
    mass_2: f64,
    spin_1: [f64; 3],
    spin_2: [f64; 3],
    luminosity_distance: f64,
    #[serde(default)]
    inclination: f64,
    #[serde(default)]
    phi_ref: f64,
    f_min: f64,
    f_ref: f64,
    delta_t: f64,
    epoch: Option<f64>,
    #[serde(default)]
    hplus: Vec<f64>,
    #[serde(default)]
    hcross: Vec<f64>,
    lal_error: Option<String>,
    #[serde(default)]
    modes: Vec<ModeData>,
    length: Option<usize>,
    sample_stride: Option<usize>,
    #[serde(default)]
    mode_array: Vec<(u32, i32)>,
    #[serde(default)]
    unlimited_extrapolation: bool,
}

#[derive(Deserialize, Debug)]
struct ModeData {
    ell: u32,
    emm: i32,
    real: Vec<f64>,
    imag: Vec<f64>,
}

impl TdCase {
    fn model<'a>(&self, data: &'a NRSur7dq4Data) -> NRSur7dq4<'a> {
        let mut model = NRSur7dq4::new(data, &self.binary())
            .unwrap()
            .with_unlimited_extrapolation(self.unlimited_extrapolation);
        if !self.mode_array.is_empty() {
            model = model.with_mode_array(&self.mode_array).unwrap();
        }
        model
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

    fn label(&self) -> String {
        format!(
            "m1={} m2={} s1={:?} s2={:?} f_min={} f_ref={} dt=1/{}{}{}",
            self.mass_1,
            self.mass_2,
            self.spin_1,
            self.spin_2,
            self.f_min,
            self.f_ref,
            (1.0 / self.delta_t).round(),
            if self.mode_array.is_empty() {
                String::new()
            } else {
                format!(" modes={:?}", self.mode_array)
            },
            if self.unlimited_extrapolation {
                " unlimited"
            } else {
                ""
            }
        )
    }
}

/// `(max |a - b| / max |b|, fraction of bit-identical samples)`.
fn compare(test: &[f64], reference: &[f64]) -> (f64, f64) {
    assert_eq!(test.len(), reference.len());
    let scale = reference.iter().fold(0.0f64, |m, v| m.max(v.abs()));
    let diff = test
        .iter()
        .zip(reference)
        .fold(0.0f64, |m, (a, b)| m.max((a - b).abs()));
    let identical = test
        .iter()
        .zip(reference)
        .filter(|(a, b)| a.to_bits() == b.to_bits())
        .count();
    (diff / scale, identical as f64 / test.len().max(1) as f64)
}

#[test]
fn nrsur7dq4_polarizations_match_lalsimulation() {
    let Some(data) = load_data() else { return };
    let reference: ReferenceFile<TdCase> = reference("nrsur7dq4");
    let mut worst: f64 = 0.0;
    let mut workspace = NRSurWorkspace::new();
    let mut out = TimeSeriesPolarizations::default();
    for case in &reference.cases {
        let model = case.model(data);
        let result = model.polarizations_with(&mut workspace, case.delta_t, case.f_min, &mut out);
        if let Some(error) = &case.lal_error {
            let err = result.expect_err(&format!("LAL fails ({error}): {}", case.label()));
            println!("error as LAL ({error}): {err} | {}", case.label());
            continue;
        }
        result.unwrap_or_else(|e| panic!("{e}: {}", case.label()));
        assert_eq!(out.hplus.len(), case.hplus.len(), "{}", case.label());
        assert_eq!(out.epoch, case.epoch.unwrap(), "{}", case.label());
        let (ep, ip) = compare(&out.hplus, &case.hplus);
        let (ec, ic) = compare(&out.hcross, &case.hcross);
        println!(
            "max rel err {ep:.2e} / {ec:.2e}, bit-identical {:.1}% / {:.1}% | {}",
            100.0 * ip,
            100.0 * ic,
            case.label()
        );
        worst = worst.max(ep).max(ec);
    }
    println!(
        "worst {worst:.3e} (lalsimulation {})",
        reference.lalsimulation_version
    );
    assert!(worst < 1e-12, "worst relative error {worst:e}");
}

#[test]
fn nrsur7dq4_modes_match_lalsimulation() {
    let Some(data) = load_data() else { return };
    let reference: ReferenceFile<TdCase> = reference("nrsur7dq4_modes");
    let mut worst: f64 = 0.0;
    for case in &reference.cases {
        let model = case.model(data);
        let modes = model.modes(case.delta_t, case.f_min).unwrap();
        assert_eq!(modes.epoch, case.epoch.unwrap(), "{}", case.label());
        assert_eq!(modes.modes.len(), case.modes.len());
        let stride = case.sample_stride.unwrap();
        for expected in &case.modes {
            let mode = modes.mode(expected.ell, expected.emm).unwrap();
            assert_eq!(mode.len(), case.length.unwrap());
            let re: Vec<f64> = mode.iter().step_by(stride).map(|h| h.re).collect();
            let im: Vec<f64> = mode.iter().step_by(stride).map(|h| h.im).collect();
            // relative to the largest |h_lm| of the mode
            let scale = expected
                .real
                .iter()
                .zip(&expected.imag)
                .fold(0.0f64, |m, (r, i)| m.max(r.hypot(*i)));
            let diff = re
                .iter()
                .zip(&im)
                .zip(expected.real.iter().zip(&expected.imag))
                .fold(0.0f64, |m, ((a, b), (c, d))| m.max((a - c).hypot(b - d)));
            let err = if scale > 0.0 { diff / scale } else { diff };
            worst = worst.max(err);
            println!(
                "({}, {:>2}) max rel err {err:.2e} | {}",
                expected.ell,
                expected.emm,
                case.label()
            );
        }
    }
    println!(
        "worst {worst:.3e} (lalsimulation {})",
        reference.lalsimulation_version
    );
    assert!(worst < 1e-12, "worst relative error {worst:e}");
}

#[derive(Deserialize)]
struct DynamicsCase {
    q: f64,
    chi_a: [f64; 3],
    chi_b: [f64; 3],
    omega_ref: f64,
    init_quat: [f64; 4],
    init_orbphase: f64,
    t: Vec<f64>,
    q0: Vec<f64>,
    qx: Vec<f64>,
    qy: Vec<f64>,
    qz: Vec<f64>,
    orbphase: Vec<f64>,
    #[serde(rename = "chiAx")]
    chi_ax: Vec<f64>,
    #[serde(rename = "chiAy")]
    chi_ay: Vec<f64>,
    #[serde(rename = "chiAz")]
    chi_az: Vec<f64>,
    #[serde(rename = "chiBx")]
    chi_bx: Vec<f64>,
    #[serde(rename = "chiBy")]
    chi_by: Vec<f64>,
    #[serde(rename = "chiBz")]
    chi_bz: Vec<f64>,
}

#[test]
fn nrsur7dq4_dynamics_match_lalsimulation() {
    let Some(data) = load_data() else { return };
    let reference: ReferenceFile<DynamicsCase> = reference("nrsur7dq4_dynamics");
    let mut worst: f64 = 0.0;
    for case in &reference.cases {
        let dynamics = data
            .dynamics(
                case.q,
                case.chi_a,
                case.chi_b,
                case.omega_ref,
                case.init_quat,
                case.init_orbphase,
                false,
            )
            .unwrap();
        assert_eq!(dynamics.t, case.t);
        let pairs = [
            (&dynamics.quaternion[0], &case.q0),
            (&dynamics.quaternion[1], &case.qx),
            (&dynamics.quaternion[2], &case.qy),
            (&dynamics.quaternion[3], &case.qz),
            (&dynamics.orbital_phase, &case.orbphase),
            (&dynamics.chi_a[0], &case.chi_ax),
            (&dynamics.chi_a[1], &case.chi_ay),
            (&dynamics.chi_a[2], &case.chi_az),
            (&dynamics.chi_b[0], &case.chi_bx),
            (&dynamics.chi_b[1], &case.chi_by),
            (&dynamics.chi_b[2], &case.chi_bz),
        ];
        let mut errors = Vec::new();
        for (test, expected) in pairs {
            // absolute for these O(1) quantities
            let (err, identical) = compare(test, expected);
            let scale = expected.iter().fold(0.0f64, |m, v| m.max(v.abs()));
            errors.push((err * scale, identical));
        }
        let max = errors.iter().fold(0.0f64, |m, e| m.max(e.0));
        let identical = errors.iter().map(|e| e.1).sum::<f64>() / errors.len() as f64;
        println!(
            "q={} omega_ref={}: max abs err {max:.2e}, bit-identical {:.1}%",
            case.q,
            case.omega_ref,
            100.0 * identical
        );
        worst = worst.max(max);
    }
    assert!(worst < 1e-12, "worst error {worst:e}");
}

/// A reused workspace (with its cache of the inertial modes) gives the same
/// waveforms as fresh evaluations.
#[test]
fn nrsur7dq4_workspace_reuse_is_exact() {
    let Some(data) = load_data() else { return };
    let base = PrecessingBinary {
        mass_1: 40.0,
        mass_2: 20.0,
        spin_1: [0.4, -0.3, 0.2],
        spin_2: [0.1, 0.5, -0.3],
        luminosity_distance: 100.0,
        inclination: 0.7,
        phi_ref: 0.1,
        f_ref: 20.0,
    };
    let variations = [
        base,
        PrecessingBinary {
            inclination: 2.1,
            phi_ref: 4.0,
            luminosity_distance: 900.0,
            ..base
        },
        PrecessingBinary {
            mass_1: 20.0,
            mass_2: 40.0,
            ..base
        },
        PrecessingBinary {
            spin_1: [0.0, 0.0, 0.3],
            ..base
        },
        base,
    ];
    let mut workspace = NRSurWorkspace::new();
    let mut out = TimeSeriesPolarizations::default();
    let mut modes = TimeSeriesModes::default();
    for (k, binary) in variations.iter().enumerate() {
        let delta_t = if k % 2 == 0 {
            1.0 / 4096.0
        } else {
            1.0 / 2048.0
        };
        let model = NRSur7dq4::new(data, binary).unwrap();
        model
            .polarizations_with(&mut workspace, delta_t, 20.0, &mut out)
            .unwrap();
        let fresh = model.polarizations(delta_t, 20.0).unwrap();
        assert_eq!(out.epoch, fresh.epoch);
        assert_eq!(out.hplus, fresh.hplus, "case {k}");
        assert_eq!(out.hcross, fresh.hcross, "case {k}");
        model
            .modes_with(&mut workspace, delta_t, 20.0, &mut modes)
            .unwrap();
        let fresh = model.modes(delta_t, 20.0).unwrap();
        for (a, b) in modes.modes.iter().zip(&fresh.modes) {
            assert_eq!((a.ell, a.m), (b.ell, b.m));
            assert_eq!(a.data, b.data, "case {k} mode ({}, {})", a.ell, a.m);
        }
        // an error does not leave a stale cache behind
        assert!(model
            .polarizations_with(&mut workspace, delta_t, 5.0, &mut out)
            .is_err());
    }
}

//! The grid evaluation paths agree with the per-frequency evaluation.
use num_complex::Complex;
use waveformrs::{
    grid::FrequencyGrid,
    imrphenomd::IMRPhenomD,
    imrphenomx::{
        IMRPhenomXAS, IMRPhenomXHM, IMRPhenomXP, IMRPhenomXPHM, Polarizations, PrecessingBinary,
        PrecessionVersion,
    },
    taylor::TaylorF2,
    waveform::Waveform,
};

const PARAMETERS: [(f64, f64, f64, f64); 4] = [
    (60.0, 0.5, 0.3, -0.2),
    (12.0, 0.25, 0.8, 0.1),
    (150.0, 0.9, -0.7, 0.5),
    (3.0, 0.7, 0.05, 0.02),
];

fn max_relative_difference(a: &[Complex<f64>], b: &[Complex<f64>]) -> f64 {
    let scale = b.iter().map(|x| x.norm()).fold(0.0, f64::max);
    a.iter()
        .zip(b)
        .map(|(x, y)| (x - y).norm())
        .fold(0.0, f64::max)
        / scale
}

fn grids() -> Vec<FrequencyGrid> {
    vec![
        FrequencyGrid::uniform(20.0, 1024.0, 0.125),
        FrequencyGrid::new(
            &(0..500)
                .map(|i| 10.0 * 1.01_f64.powi(i))
                .collect::<Vec<_>>(),
        ),
    ]
}

fn check<W: Waveform>(name: &str, model: W, phi_c: f64, tolerance: f64) {
    for grid in grids() {
        let (plus, cross) = model.waveform_modes(grid.frequencies(), phi_c, 0.7);
        let zero = Complex::new(0.0, 0.0);
        let (mut grid_plus, mut grid_cross) = (vec![zero; grid.len()], vec![zero; grid.len()]);
        model.waveform_modes_on_grid(&grid, phi_c, 0.7, &mut grid_plus, &mut grid_cross);
        let error = max_relative_difference(&grid_plus, &plus)
            .max(max_relative_difference(&grid_cross, &cross));
        assert!(error < tolerance, "{name}: {error:e}");
    }
}

#[test]
fn aligned_spin_models_on_grids() {
    for (m, q, c1, c2) in PARAMETERS {
        check(
            "TaylorF2",
            TaylorF2::new(m, q, c1, c2, 400.0).with_reference_frequency(20.0),
            0.3,
            1e-10,
        );
        check(
            "IMRPhenomD",
            IMRPhenomD::new(m, q, c1, c2, 400.0).with_reference_frequency(20.0),
            0.3,
            1e-10,
        );
        check(
            "IMRPhenomXAS",
            IMRPhenomXAS::new(m, q, c1, c2, 400.0).with_reference_frequency(20.0),
            0.3,
            1e-10,
        );
        // the default (reference at the first frequency) path
        check(
            "IMRPhenomXAS, f_ref = f_min",
            IMRPhenomXAS::new(m, q, c1, c2, 400.0),
            0.3,
            1e-10,
        );
    }
}

#[test]
fn imrphenomxhm_on_grids() {
    for (m, q, c1, c2) in PARAMETERS {
        let model = IMRPhenomXHM::new(m, q, c1, c2, 400.0).with_reference_frequency(20.0);
        for grid in grids() {
            let (plus, cross) = model.polarizations(grid.frequencies(), 0.3, 0.7);
            let zero = Complex::new(0.0, 0.0);
            let (mut gp, mut gc) = (vec![zero; grid.len()], vec![zero; grid.len()]);
            model.polarizations_on_grid(&grid, 0.3, 0.7, &mut gp, &mut gc);
            let error =
                max_relative_difference(&gp, &plus).max(max_relative_difference(&gc, &cross));
            assert!(error < 1e-10, "IMRPhenomXHM M = {m}: {error:e}");
        }
    }
}

fn precessing_binaries() -> [PrecessingBinary; 2] {
    [
        // the binary of examples/benchmark.rs
        PrecessingBinary {
            mass_1: 40.0,
            mass_2: 20.0,
            spin_1: [0.4, -0.3, 0.2],
            spin_2: [0.1, 0.5, -0.3],
            luminosity_distance: 400.0,
            inclination: 0.4,
            phi_ref: 0.1,
            f_ref: 20.0,
        },
        // reference at the first frequency of the grid
        PrecessingBinary {
            mass_1: 29.0,
            mass_2: 36.0,
            spin_1: [-0.2, 0.3, -0.1],
            spin_2: [0.4, 0.2, 0.3],
            luminosity_distance: 500.0,
            inclination: 2.1,
            phi_ref: 1.2,
            f_ref: 0.0,
        },
    ]
}

/// A binary whose SpinTaylor integration does not reach the first
/// frequency, so that version 320 falls back to the MSA angles (for
/// IMRPhenomXPHM only with the (2,2) mode alone), with grids starting at
/// its reference frequency.
fn heavy_binary() -> (PrecessingBinary, Vec<FrequencyGrid>) {
    let binary = PrecessingBinary {
        mass_1: 180.0,
        mass_2: 65.0,
        spin_1: [0.2, 0.1, -0.1],
        spin_2: [-0.1, -0.2, 0.1],
        luminosity_distance: 500.0,
        inclination: 0.4,
        phi_ref: 0.3,
        f_ref: 15.5,
    };
    let grids = vec![
        FrequencyGrid::uniform(15.5, 200.0, 0.25),
        FrequencyGrid::new(
            &(0..300)
                .map(|i| 15.5 * 1.01_f64.powi(i))
                .collect::<Vec<_>>(),
        ),
    ];
    (binary, grids)
}

const VERSIONS: [PrecessionVersion; 3] = [
    PrecessionVersion::Nnlo,
    PrecessionVersion::Msa,
    PrecessionVersion::SpinTaylor,
];

/// The grid path `on_grid` agrees with the per-frequency polarisations
/// `expected` and vanishes where they do.
fn check_precessing(
    name: &str,
    grid: &FrequencyGrid,
    expected: Polarizations,
    on_grid: impl Fn(&FrequencyGrid, &mut [Complex<f64>], &mut [Complex<f64>]),
) {
    let (plus, cross) = expected;
    let zero = Complex::new(0.0, 0.0);
    let (mut gp, mut gc) = (vec![zero; grid.len()], vec![zero; grid.len()]);
    on_grid(grid, &mut gp, &mut gc);
    let error = max_relative_difference(&gp, &plus).max(max_relative_difference(&gc, &cross));
    assert!(error < 1e-10, "{name} on {} points: {error:e}", grid.len());
    for (index, h) in plus.iter().enumerate() {
        if h.norm() == 0.0 {
            assert_eq!(gp[index], zero, "{name}");
        }
    }
}

#[test]
fn imrphenomxp_on_grids() {
    for binary in precessing_binaries() {
        for version in VERSIONS {
            let model = IMRPhenomXP::with_version(&binary, version).unwrap();
            for grid in grids() {
                check_precessing(
                    &format!("IMRPhenomXP {version:?}"),
                    &grid,
                    model.polarizations(grid.frequencies()).unwrap(),
                    |grid, hp, hc| model.polarizations_on_grid(grid, hp, hc).unwrap(),
                );
            }
        }
    }
}

#[test]
fn imrphenomxphm_on_grids() {
    for binary in precessing_binaries() {
        for version in VERSIONS {
            let model = IMRPhenomXPHM::with_version(&binary, version).unwrap();
            let restricted = model
                .clone()
                .with_mode_array(&[(2, 1), (3, 3), (4, 4)])
                .unwrap();
            for grid in grids() {
                for model in [&model, &restricted] {
                    check_precessing(
                        &format!("IMRPhenomXPHM {version:?} {:?}", model.modes()),
                        &grid,
                        model.polarizations(grid.frequencies()).unwrap(),
                        |grid, hp, hc| model.polarizations_on_grid(grid, hp, hc).unwrap(),
                    );
                }
            }
        }
    }
}

/// The SpinTaylor angles fall back to the MSA angles on the grid path as on
/// the per-frequency path.
#[test]
fn spintaylor_fallback_on_grids() {
    let (binary, grids) = heavy_binary();
    let xp = IMRPhenomXP::new(&binary).unwrap();
    let xp_msa = IMRPhenomXP::with_version(&binary, PrecessionVersion::Msa).unwrap();
    let xphm = IMRPhenomXPHM::new(&binary).unwrap();
    let xphm_22 = xphm.clone().with_mode_array(&[(2, 2)]).unwrap();
    let xphm_22_msa = IMRPhenomXPHM::with_version(&binary, PrecessionVersion::Msa)
        .unwrap()
        .with_mode_array(&[(2, 2)])
        .unwrap();
    let zero = Complex::new(0.0, 0.0);
    for grid in &grids {
        let n = grid.len();
        check_precessing(
            "IMRPhenomXP 320 -> 223",
            grid,
            xp.polarizations(grid.frequencies()).unwrap(),
            |grid, hp, hc| xp.polarizations_on_grid(grid, hp, hc).unwrap(),
        );
        let (mut hp, mut hc) = (vec![zero; n], vec![zero; n]);
        let (mut hp_msa, mut hc_msa) = (vec![zero; n], vec![zero; n]);
        xp.polarizations_on_grid(grid, &mut hp, &mut hc).unwrap();
        xp_msa
            .polarizations_on_grid(grid, &mut hp_msa, &mut hc_msa)
            .unwrap();
        assert_eq!((&hp, &hc), (&hp_msa, &hc_msa));

        for model in [&xphm, &xphm_22] {
            check_precessing(
                &format!("IMRPhenomXPHM 320 {:?}", model.modes()),
                grid,
                model.polarizations(grid.frequencies()).unwrap(),
                |grid, hp, hc| model.polarizations_on_grid(grid, hp, hc).unwrap(),
            );
        }
        xphm_22
            .polarizations_on_grid(grid, &mut hp, &mut hc)
            .unwrap();
        xphm_22_msa
            .polarizations_on_grid(grid, &mut hp_msa, &mut hc_msa)
            .unwrap();
        assert_eq!((&hp, &hc), (&hp_msa, &hc_msa));
    }
}

#[test]
fn precessing_grid_paths_reject_invalid_grids() {
    let binary = precessing_binaries()[0];
    let xp = IMRPhenomXP::new(&binary).unwrap();
    let xphm = IMRPhenomXPHM::new(&binary).unwrap();
    let zero = Complex::new(0.0, 0.0);
    for frequencies in [vec![30.0, 20.0, 40.0], vec![-1.0, 20.0]] {
        let grid = FrequencyGrid::new(&frequencies);
        let (mut hp, mut hc) = (vec![zero; grid.len()], vec![zero; grid.len()]);
        assert_eq!(
            xp.polarizations_on_grid(&grid, &mut hp, &mut hc),
            Err(xp.polarizations(&frequencies).unwrap_err())
        );
        assert_eq!(
            xphm.polarizations_on_grid(&grid, &mut hp, &mut hc),
            Err(xphm.polarizations(&frequencies).unwrap_err())
        );
    }
    let empty = FrequencyGrid::new(&[]);
    assert_eq!(xp.polarizations_on_grid(&empty, &mut [], &mut []), Ok(()));
    assert_eq!(xphm.polarizations_on_grid(&empty, &mut [], &mut []), Ok(()));
}

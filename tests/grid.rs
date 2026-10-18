//! The grid evaluation paths agree with the per-frequency evaluation.
use num_complex::Complex;
use waveformrs::{
    grid::FrequencyGrid,
    imrphenomd::IMRPhenomD,
    imrphenomx::{IMRPhenomXAS, IMRPhenomXHM},
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
            0.6,
            1e-10,
        );
        check(
            "IMRPhenomXAS",
            IMRPhenomXAS::new(m, q, c1, c2, 400.0).with_reference_frequency(20.0),
            0.6,
            1e-10,
        );
        // the default (reference at the first frequency) path
        check(
            "IMRPhenomXAS, f_ref = f_min",
            IMRPhenomXAS::new(m, q, c1, c2, 400.0),
            0.6,
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

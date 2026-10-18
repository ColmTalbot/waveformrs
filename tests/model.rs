//! The common model interface (`waveformrs::model`): the trait methods
//! reproduce LALSimulation's `SimInspiralChooseFDWaveformSequence` (with the
//! reference data of `tests/lalsimulation.rs`) and agree with the
//! model-specific methods.
mod common;

use num_complex::Complex;
use waveformrs::{
    choose_fd_model, choose_td_model,
    grid::FrequencyGrid,
    imrphenomd::IMRPhenomD,
    imrphenomx::{
        IMRPhenomXAS, IMRPhenomXError, IMRPhenomXHM, IMRPhenomXP, IMRPhenomXPHM, PrecessingBinary,
        PrecessionVersion,
    },
    taylor::TaylorF2,
    waveform::Waveform,
    Approximant, Binary, Error, FrequencyDomainModel, WaveformOptions,
};

fn binary(case: &common::Case) -> Binary {
    Binary {
        mass_1: case.mass_1,
        mass_2: case.mass_2,
        spin_1: case.spin_1,
        spin_2: case.spin_2,
        luminosity_distance: case.luminosity_distance,
        inclination: case.inclination,
        phi_ref: case.phi_ref,
        f_ref: case.f_ref,
    }
}

fn options(case: &common::Case) -> WaveformOptions {
    let version = match case.lal_flags.get("PhenomXPrecVersion").copied() {
        Some(102.0) => PrecessionVersion::Nnlo,
        Some(223.0) => PrecessionVersion::Msa,
        _ => PrecessionVersion::SpinTaylor,
    };
    let mut options = WaveformOptions::default().with_precession_version(version);
    if !case.mode_array.is_empty() {
        options = options.with_mode_array(&case.mode_array);
    }
    options
}

/// `choose_fd_model(...).polarizations` against
/// `SimInspiralChooseFDWaveformSequence` for every suite on non-uniform
/// grids, at the tolerances of `tests/lalsimulation.rs`.
#[test]
fn choose_fd_model_matches_lalsimulation() {
    let suites = [
        ("taylorf2", Approximant::TaylorF2, 1e-10, 1e-14),
        ("imrphenomd", Approximant::IMRPhenomD, 1e-10, 1e-14),
        ("imrphenomxas", Approximant::IMRPhenomXAS, 1e-9, 1e-12),
        ("imrphenomxhm", Approximant::IMRPhenomXHM, 1e-9, 1e-12),
        ("imrphenomxp_nnlo", Approximant::IMRPhenomXP, 1e-9, 1e-12),
        (
            "imrphenomxp_spintaylor",
            Approximant::IMRPhenomXP,
            1e-9,
            1e-12,
        ),
        ("imrphenomxp_msa", Approximant::IMRPhenomXP, 1e-9, 1e-12),
        (
            "imrphenomxphm_nnlo",
            Approximant::IMRPhenomXPHM,
            1e-9,
            1e-12,
        ),
        (
            "imrphenomxphm_spintaylor",
            Approximant::IMRPhenomXPHM,
            1e-9,
            1e-12,
        ),
        ("imrphenomxphm_msa", Approximant::IMRPhenomXPHM, 1e-9, 1e-12),
        (
            "imrphenomxphm_modes",
            Approximant::IMRPhenomXPHM,
            1e-9,
            1e-12,
        ),
    ];
    for (suite, approximant, rel, mismatch) in suites {
        let reference = common::load(suite);
        assert!(reference
            .cases
            .iter()
            .all(|case| case.approximant == approximant.to_string()));
        common::check_all(suite, rel, mismatch, |case| {
            if case.frequencies[0] == 0.0 {
                // uniform-grid cases (SimInspiralChooseFDWaveform) are
                // covered by tests/lalsimulation.rs
                return (case.hplus(), case.hcross());
            }
            choose_fd_model(approximant, &binary(case), &options(case))
                .unwrap()
                .polarizations(&case.frequencies)
                .unwrap()
        });
    }
}

fn test_binary() -> Binary {
    Binary {
        mass_1: 36.0,
        mass_2: 29.0,
        spin_1: [0.0, 0.0, 0.3],
        spin_2: [0.0, 0.0, -0.4],
        luminosity_distance: 500.0,
        inclination: 0.7,
        phi_ref: 1.2,
        f_ref: 25.0,
    }
}

fn precessing_binary() -> Binary {
    Binary {
        spin_1: [0.4, 0.2, 0.3],
        spin_2: [-0.2, 0.3, -0.1],
        ..test_binary()
    }
}

fn frequencies() -> Vec<f64> {
    (0..300).map(|i| 20.0 * 1.01f64.powi(i)).collect()
}

/// The trait methods are the model-specific methods with the parameters of
/// the binary, bit for bit.
#[test]
fn trait_methods_are_the_model_methods() {
    let b = test_binary();
    let f = frequencies();
    let (m, q) = (b.total_mass(), b.mass_2 / b.mass_1);
    let (c1, c2) = (b.spin_1[2], b.spin_2[2]);
    let d = b.luminosity_distance;

    let tf2 = TaylorF2::new(m, q, c1, c2, d).with_reference_frequency(b.f_ref);
    assert_eq!(
        FrequencyDomainModel::polarizations(&TaylorF2::from_binary(&b).unwrap(), &f).unwrap(),
        tf2.waveform_modes(&f, b.phi_ref, b.inclination)
    );
    let imrd = IMRPhenomD::new(m, q, c1, c2, d).with_reference_frequency(b.f_ref);
    assert_eq!(
        FrequencyDomainModel::polarizations(&IMRPhenomD::from_binary(&b).unwrap(), &f).unwrap(),
        imrd.waveform_modes(&f, b.phi_ref, b.inclination)
    );
    let xas = IMRPhenomXAS::from_components(b.mass_1, b.mass_2, c1, c2, d)
        .unwrap()
        .with_reference_frequency(b.f_ref);
    assert_eq!(
        FrequencyDomainModel::polarizations(&IMRPhenomXAS::from_binary(&b).unwrap(), &f).unwrap(),
        xas.waveform_modes(&f, b.phi_ref, b.inclination)
    );
    let xhm = IMRPhenomXHM::from_components(b.mass_1, b.mass_2, c1, c2, d)
        .unwrap()
        .with_reference_frequency(b.f_ref);
    assert_eq!(
        FrequencyDomainModel::polarizations(&IMRPhenomXHM::from_binary(&b).unwrap(), &f).unwrap(),
        xhm.polarizations(&f, b.phi_ref, b.inclination)
    );

    let p = precessing_binary();
    assert_eq!(
        FrequencyDomainModel::polarizations(&IMRPhenomXP::from_binary(&p).unwrap(), &f).unwrap(),
        IMRPhenomXP::new(&p).unwrap().polarizations(&f).unwrap()
    );
    assert_eq!(
        FrequencyDomainModel::polarizations(&IMRPhenomXPHM::from_binary(&p).unwrap(), &f).unwrap(),
        IMRPhenomXPHM::new(&p).unwrap().polarizations(&f).unwrap()
    );
}

fn every_model(b: &Binary, p: &Binary) -> Vec<Box<dyn FrequencyDomainModel>> {
    Approximant::ALL
        .into_iter()
        .filter(Approximant::is_frequency_domain)
        .map(|approximant| {
            let binary = if approximant.is_precessing() { p } else { b };
            choose_fd_model(approximant, binary, &WaveformOptions::default()).unwrap()
        })
        .collect()
}

fn max_relative_difference(a: &[Complex<f64>], b: &[Complex<f64>]) -> f64 {
    let scale = b.iter().map(|x| x.norm()).fold(0.0, f64::max);
    a.iter()
        .zip(b)
        .map(|(x, y)| (x - y).norm())
        .fold(0.0, f64::max)
        / scale
}

/// `polarizations_into` is `polarizations` bit for bit, and the grid path
/// agrees to rounding, for every model (with `f_ref` set and with `f_ref =
/// 0`).
#[test]
fn into_and_grid_paths_agree() {
    let grid = FrequencyGrid::uniform(20.0, 1024.0, 0.25);
    for f_ref in [25.0, 0.0] {
        let b = Binary {
            f_ref,
            ..test_binary()
        };
        let p = Binary {
            f_ref,
            ..precessing_binary()
        };
        for model in every_model(&b, &p) {
            let f = grid.frequencies();
            let (hp, hc) = model.polarizations(f).unwrap();
            let zero = Complex::new(0.0, 0.0);
            let (mut hp_into, mut hc_into) = (vec![zero; f.len()], vec![zero; f.len()]);
            model
                .polarizations_into(f, &mut hp_into, &mut hc_into)
                .unwrap();
            assert_eq!((&hp_into, &hc_into), (&hp, &hc), "{}", model.approximant());
            let (mut hp_grid, mut hc_grid) = (vec![zero; f.len()], vec![zero; f.len()]);
            model
                .polarizations_on_grid(&grid, &mut hp_grid, &mut hc_grid)
                .unwrap();
            let error =
                max_relative_difference(&hp_grid, &hp).max(max_relative_difference(&hc_grid, &hc));
            assert!(
                error < 1e-10,
                "{} f_ref = {f_ref}: {error:e}",
                model.approximant()
            );
        }
    }
}

/// `f_ref = 0` references the phase to the first frequency of the grid
/// (to coalescence for TaylorF2), as in LALSimulation.
#[test]
fn zero_reference_frequency() {
    let b = Binary {
        f_ref: 0.0,
        ..test_binary()
    };
    let f = frequencies();
    let explicit = Binary { f_ref: f[0], ..b };
    for approximant in [
        Approximant::IMRPhenomD,
        Approximant::IMRPhenomXAS,
        Approximant::IMRPhenomXHM,
    ] {
        let options = WaveformOptions::default();
        assert_eq!(
            choose_fd_model(approximant, &b, &options)
                .unwrap()
                .polarizations(&f)
                .unwrap(),
            choose_fd_model(approximant, &explicit, &options)
                .unwrap()
                .polarizations(&f)
                .unwrap(),
            "{approximant}"
        );
    }
    let (m, q) = (b.total_mass(), b.mass_2 / b.mass_1);
    let coalescence = TaylorF2::new(m, q, 0.3, -0.4, 500.0).waveform_modes(&f, 1.2, 0.7);
    assert_eq!(
        choose_fd_model(Approximant::TaylorF2, &b, &WaveformOptions::default())
            .unwrap()
            .polarizations(&f)
            .unwrap(),
        coalescence
    );
}

/// The labels of the bodies do not matter (the models order them as LAL).
#[test]
fn body_labels_do_not_matter() {
    let b = test_binary();
    let swapped = Binary {
        mass_1: b.mass_2,
        mass_2: b.mass_1,
        spin_1: b.spin_2,
        spin_2: b.spin_1,
        ..b
    };
    let f = frequencies();
    for approximant in [Approximant::TaylorF2, Approximant::IMRPhenomD] {
        let options = WaveformOptions::default();
        let (hp, hc) = choose_fd_model(approximant, &b, &options)
            .unwrap()
            .polarizations(&f)
            .unwrap();
        let (hp_s, hc_s) = choose_fd_model(approximant, &swapped, &options)
            .unwrap()
            .polarizations(&f)
            .unwrap();
        assert!(max_relative_difference(&hp_s, &hp) < 1e-12, "{approximant}");
        assert!(max_relative_difference(&hc_s, &hc) < 1e-12, "{approximant}");
    }
}

#[test]
fn aligned_spin_models_reject_in_plane_spins() {
    let p = precessing_binary();
    for approximant in [
        Approximant::TaylorF2,
        Approximant::IMRPhenomD,
        Approximant::IMRPhenomXAS,
        Approximant::IMRPhenomXHM,
    ] {
        let err = choose_fd_model(approximant, &p, &WaveformOptions::default()).unwrap_err();
        assert!(
            matches!(err, Error::InPlaneSpins(a) if a == approximant),
            "{err}"
        );
        assert!(err.to_string().contains("non-zero transverse spins"));
    }
    for approximant in [Approximant::IMRPhenomXP, Approximant::IMRPhenomXPHM] {
        assert!(choose_fd_model(approximant, &p, &WaveformOptions::default()).is_ok());
    }
}

#[test]
fn invalid_parameters_are_errors() {
    let b = test_binary();
    let negative_mass = Binary { mass_2: -1.0, ..b };
    for approximant in [Approximant::TaylorF2, Approximant::IMRPhenomD] {
        assert!(matches!(
            choose_fd_model(approximant, &negative_mass, &WaveformOptions::default()),
            Err(Error::InvalidParameter { .. })
        ));
    }
    assert!(matches!(
        IMRPhenomXAS::from_binary(&negative_mass),
        Err(Error::IMRPhenomX(IMRPhenomXError::NonPositiveMass))
    ));
    let negative_f_ref = Binary { f_ref: -1.0, ..b };
    for approximant in Approximant::ALL
        .into_iter()
        .filter(Approximant::is_frequency_domain)
    {
        assert!(
            choose_fd_model(approximant, &negative_f_ref, &WaveformOptions::default()).is_err(),
            "{approximant}"
        );
    }
    let spin = Binary {
        spin_1: [0.0, 0.0, 1.2],
        ..b
    };
    assert!(IMRPhenomD::from_binary(&spin).is_err());
    assert!(IMRPhenomXHM::from_binary(&spin).is_err());
}

#[test]
fn options_are_checked() {
    let b = test_binary();
    let modes = WaveformOptions::default().with_mode_array(&[(2, 2), (2, -2)]);
    for approximant in [
        Approximant::TaylorF2,
        Approximant::IMRPhenomD,
        Approximant::IMRPhenomXAS,
        Approximant::IMRPhenomXP,
    ] {
        assert!(matches!(
            choose_fd_model(approximant, &b, &modes),
            Err(Error::Unsupported { .. })
        ));
    }
    let bad = WaveformOptions::default().with_mode_array(&[(5, 5)]);
    assert!(matches!(
        choose_fd_model(Approximant::IMRPhenomXHM, &b, &bad),
        Err(Error::IMRPhenomX(IMRPhenomXError::InvalidModeArray))
    ));
    assert!(matches!(
        choose_fd_model(Approximant::IMRPhenomXPHM, &precessing_binary(), &bad),
        Err(Error::IMRPhenomX(IMRPhenomXError::InvalidModeArray))
    ));
    assert!(matches!(
        choose_fd_model(Approximant::NRSur7dq4, &b, &WaveformOptions::default()),
        Err(Error::Unsupported { .. })
    ));
    assert!(matches!(
        choose_td_model(Approximant::IMRPhenomXAS, &b, &WaveformOptions::default()),
        Err(Error::Unsupported { .. })
    ));
}

/// A mode array of IMRPhenomXHM (with signs) and of IMRPhenomXPHM (either
/// sign selects both) gives the model-specific result.
#[test]
fn mode_arrays_select_modes() {
    let f = frequencies();
    let b = test_binary();
    let options = WaveformOptions::default().with_mode_array(&[(2, -2), (3, 3)]);
    let xhm = IMRPhenomXHM::from_components(b.mass_1, b.mass_2, 0.3, -0.4, 500.0)
        .unwrap()
        .with_reference_frequency(b.f_ref)
        .with_mode_array(&[(2, -2), (3, 3)]);
    assert_eq!(
        choose_fd_model(Approximant::IMRPhenomXHM, &b, &options)
            .unwrap()
            .polarizations(&f)
            .unwrap(),
        xhm.polarizations(&f, b.phi_ref, b.inclination)
    );
    let p = precessing_binary();
    let xphm = IMRPhenomXPHM::new(&p)
        .unwrap()
        .with_mode_array(&[(2, 2), (3, 3)])
        .unwrap();
    assert_eq!(
        choose_fd_model(Approximant::IMRPhenomXPHM, &p, &options)
            .unwrap()
            .polarizations(&f)
            .unwrap(),
        xphm.polarizations(&f).unwrap()
    );
}

#[test]
fn approximant_names() {
    for approximant in Approximant::ALL {
        let name = approximant.to_string();
        assert_eq!(name.parse::<Approximant>().unwrap(), approximant);
    }
    let err = "IMRPhenomXPHMM".parse::<Approximant>().unwrap_err();
    assert!(matches!(&err, Error::UnknownApproximant(name) if name == "IMRPhenomXPHMM"));
    assert!(err.to_string().contains("IMRPhenomXPHM"));
    assert!(Approximant::IMRPhenomXPHM.is_precessing());
    assert!(!Approximant::IMRPhenomXHM.is_precessing());
    assert!(!Approximant::NRSur7dq4.is_frequency_domain());
}

/// `PrecessingBinary` is `Binary`.
#[test]
fn precessing_binary_is_binary() {
    let p: PrecessingBinary = precessing_binary();
    let b: Binary = p;
    assert_eq!(b, p);
}

#[test]
fn errors_wrap_the_model_errors() {
    let err: Error = IMRPhenomXError::NonPositiveMass.into();
    assert_eq!(
        err.to_string(),
        IMRPhenomXError::NonPositiveMass.to_string()
    );
    assert!(std::error::Error::source(&err).is_some());
    // a boxed model is built with choose_fd_model
    assert!(matches!(
        <Box<dyn FrequencyDomainModel>>::from_binary(&test_binary()),
        Err(Error::InvalidRequest(_))
    ));
}

/// NRSur7dq4 through the time-domain interface is the model-specific API.
#[test]
fn nrsur7dq4_through_the_time_domain_interface() {
    if common::nrsur7dq4_data().is_none() {
        return;
    }
    let b = Binary {
        mass_1: 40.0,
        mass_2: 20.0,
        spin_1: [0.4, -0.3, 0.2],
        spin_2: [0.1, 0.5, -0.3],
        luminosity_distance: 100.0,
        inclination: 0.7,
        phi_ref: 0.1,
        f_ref: 20.0,
    };
    let model = choose_td_model(Approximant::NRSur7dq4, &b, &WaveformOptions::default()).unwrap();
    assert_eq!(model.approximant(), Approximant::NRSur7dq4);
    let h = model.polarizations(1.0 / 2048.0, 20.0).unwrap();
    let direct =
        waveformrs::nrsur::NRSur7dq4::new(waveformrs::nrsur::NRSur7dq4Data::global().unwrap(), &b)
            .unwrap();
    let expected = direct.polarizations(1.0 / 2048.0, 20.0).unwrap();
    assert_eq!(h.hplus, expected.hplus);
    assert_eq!(h.hcross, expected.hcross);
    assert_eq!(h.epoch, expected.epoch);
    let modes = model.modes(1.0 / 2048.0, 20.0).unwrap();
    let expected = direct.modes(1.0 / 2048.0, 20.0).unwrap();
    assert_eq!(modes.modes.len(), expected.modes.len());
    for (a, e) in modes.modes.iter().zip(&expected.modes) {
        assert_eq!((a.ell, a.m), (e.ell, e.m));
        assert_eq!(a.data, e.data);
    }
}

//! The switch for the automatic changes of the Euler-angle prescription of
//! IMRPhenomXP and IMRPhenomXPHM (`with_version_fallback`): by default the
//! models change the prescription as LALSimulation does (checked against
//! LAL's reference data); without fallbacks they use the requested
//! prescription or return `PrecessionVersionUnavailable`.
mod common;

use num_complex::Complex;
use waveformrs::{
    choose_fd_model,
    grid::FrequencyGrid,
    imrphenomx::{IMRPhenomXError, IMRPhenomXP, IMRPhenomXPHM, PrecessionVersion},
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

fn generic() -> Binary {
    Binary {
        mass_1: 36.0,
        mass_2: 29.0,
        spin_1: [0.4, 0.2, 0.3],
        spin_2: [-0.2, 0.3, -0.1],
        luminosity_distance: 500.0,
        inclination: 0.7,
        phi_ref: 1.2,
        f_ref: 20.0,
    }
}

fn frequencies() -> Vec<f64> {
    (0..300).map(|i| 20.0 * 1.01f64.powi(i)).collect()
}

fn unavailable(requested: u32, fallback: u32) -> impl Fn(&IMRPhenomXError) -> bool {
    move |err| {
        matches!(err, IMRPhenomXError::PrecessionVersionUnavailable { requested: r, fallback: f, .. }
            if *r == requested && *f == fallback)
    }
}

/// Both models with and without fallbacks.
fn models(b: &Binary, version: PrecessionVersion, allow: bool) -> (IMRPhenomXP, IMRPhenomXPHM) {
    (
        IMRPhenomXP::with_version(b, version)
            .unwrap()
            .with_version_fallback(allow),
        IMRPhenomXPHM::with_version(b, version)
            .unwrap()
            .with_version_fallback(allow),
    )
}

/// Without fallbacks, a binary for which the requested prescription works
/// gives the same polarizations, grid polarizations, uniform polarizations
/// and modes.
#[test]
fn no_fallback_needed_gives_the_same_waveforms() {
    let f = frequencies();
    let grid = FrequencyGrid::new(&f);
    let zero = Complex::new(0.0, 0.0);
    for version in [
        PrecessionVersion::SpinTaylor,
        PrecessionVersion::Msa,
        PrecessionVersion::Nnlo,
    ] {
        for f_ref in [20.0, 0.0] {
            let b = Binary { f_ref, ..generic() };
            let (xp, xphm) = models(&b, version, true);
            let (xp_strict, xphm_strict) = models(&b, version, false);
            assert!(xp.version_fallback() && !xp_strict.version_fallback());
            assert!(xphm.version_fallback() && !xphm_strict.version_fallback());
            assert_eq!(
                xp_strict.polarizations(&f).unwrap(),
                xp.polarizations(&f).unwrap()
            );
            assert_eq!(
                xphm_strict.polarizations(&f).unwrap(),
                xphm.polarizations(&f).unwrap()
            );
            assert_eq!(
                xp_strict.polarizations_uniform(1.0, 20.0, 256.0).unwrap(),
                xp.polarizations_uniform(1.0, 20.0, 256.0).unwrap()
            );
            assert_eq!(
                xphm_strict.polarizations_uniform(1.0, 20.0, 256.0).unwrap(),
                xphm.polarizations_uniform(1.0, 20.0, 256.0).unwrap()
            );
            let (mut hp, mut hc) = (vec![zero; f.len()], vec![zero; f.len()]);
            let (mut hp_strict, mut hc_strict) = (vec![zero; f.len()], vec![zero; f.len()]);
            xphm.polarizations_on_grid(&grid, &mut hp, &mut hc).unwrap();
            xphm_strict
                .polarizations_on_grid(&grid, &mut hp_strict, &mut hc_strict)
                .unwrap();
            assert_eq!((hp_strict, hc_strict), (hp, hc));
            assert_eq!(
                FrequencyDomainModel::modes(&xphm_strict, &f).unwrap(),
                FrequencyDomainModel::modes(&xphm, &f).unwrap()
            );
            assert_eq!(
                FrequencyDomainModel::modes(&xp_strict, &f).unwrap(),
                FrequencyDomainModel::modes(&xp, &f).unwrap()
            );
        }
    }
}

/// SpinTaylor angles with negligible in-plane spins: LAL (and the default)
/// uses the NNLO angles; without fallbacks this is an error.
#[test]
fn negligible_in_plane_spins() {
    let f = frequencies();
    for in_plane in [0.0, 5e-8] {
        let b = Binary {
            spin_1: [in_plane, 0.0, 0.5],
            spin_2: [0.0, 0.0, -0.3],
            ..generic()
        };
        let (xp, xphm) = models(&b, PrecessionVersion::SpinTaylor, true);
        let (xp_nnlo, xphm_nnlo) = models(&b, PrecessionVersion::Nnlo, true);
        assert_eq!(
            xp.polarizations(&f).unwrap(),
            xp_nnlo.polarizations(&f).unwrap()
        );
        assert_eq!(
            xphm.polarizations(&f).unwrap(),
            xphm_nnlo.polarizations(&f).unwrap()
        );
        let (xp, xphm) = models(&b, PrecessionVersion::SpinTaylor, false);
        let is = unavailable(320, 102);
        assert!(is(&xp.polarizations(&f).unwrap_err()));
        assert!(is(&xp.polarizations_uniform(1.0, 20.0, 256.0).unwrap_err()));
        assert!(is(&xphm.polarizations(&f).unwrap_err()));
        assert!(is(&xphm
            .polarizations_uniform(1.0, 20.0, 256.0)
            .unwrap_err()));
        let zero = Complex::new(0.0, 0.0);
        let grid = FrequencyGrid::new(&f);
        let (mut hp, mut hc) = (vec![zero; f.len()], vec![zero; f.len()]);
        assert!(is(&xp
            .polarizations_on_grid(&grid, &mut hp, &mut hc)
            .unwrap_err()));
        assert!(is(&xphm
            .polarizations_on_grid(&grid, &mut hp, &mut hc)
            .unwrap_err()));
        for err in [
            FrequencyDomainModel::modes(&xp, &f).unwrap_err(),
            FrequencyDomainModel::modes(&xphm, &f).unwrap_err(),
        ] {
            assert!(matches!(&err, Error::IMRPhenomX(e) if is(e)), "{err}");
            assert!(err.to_string().contains("fallbacks are disabled"));
        }
    }
}

/// Heavy binaries whose SpinTaylor angles cannot be computed: LAL (and the
/// default) uses the MSA angles (reference data of
/// `imrphenomx*_spintaylor_fallback`); without fallbacks this is an error.
#[test]
fn spintaylor_angles_unavailable() {
    for (suite, approximant) in [
        ("imrphenomxp_spintaylor_fallback", Approximant::IMRPhenomXP),
        (
            "imrphenomxphm_spintaylor_fallback",
            Approximant::IMRPhenomXPHM,
        ),
    ] {
        // the default reproduces LAL's fallback
        common::check_all(suite, 1e-9, 1e-12, |case| {
            if case.frequencies[0] == 0.0 {
                return (case.hplus(), case.hcross());
            }
            let options = WaveformOptions::default().with_version_fallback(true);
            choose_fd_model(approximant, &binary(case), &options)
                .unwrap()
                .polarizations(&case.frequencies)
                .unwrap()
        });
        let reference = common::load(suite);
        for case in &reference.cases {
            let b = binary(case);
            let options = WaveformOptions::default().with_version_fallback(false);
            let model = choose_fd_model(approximant, &b, &options).unwrap();
            let is = unavailable(320, 223);
            let err = if case.frequencies[0] == 0.0 {
                let delta_f = case.frequencies[1];
                let result = match approximant {
                    Approximant::IMRPhenomXP => IMRPhenomXP::new(&b)
                        .unwrap()
                        .with_version_fallback(false)
                        .polarizations_uniform(delta_f, case.f_min, case.f_max),
                    _ => IMRPhenomXPHM::new(&b)
                        .unwrap()
                        .with_version_fallback(false)
                        .polarizations_uniform(delta_f, case.f_min, case.f_max),
                };
                Error::IMRPhenomX(result.unwrap_err())
            } else {
                let grid = FrequencyGrid::new(&case.frequencies);
                let zero = Complex::new(0.0, 0.0);
                let n = grid.len();
                let (mut hp, mut hc) = (vec![zero; n], vec![zero; n]);
                let grid_err = model
                    .polarizations_on_grid(&grid, &mut hp, &mut hc)
                    .unwrap_err();
                assert!(matches!(&grid_err, Error::IMRPhenomX(e) if is(e)));
                model.polarizations(&case.frequencies).unwrap_err()
            };
            assert!(
                matches!(&err, Error::IMRPhenomX(e) if is(e)),
                "{err} | {}",
                case.label()
            );
        }
    }
}

/// A binary whose MSA system cannot be initialised (in the reference data
/// of `imrphenomxp_msa`): LAL (and the default) uses the NNLO angles;
/// without fallbacks this is an error.
#[test]
fn msa_system_unavailable() {
    let reference = common::load("imrphenomxp_msa");
    let mut failing = 0;
    for case in &reference.cases {
        if case.frequencies[0] == 0.0 {
            continue;
        }
        let b = binary(case);
        let (xp, xphm) = models(&b, PrecessionVersion::Msa, false);
        match xp.polarizations(&case.frequencies) {
            Ok(h) => {
                let (default, _) = models(&b, PrecessionVersion::Msa, true);
                assert_eq!(h, default.polarizations(&case.frequencies).unwrap());
            }
            Err(err) => {
                assert!(unavailable(223, 102)(&err), "{err}");
                assert!(unavailable(223, 102)(
                    &xphm.polarizations(&case.frequencies).unwrap_err()
                ));
                // the default falls back to the NNLO angles, as LAL
                let (xp_default, _) = models(&b, PrecessionVersion::Msa, true);
                let (xp_nnlo, _) = models(&b, PrecessionVersion::Nnlo, true);
                assert_eq!(
                    xp_default.polarizations(&case.frequencies).unwrap(),
                    xp_nnlo.polarizations(&case.frequencies).unwrap()
                );
                failing += 1;
            }
        }
    }
    assert_eq!(failing, 1);
}

#[test]
fn options_default_to_lal_behaviour() {
    let options = WaveformOptions::default();
    assert!(options.version_fallback);
    assert!(WaveformOptions::lalsimulation_defaults().version_fallback);
    let model = IMRPhenomXPHM::new(&generic()).unwrap();
    assert!(model.version_fallback());
    // the flag survives a change of the mode array
    let strict = model
        .with_version_fallback(false)
        .with_mode_array(&[(2, 2), (3, 3)])
        .unwrap();
    assert!(!strict.version_fallback());
}

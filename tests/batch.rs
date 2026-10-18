//! Batched evaluation: one thread is exactly a serial loop, and the results
//! do not depend on the number of threads.
mod common;

use num_complex::Complex;
use waveformrs::{
    batch::Batch,
    choose_fd_model,
    grid::FrequencyGrid,
    imrphenomx::{IMRPhenomXAS, IMRPhenomXPHM, PrecessionVersion},
    model::TimeSeriesPolarizations,
    nrsur::{NRSur7dq4, NRSurWorkspace},
    Approximant, Binary, Error, FrequencyDomainModel, TimeDomainModel, WaveformOptions,
};

type C64 = Complex<f64>;

/// Deterministic pseudo-random binaries.
fn binaries(n: usize, precessing: bool) -> Vec<Binary> {
    let mut state = 0x2545_f491_4f6c_dd1du64;
    let mut uniform = move || {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        (state >> 11) as f64 / (1u64 << 53) as f64
    };
    (0..n)
        .map(|_| {
            let mass_1 = 20.0 + 40.0 * uniform();
            let mass_2 = mass_1 * (0.3 + 0.7 * uniform());
            let mut spin = || {
                let s = [
                    0.6 * uniform() - 0.3,
                    0.6 * uniform() - 0.3,
                    1.2 * uniform() - 0.6,
                ];
                if precessing {
                    s
                } else {
                    [0.0, 0.0, s[2]]
                }
            };
            Binary {
                mass_1,
                mass_2,
                spin_1: spin(),
                spin_2: spin(),
                luminosity_distance: 100.0 + 900.0 * uniform(),
                inclination: 3.0 * uniform(),
                phi_ref: 6.0 * uniform(),
                f_ref: 20.0,
            }
        })
        .collect()
}

fn buffers(n: usize, grid: &FrequencyGrid) -> (Vec<C64>, Vec<C64>) {
    let zero = C64::new(0.0, 0.0);
    (vec![zero; n * grid.len()], vec![zero; n * grid.len()])
}

/// One thread equals a serial loop over the models, and 2, 3 and 8 threads
/// equal one thread, bit for bit.
#[test]
fn grid_batches_do_not_depend_on_threads() {
    let grid = FrequencyGrid::uniform(20.0, 512.0, 0.5);
    let aligned = binaries(13, false);
    let precessing = binaries(7, true);
    let n_grid = grid.len();

    // serial loop
    let (mut hp_loop, mut hc_loop) = buffers(aligned.len(), &grid);
    for (i, b) in aligned.iter().enumerate() {
        IMRPhenomXAS::from_binary(b)
            .unwrap()
            .polarizations_on_grid(
                &grid,
                &mut hp_loop[i * n_grid..(i + 1) * n_grid],
                &mut hc_loop[i * n_grid..(i + 1) * n_grid],
            )
            .unwrap();
    }
    for threads in [1, 2, 3, 8] {
        let (mut hp, mut hc) = buffers(aligned.len(), &grid);
        let results = Batch::new(threads).polarizations_on_grid(
            &aligned,
            IMRPhenomXAS::from_binary,
            &grid,
            &mut hp,
            &mut hc,
        );
        assert!(results.iter().all(Result::is_ok));
        assert_eq!((&hp, &hc), (&hp_loop, &hc_loop), "{threads} threads");
    }

    let options = WaveformOptions::default();
    let build = |b: &Binary| choose_fd_model(Approximant::IMRPhenomXPHM, b, &options);
    let (mut hp_serial, mut hc_serial) = buffers(precessing.len(), &grid);
    Batch::serial().polarizations_on_grid(
        &precessing,
        build,
        &grid,
        &mut hp_serial,
        &mut hc_serial,
    );
    let models: Vec<IMRPhenomXPHM> = precessing
        .iter()
        .map(|b| IMRPhenomXPHM::new(b).unwrap())
        .collect();
    for threads in [2, 3, 8] {
        let (mut hp, mut hc) = buffers(precessing.len(), &grid);
        let results =
            Batch::new(threads).polarizations_on_grid(&precessing, build, &grid, &mut hp, &mut hc);
        assert!(results.iter().all(Result::is_ok));
        assert_eq!((&hp, &hc), (&hp_serial, &hc_serial), "{threads} threads");
        let (mut hp, mut hc) = buffers(precessing.len(), &grid);
        Batch::new(threads).models_on_grid(&models, &grid, &mut hp, &mut hc);
        assert_eq!((&hp, &hc), (&hp_serial, &hc_serial), "{threads} threads");
    }
}

#[test]
fn polarization_batches_do_not_depend_on_threads() {
    let frequencies: Vec<f64> = (0..200).map(|i| 20.0 * 1.012f64.powi(i)).collect();
    let precessing = binaries(9, true);
    let build = |b: &Binary| {
        let options = WaveformOptions::default().with_precession_version(PrecessionVersion::Msa);
        choose_fd_model(Approximant::IMRPhenomXP, b, &options)
    };
    let serial: Vec<_> = Batch::serial()
        .polarizations(&precessing, build, &frequencies)
        .into_iter()
        .map(Result::unwrap)
        .collect();
    for (b, h) in precessing.iter().zip(&serial) {
        assert_eq!(&build(b).unwrap().polarizations(&frequencies).unwrap(), h);
    }
    for threads in [2, 4, 16] {
        let parallel: Vec<_> = Batch::new(threads)
            .polarizations(&precessing, build, &frequencies)
            .into_iter()
            .map(Result::unwrap)
            .collect();
        assert_eq!(parallel, serial, "{threads} threads");
    }
}

/// A binary LAL rejects fails on its own; its row is zero.
#[test]
fn failures_are_per_binary() {
    let grid = FrequencyGrid::uniform(20.0, 256.0, 1.0);
    let mut aligned = binaries(6, false);
    aligned[3].mass_2 = -1.0;
    aligned[4].spin_1[0] = 0.1; // in-plane spin for an aligned-spin model
    for threads in [1, 3] {
        let (mut hp, mut hc) = buffers(aligned.len(), &grid);
        hp.fill(C64::new(1.0, 1.0));
        let results = Batch::new(threads).polarizations_on_grid(
            &aligned,
            IMRPhenomXAS::from_binary,
            &grid,
            &mut hp,
            &mut hc,
        );
        for (i, result) in results.iter().enumerate() {
            match i {
                3 => assert!(matches!(result, Err(Error::IMRPhenomX(_)))),
                4 => assert!(matches!(result, Err(Error::InPlaneSpins(_)))),
                _ => assert!(result.is_ok()),
            }
        }
        let n = grid.len();
        assert!(hp[3 * n..5 * n].iter().all(|h| h.norm() == 0.0));
        assert!(hp[2 * n..3 * n].iter().all(|h| h.norm() > 0.0));
    }
}

#[test]
fn thread_counts() {
    assert_eq!(Batch::serial().threads(), 1);
    assert_eq!(Batch::default().threads(), 1);
    assert_eq!(Batch::new(5).threads(), 5);
    assert!(Batch::new(0).threads() >= 1);
}

/// The time-domain batch of NRSur7dq4 (one workspace per thread) does not
/// depend on the number of threads and equals the direct evaluation.
#[test]
fn nrsur7dq4_batches_do_not_depend_on_threads() {
    let Some(data) = common::nrsur7dq4_data() else {
        return;
    };
    // from the start of the surrogate
    let mut precessing: Vec<Binary> = binaries(6, true)
        .into_iter()
        .map(|b| Binary { f_ref: 0.0, ..b })
        .collect();
    // the same intrinsic parameters twice in a row (the workspace cache)
    precessing[2] = Binary {
        inclination: 0.3,
        ..precessing[1]
    };
    let build = |b: &Binary| Ok(NRSur7dq4::new(data, b)?);
    let (delta_t, f_min) = (1.0 / 2048.0, 0.0);
    let mut serial = vec![TimeSeriesPolarizations::default(); precessing.len()];
    let results =
        Batch::serial().time_domain_polarizations(&precessing, build, delta_t, f_min, &mut serial);
    for result in &results {
        assert!(result.is_ok(), "{result:?}");
    }
    for (b, h) in precessing.iter().zip(&serial) {
        let direct = NRSur7dq4::new(data, b)
            .unwrap()
            .polarizations(delta_t, f_min)
            .unwrap();
        assert_eq!(
            (&h.hplus, &h.hcross, h.epoch),
            (&direct.hplus, &direct.hcross, direct.epoch)
        );
    }
    let models: Vec<NRSur7dq4> = precessing
        .iter()
        .map(|b| NRSur7dq4::new(data, b).unwrap())
        .collect();
    for threads in [2, 3, 8] {
        let mut out = vec![TimeSeriesPolarizations::default(); precessing.len()];
        Batch::new(threads).time_domain_polarizations(
            &precessing,
            |b: &Binary| NRSur7dq4::from_binary(b),
            delta_t,
            f_min,
            &mut out,
        );
        let mut from_models = vec![TimeSeriesPolarizations::default(); precessing.len()];
        Batch::new(threads).time_domain_models(&models, delta_t, f_min, &mut from_models);
        for ((a, b), s) in out.iter().zip(&from_models).zip(&serial) {
            assert_eq!(
                (&a.hplus, &a.hcross, a.epoch),
                (&s.hplus, &s.hcross, s.epoch)
            );
            assert_eq!(
                (&b.hplus, &b.hcross, b.epoch),
                (&s.hplus, &s.hcross, s.epoch)
            );
        }
    }
    // the workspace type of the trait
    let _: <NRSur7dq4 as TimeDomainModel>::Workspace = NRSurWorkspace::new();
}

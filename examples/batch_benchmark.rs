//! Throughput of the batched evaluation (`waveformrs::batch::Batch`) for 1,
//! 2, 4, 8 and all available threads.
//!
//! Each model is built from a binary and evaluated on its grid path inside
//! the batch, on the 8033 frequencies from 20 to 1024 Hz with
//! `delta_f = 0.125 Hz` (NRSur7dq4: time domain at 4096 Hz from the start
//! of the surrogate, one workspace per thread). The binaries are drawn as in
//! `examples/serial_benchmark.rs`: total mass log-uniform in [20, 150]
//! solar masses, mass ratio in [1/6, 1], isotropic spins with
//! `|chi| <= 0.9` (only the aligned components for the aligned-spin
//! models), isotropic inclination, distance 100-3000 Mpc, `f_ref = 20 Hz`.
//!
//! The numbers depend on the load of the machine; on a laptop the
//! efficiency cores and thermal limits make the scaling sub-linear.
//!
//!     cargo run --release --example batch_benchmark [-- scale]
//!
//! `scale` (default 1) multiplies the number of waveforms per model.
//! NRSur7dq4 needs its data file (`NRSUR7DQ4_DATA` or `LAL_DATA_PATH`).
use std::time::Instant;

use num_complex::Complex;
use waveformrs::{
    batch::Batch,
    grid::FrequencyGrid,
    imrphenomx::{IMRPhenomXAS, IMRPhenomXHM, IMRPhenomXP, IMRPhenomXPHM, PrecessionVersion},
    model::TimeSeriesPolarizations,
    nrsur::{NRSur7dq4, NRSur7dq4Data},
    Binary, FrequencyDomainModel, Result,
};

/// xorshift64* uniform deviates in [0, 1).
struct Rng(u64);

impl Rng {
    fn uniform(&mut self) -> f64 {
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        let x = self.0.wrapping_mul(0x2545_f491_4f6c_dd1d);
        (x >> 11) as f64 / (1u64 << 53) as f64
    }

    /// An isotropic vector of magnitude uniform in [0, max).
    fn spin(&mut self, max: f64) -> [f64; 3] {
        let magnitude = max * self.uniform();
        let cos_theta = 2.0 * self.uniform() - 1.0;
        let sin_theta = (1.0 - cos_theta * cos_theta).sqrt();
        let phi = 2.0 * std::f64::consts::PI * self.uniform();
        [
            magnitude * sin_theta * phi.cos(),
            magnitude * sin_theta * phi.sin(),
            magnitude * cos_theta,
        ]
    }
}

fn draws(n: usize, precessing: bool, seed: u64) -> Vec<Binary> {
    let mut rng = Rng(seed);
    (0..n)
        .map(|_| {
            let total_mass = 20.0 * (150.0f64 / 20.0).powf(rng.uniform());
            let q = 1.0 / 6.0 + (1.0 - 1.0 / 6.0) * rng.uniform();
            let mass_1 = total_mass / (1.0 + q);
            let (mut spin_1, mut spin_2) = (rng.spin(0.9), rng.spin(0.9));
            if !precessing {
                spin_1 = [0.0, 0.0, spin_1[2]];
                spin_2 = [0.0, 0.0, spin_2[2]];
            }
            Binary {
                mass_1,
                mass_2: total_mass - mass_1,
                spin_1,
                spin_2,
                luminosity_distance: 100.0 + 2900.0 * rng.uniform(),
                inclination: (2.0 * rng.uniform() - 1.0).acos(),
                phi_ref: 2.0 * std::f64::consts::PI * rng.uniform(),
                f_ref: 20.0,
            }
        })
        .collect()
}

/// Waveforms per second of `run` (which evaluates `n` waveforms and returns
/// the number of failures), the best of `repeats` runs.
fn throughput(n: usize, repeats: usize, mut run: impl FnMut() -> usize) -> (f64, usize) {
    let mut best = f64::INFINITY;
    let mut failures = 0;
    for _ in 0..repeats {
        let start = Instant::now();
        failures = run();
        best = best.min(start.elapsed().as_secs_f64());
    }
    (n as f64 / best, failures)
}

fn report(name: &str, n: usize, rates: &[(f64, usize)]) {
    let base = rates[0].0;
    let cells: Vec<String> = rates
        .iter()
        .map(|&(rate, _)| format!("{rate:>9.0} ({:>4.1}x)", rate / base))
        .collect();
    let failures: usize = rates.iter().map(|&(_, f)| f).max().unwrap_or(0);
    println!(
        "{name:<26} {n:>5} {}{}",
        cells.join(" "),
        if failures > 0 {
            format!("  [{failures} failed]")
        } else {
            String::new()
        }
    );
}

fn count_failures(results: &[Result<()>]) -> usize {
    results.iter().filter(|r| r.is_err()).count()
}

fn main() {
    let scale: usize = std::env::args()
        .nth(1)
        .map(|s| s.parse().expect("scale must be an integer"))
        .unwrap_or(1);
    let grid = FrequencyGrid::uniform(20.0, 1024.0, 0.125);
    assert_eq!(grid.len(), 8033);
    let all = Batch::new(0).threads();
    let mut threads = vec![1, 2, 4, 8, all];
    threads.sort_unstable();
    threads.dedup();

    println!(
        "waveforms per second (speed-up over 1 thread) on {} frequencies, {all} threads available",
        grid.len()
    );
    let header: Vec<String> = threads
        .iter()
        .map(|t| {
            format!(
                "{:>17}",
                format!("{t} thread{}", if *t == 1 { "" } else { "s" })
            )
        })
        .collect();
    println!("{:<26} {:>5} {}", "model", "n", header.join(""));

    fn fd<M: FrequencyDomainModel>(
        name: &str,
        n: usize,
        precessing: bool,
        grid: &FrequencyGrid,
        threads: &[usize],
        build: impl Fn(&Binary) -> Result<M> + Sync,
    ) {
        let binaries = draws(n, precessing, 0x5eed_0001);
        let zero = Complex::new(0.0, 0.0);
        let mut hp = vec![zero; n * grid.len()];
        let mut hc = vec![zero; n * grid.len()];
        // warm up
        let warm = n.min(8);
        Batch::serial().polarizations_on_grid(
            &binaries[..warm],
            &build,
            grid,
            &mut hp[..warm * grid.len()],
            &mut hc[..warm * grid.len()],
        );
        let rates: Vec<(f64, usize)> = threads
            .iter()
            .map(|&t| {
                throughput(n, 3, || {
                    count_failures(
                        &Batch::new(t)
                            .polarizations_on_grid(&binaries, &build, grid, &mut hp, &mut hc),
                    )
                })
            })
            .collect();
        report(name, n, &rates);
    }

    fd("IMRPhenomXAS", 2000 * scale, false, &grid, &threads, |b| {
        IMRPhenomXAS::from_binary(b)
    });
    fd("IMRPhenomXHM", 400 * scale, false, &grid, &threads, |b| {
        IMRPhenomXHM::from_binary(b)
    });
    fd(
        "IMRPhenomXP (SpinTaylor)",
        400 * scale,
        true,
        &grid,
        &threads,
        IMRPhenomXP::from_binary,
    );
    fd(
        "IMRPhenomXPHM (SpinTaylor)",
        200 * scale,
        true,
        &grid,
        &threads,
        IMRPhenomXPHM::from_binary,
    );
    fd(
        "IMRPhenomXPHM (MSA)",
        200 * scale,
        true,
        &grid,
        &threads,
        |b| Ok(IMRPhenomXPHM::with_version(b, PrecessionVersion::Msa)?),
    );

    match NRSur7dq4Data::global() {
        Ok(data) => {
            let n = 200 * scale;
            // reference at the start of the surrogate (f_ref = f_min = 0)
            let binaries: Vec<Binary> = draws(n, true, 0x5eed_0002)
                .into_iter()
                .map(|b| Binary { f_ref: 0.0, ..b })
                .collect();
            let mut out = vec![TimeSeriesPolarizations::default(); n];
            let build = |b: &Binary| Ok(NRSur7dq4::new(data, b)?);
            let rates: Vec<(f64, usize)> = threads
                .iter()
                .map(|&t| {
                    throughput(n, 3, || {
                        count_failures(&Batch::new(t).time_domain_polarizations(
                            &binaries,
                            build,
                            1.0 / 4096.0,
                            0.0,
                            &mut out,
                        ))
                    })
                })
                .collect();
            report("NRSur7dq4 (time domain)", n, &rates);
        }
        Err(err) => println!("NRSur7dq4 skipped: {err}"),
    }
}

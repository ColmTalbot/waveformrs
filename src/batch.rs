// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2024-2026 Colm Talbot. Parts of waveformrs are translated
// from LALSuite, GSL, glibc and fdlibm: see NOTICE.md for their authors.
//! Many waveforms at once, optionally on several threads.
//!
//! [`Batch`] evaluates a model for many binaries (or many prepared models)
//! on one [`FrequencyGrid`] (or at one sampling interval for the
//! time-domain models), writing into caller-provided buffers or returning a
//! `Vec` of results. With [`Batch::new`]`(threads)` the work is split into
//! chunks of consecutive binaries that `threads` scoped threads
//! (`std::thread::scope`, no extra dependency) take in turn; each thread has
//! its own buffers or workspace (one [`crate::nrsur::NRSurWorkspace`] per
//! thread for NRSur7dq4).
//!
//! The results do not depend on the number of threads: every waveform is
//! computed by the same code from its own inputs, and with one thread the
//! binaries are evaluated in order on the calling thread, exactly as a
//! serial loop.
//!
//! ```
//! use num_complex::Complex;
//! use waveformrs::{
//!     batch::Batch, grid::FrequencyGrid, imrphenomx::IMRPhenomXAS, Binary, FrequencyDomainModel,
//! };
//!
//! let grid = FrequencyGrid::uniform(20.0, 512.0, 0.25);
//! let binaries: Vec<Binary> = (0..8)
//!     .map(|i| Binary {
//!         mass_1: 30.0 + i as f64,
//!         mass_2: 25.0,
//!         spin_1: [0.0, 0.0, 0.2],
//!         spin_2: [0.0, 0.0, -0.1],
//!         luminosity_distance: 400.0,
//!         inclination: 0.4,
//!         phi_ref: 0.3,
//!         f_ref: 20.0,
//!     })
//!     .collect();
//! let zero = Complex::new(0.0, 0.0);
//! let mut hplus = vec![zero; binaries.len() * grid.len()];
//! let mut hcross = hplus.clone();
//! // waveform i in hplus[i * grid.len()..(i + 1) * grid.len()]
//! let results = Batch::new(4).polarizations_on_grid(
//!     &binaries,
//!     IMRPhenomXAS::from_binary,
//!     &grid,
//!     &mut hplus,
//!     &mut hcross,
//! );
//! assert!(results.iter().all(Result::is_ok));
//! ```

use std::sync::Mutex;

use num_complex::Complex;

use crate::error::{Error, Result};
use crate::grid::FrequencyGrid;
use crate::imrphenomx::Polarizations;
use crate::model::{Binary, FrequencyDomainModel, TimeDomainModel, TimeSeriesPolarizations};

/// Settings of a batch evaluation: the number of threads.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Batch {
    threads: usize,
}

impl Default for Batch {
    /// One thread.
    fn default() -> Self {
        Self::serial()
    }
}

impl Batch {
    /// Evaluate on `threads` threads; `0` selects
    /// [`std::thread::available_parallelism`].
    pub fn new(threads: usize) -> Self {
        let threads = if threads == 0 {
            std::thread::available_parallelism().map_or(1, |n| n.get())
        } else {
            threads
        };
        Self { threads }
    }

    /// Evaluate in order on the calling thread.
    pub fn serial() -> Self {
        Self { threads: 1 }
    }

    /// The number of threads.
    pub fn threads(&self) -> usize {
        self.threads
    }

    /// Run `work` on every item, with one `init()` state per thread. With
    /// one thread (or one item) the items are processed in order on the
    /// calling thread; otherwise the threads take the items in turn.
    fn run<T, W, I, F>(&self, items: Vec<T>, init: I, work: F)
    where
        T: Send,
        I: Fn() -> W + Sync,
        F: Fn(&mut W, T) + Sync,
    {
        let threads = self.threads.min(items.len());
        if threads <= 1 {
            let mut state = init();
            for item in items {
                work(&mut state, item);
            }
            return;
        }
        // pop() takes the items from the back: reverse to start at the front
        let queue = Mutex::new(items.into_iter().rev().collect::<Vec<T>>());
        std::thread::scope(|scope| {
            for _ in 0..threads {
                scope.spawn(|| {
                    let mut state = init();
                    loop {
                        let item = queue.lock().expect("batch queue").pop();
                        match item {
                            Some(item) => work(&mut state, item),
                            None => break,
                        }
                    }
                });
            }
        });
    }

    /// The number of consecutive waveforms per work item: about four items
    /// per thread, to balance waveforms of different cost.
    fn chunk_size(&self, n: usize) -> usize {
        if self.threads <= 1 {
            n.max(1)
        } else {
            n.div_ceil(4 * self.threads).max(1)
        }
    }

    /// Build a model for each binary with `build` (for example
    /// `IMRPhenomXPHM::from_binary`, or a closure around
    /// [`crate::choose_fd_model`]) and write its polarizations on `grid`
    /// ([`FrequencyDomainModel::polarizations_on_grid`]) into row `i` of
    /// `hplus` and `hcross`, `[i * grid.len()..(i + 1) * grid.len()]`.
    ///
    /// Returns the result of each binary; the rows of failed binaries are
    /// zero.
    ///
    /// # Panics
    ///
    /// If `hplus` or `hcross` do not have `binaries.len() * grid.len()`
    /// elements.
    pub fn polarizations_on_grid<M, B>(
        &self,
        binaries: &[Binary],
        build: B,
        grid: &FrequencyGrid,
        hplus: &mut [Complex<f64>],
        hcross: &mut [Complex<f64>],
    ) -> Vec<Result<()>>
    where
        M: FrequencyDomainModel,
        B: Fn(&Binary) -> Result<M> + Sync,
    {
        self.on_grid(binaries.len(), grid, hplus, hcross, |i, hp, hc| {
            build(&binaries[i])?.polarizations_on_grid(grid, hp, hc)
        })
    }

    /// The polarizations of prepared `models` on `grid`, written into row
    /// `i` of `hplus` and `hcross` as for [`Batch::polarizations_on_grid`].
    /// `models` can be concrete models or the boxed models of
    /// [`crate::choose_fd_model`].
    ///
    /// # Panics
    ///
    /// If `hplus` or `hcross` do not have `models.len() * grid.len()`
    /// elements.
    pub fn models_on_grid<M: FrequencyDomainModel>(
        &self,
        models: &[M],
        grid: &FrequencyGrid,
        hplus: &mut [Complex<f64>],
        hcross: &mut [Complex<f64>],
    ) -> Vec<Result<()>> {
        self.on_grid(models.len(), grid, hplus, hcross, |i, hp, hc| {
            models[i].polarizations_on_grid(grid, hp, hc)
        })
    }

    /// Row `i` of `hplus`, `hcross` from `evaluate(i, row_plus, row_cross)`
    /// for `i < n`.
    fn on_grid<E>(
        &self,
        n: usize,
        grid: &FrequencyGrid,
        hplus: &mut [Complex<f64>],
        hcross: &mut [Complex<f64>],
        evaluate: E,
    ) -> Vec<Result<()>>
    where
        E: Fn(usize, &mut [Complex<f64>], &mut [Complex<f64>]) -> Result<()> + Sync,
    {
        let len = grid.len();
        assert_eq!(hplus.len(), n * len, "hplus must hold one row per waveform");
        assert_eq!(
            hcross.len(),
            n * len,
            "hcross must hold one row per waveform"
        );
        let mut results: Vec<Result<()>> = (0..n).map(|_| Ok(())).collect();
        if len == 0 {
            // empty rows: evaluate for the errors only
            for (i, result) in results.iter_mut().enumerate() {
                *result = evaluate(i, &mut [], &mut []);
            }
            return results;
        }
        let chunk = self.chunk_size(n);
        let items: Vec<_> = hplus
            .chunks_mut(chunk * len)
            .zip(hcross.chunks_mut(chunk * len))
            .zip(results.chunks_mut(chunk))
            .enumerate()
            .map(|(k, ((hp, hc), results))| (k * chunk, hp, hc, results))
            .collect();
        self.run(
            items,
            || (),
            |_, (start, hp, hc, results)| {
                for (j, ((hp, hc), result)) in hp
                    .chunks_mut(len)
                    .zip(hc.chunks_mut(len))
                    .zip(results.iter_mut())
                    .enumerate()
                {
                    *result = evaluate(start + j, hp, hc);
                    if result.is_err() {
                        hp.fill(Complex::new(0.0, 0.0));
                        hc.fill(Complex::new(0.0, 0.0));
                    }
                }
            },
        );
        results
    }

    /// The polarizations at `frequencies` of the model built by `build` for
    /// each binary ([`FrequencyDomainModel::polarizations`]).
    pub fn polarizations<M, B>(
        &self,
        binaries: &[Binary],
        build: B,
        frequencies: &[f64],
    ) -> Vec<Result<Polarizations>>
    where
        M: FrequencyDomainModel,
        B: Fn(&Binary) -> Result<M> + Sync,
    {
        let n = binaries.len();
        let mut results: Vec<Result<Polarizations>> = (0..n)
            .map(|_| Err(Error::InvalidRequest("not evaluated")))
            .collect();
        let chunk = self.chunk_size(n);
        let items: Vec<_> = results
            .chunks_mut(chunk)
            .enumerate()
            .map(|(k, results)| (k * chunk, results))
            .collect();
        self.run(
            items,
            || (),
            |_, (start, results)| {
                for (j, result) in results.iter_mut().enumerate() {
                    *result = build(&binaries[start + j])
                        .and_then(|model| model.polarizations(frequencies));
                }
            },
        );
        results
    }

    /// The time-domain polarizations (sampling interval `delta_t`, from the
    /// frequency `f_min`) of the model built by `build` for each binary
    /// (for example `NRSur7dq4::from_binary`), written into `out[i]` and
    /// reusing its buffers. Each thread has its own
    /// [`TimeDomainModel::Workspace`].
    ///
    /// # Panics
    ///
    /// If `out` does not have the length of `binaries`.
    pub fn time_domain_polarizations<M, B>(
        &self,
        binaries: &[Binary],
        build: B,
        delta_t: f64,
        f_min: f64,
        out: &mut [TimeSeriesPolarizations],
    ) -> Vec<Result<()>>
    where
        M: TimeDomainModel,
        B: Fn(&Binary) -> Result<M> + Sync,
    {
        self.time_domain(binaries.len(), out, |workspace, i, out| {
            build(&binaries[i])?.polarizations_into(workspace, delta_t, f_min, out)
        })
    }

    /// The time-domain polarizations of prepared `models`, written into
    /// `out[i]` as for [`Batch::time_domain_polarizations`].
    ///
    /// # Panics
    ///
    /// If `out` does not have the length of `models`.
    pub fn time_domain_models<M: TimeDomainModel>(
        &self,
        models: &[M],
        delta_t: f64,
        f_min: f64,
        out: &mut [TimeSeriesPolarizations],
    ) -> Vec<Result<()>> {
        self.time_domain(models.len(), out, |workspace, i, out| {
            models[i].polarizations_into(workspace, delta_t, f_min, out)
        })
    }

    fn time_domain<W, E>(
        &self,
        n: usize,
        out: &mut [TimeSeriesPolarizations],
        evaluate: E,
    ) -> Vec<Result<()>>
    where
        W: Default,
        E: Fn(&mut W, usize, &mut TimeSeriesPolarizations) -> Result<()> + Sync,
    {
        assert_eq!(out.len(), n, "out must hold one waveform per binary");
        let mut results: Vec<Result<()>> = (0..n).map(|_| Ok(())).collect();
        let chunk = self.chunk_size(n);
        let items: Vec<_> = out
            .chunks_mut(chunk)
            .zip(results.chunks_mut(chunk))
            .enumerate()
            .map(|(k, (out, results))| (k * chunk, out, results))
            .collect();
        self.run(items, W::default, |workspace, (start, out, results)| {
            for (j, (out, result)) in out.iter_mut().zip(results.iter_mut()).enumerate() {
                *result = evaluate(workspace, start + j, out);
            }
        });
        results
    }
}

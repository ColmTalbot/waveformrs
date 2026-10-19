// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2024-2026 Colm Talbot. Parts of waveformrs are translated
// from LALSuite, GSL, glibc and fdlibm: see NOTICE.md for their authors.
//! Adapters of the crate's models for the bindings.

use std::cell::RefCell;
use std::panic::{catch_unwind, AssertUnwindSafe};

use num_complex::Complex64;
use waveformrs::conversions::FdWaveform;
use waveformrs::grid::FrequencyGrid;
use waveformrs::model::TimeSeriesModes;
use waveformrs::nrsur::{NRSurWorkspace, TimeSeriesPolarizations};
use waveformrs::{
    choose_td_model, Approximant, Binary, Error, FrequencyDomainModel, FrequencyModes, Result,
    TimeDomainModel, WaveformOptions,
};

/// A frequency-domain model whose grid evaluation is its evaluation at the
/// frequencies of the grid (`SimInspiralChooseFDWaveformSequence`, bit for
/// bit equal to `polarizations`), for batches with `exact=True`.
#[derive(Debug)]
pub struct Sequence(pub Box<dyn FrequencyDomainModel>);

impl FrequencyDomainModel for Sequence {
    fn from_binary(_binary: &Binary) -> Result<Self> {
        Err(Error::InvalidRequest("use choose_fd_model"))
    }

    fn approximant(&self) -> Approximant {
        self.0.approximant()
    }

    fn polarizations_into(
        &self,
        frequencies: &[f64],
        hplus: &mut [Complex64],
        hcross: &mut [Complex64],
    ) -> Result<()> {
        self.0.polarizations_into(frequencies, hplus, hcross)
    }

    fn polarizations_on_grid(
        &self,
        grid: &FrequencyGrid,
        hplus: &mut [Complex64],
        hcross: &mut [Complex64],
    ) -> Result<()> {
        self.0.polarizations_into(grid.frequencies(), hplus, hcross)
    }

    fn modes(&self, frequencies: &[f64]) -> Result<FrequencyModes> {
        self.0.modes(frequencies)
    }
}

/// `SimInspiralChooseTDWaveform` of a frequency-domain approximant
/// (`XLALSimInspiralTDFromFD`, [`FdWaveform::choose_td`]) as a
/// [`TimeDomainModel`].
#[derive(Debug)]
struct ChooseTd(FdWaveform);

impl TimeDomainModel for ChooseTd {
    type Workspace = NRSurWorkspace;

    fn from_binary(_binary: &Binary) -> Result<Self> {
        Err(Error::InvalidRequest("use td_model"))
    }

    fn approximant(&self) -> Approximant {
        self.0.approximant
    }

    fn polarizations_into(
        &self,
        _workspace: &mut NRSurWorkspace,
        delta_t: f64,
        f_min: f64,
        out: &mut TimeSeriesPolarizations,
    ) -> Result<()> {
        let h = self.0.choose_td(delta_t, f_min, self.0.binary.f_ref)?;
        out.epoch = h.hplus.epoch.seconds();
        out.delta_t = h.hplus.delta_t;
        out.hplus = h.hplus.data;
        out.hcross = h.hcross.data;
        Ok(())
    }

    fn modes_into(
        &self,
        _workspace: &mut NRSurWorkspace,
        _delta_t: f64,
        _f_min: f64,
        _out: &mut TimeSeriesModes,
    ) -> Result<()> {
        Err(Error::Unsupported {
            approximant: self.0.approximant,
            what: "time-domain modes",
        })
    }
}

/// Any approximant as a time-domain model: NRSur7dq4
/// ([`choose_td_model`]) or the `SimInspiralChooseTDWaveform` of a
/// frequency-domain approximant.
#[derive(Debug)]
pub struct AnyTd(Box<dyn TimeDomainModel<Workspace = NRSurWorkspace>>);

/// The time-domain model of `approximant` for `binary`.
pub fn td_model(
    approximant: Approximant,
    binary: &Binary,
    options: &WaveformOptions,
) -> Result<AnyTd> {
    if approximant.is_frequency_domain() {
        let waveform = FdWaveform::new(approximant, *binary, options.clone())?;
        Ok(AnyTd(Box::new(ChooseTd(waveform))))
    } else {
        Ok(AnyTd(choose_td_model(approximant, binary, options)?))
    }
}

impl TimeDomainModel for AnyTd {
    type Workspace = NRSurWorkspace;

    fn from_binary(_binary: &Binary) -> Result<Self> {
        Err(Error::InvalidRequest("use td_model"))
    }

    fn approximant(&self) -> Approximant {
        self.0.approximant()
    }

    fn polarizations_into(
        &self,
        workspace: &mut NRSurWorkspace,
        delta_t: f64,
        f_min: f64,
        out: &mut TimeSeriesPolarizations,
    ) -> Result<()> {
        self.0.polarizations_into(workspace, delta_t, f_min, out)
    }

    fn modes_into(
        &self,
        workspace: &mut NRSurWorkspace,
        delta_t: f64,
        f_min: f64,
        out: &mut TimeSeriesModes,
    ) -> Result<()> {
        self.0.modes_into(workspace, delta_t, f_min, out)
    }
}

thread_local! {
    /// The NRSur7dq4 workspace of each thread for single waveforms: no
    /// allocations after the first call, and the surrogate is not evaluated
    /// again when only extrinsic parameters change.
    static WORKSPACE: RefCell<NRSurWorkspace> = RefCell::new(NRSurWorkspace::new());
}

/// Run `f` with the workspace of this thread.
pub fn with_workspace<T>(f: impl FnOnce(&mut NRSurWorkspace) -> T) -> T {
    WORKSPACE.with_borrow_mut(f)
}

/// The error of a binary whose evaluation panicked (raised as `ModelError`).
pub const PANICKED: &str =
    "the model panicked (a bug of waveformrs; the panic message is on stderr)";

/// `f()`, with a panic turned into an error: one binary of a batch cannot
/// abort the others.
pub fn guard<T>(f: impl FnOnce() -> Result<T>) -> Result<T> {
    catch_unwind(AssertUnwindSafe(f)).unwrap_or(Err(Error::InvalidRequest(PANICKED)))
}

/// A model whose evaluation cannot panic (see [`guard`]), for batches.
#[derive(Debug)]
pub struct Guarded<M>(pub M);

impl<M: FrequencyDomainModel> FrequencyDomainModel for Guarded<M> {
    fn from_binary(binary: &Binary) -> Result<Self> {
        guard(|| M::from_binary(binary)).map(Self)
    }

    fn approximant(&self) -> Approximant {
        self.0.approximant()
    }

    fn polarizations_into(
        &self,
        frequencies: &[f64],
        hplus: &mut [Complex64],
        hcross: &mut [Complex64],
    ) -> Result<()> {
        guard(|| self.0.polarizations_into(frequencies, hplus, hcross))
    }

    fn polarizations_on_grid(
        &self,
        grid: &FrequencyGrid,
        hplus: &mut [Complex64],
        hcross: &mut [Complex64],
    ) -> Result<()> {
        guard(|| self.0.polarizations_on_grid(grid, hplus, hcross))
    }

    fn modes(&self, frequencies: &[f64]) -> Result<FrequencyModes> {
        guard(|| self.0.modes(frequencies))
    }
}

impl<M: TimeDomainModel<Workspace = NRSurWorkspace>> TimeDomainModel for Guarded<M> {
    type Workspace = NRSurWorkspace;

    fn from_binary(binary: &Binary) -> Result<Self> {
        guard(|| M::from_binary(binary)).map(Self)
    }

    fn approximant(&self) -> Approximant {
        self.0.approximant()
    }

    fn polarizations_into(
        &self,
        workspace: &mut NRSurWorkspace,
        delta_t: f64,
        f_min: f64,
        out: &mut TimeSeriesPolarizations,
    ) -> Result<()> {
        guard(|| self.0.polarizations_into(workspace, delta_t, f_min, out))
    }

    fn modes_into(
        &self,
        workspace: &mut NRSurWorkspace,
        delta_t: f64,
        f_min: f64,
        out: &mut TimeSeriesModes,
    ) -> Result<()> {
        guard(|| self.0.modes_into(workspace, delta_t, f_min, out))
    }
}

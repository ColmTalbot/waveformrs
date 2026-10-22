// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2024-2026 Colm Talbot. Parts of waveformrs are translated
// from LALSuite, GSL and fdlibm: see NOTICE.md for their authors.
//! Conditioned waveforms and conversions between the time and frequency
//! domains, following `XLALSimInspiralTD` and `XLALSimInspiralFD` of
//! LALSimulation 6.2.0 (`LALSimInspiralGeneratorConditioning.c`,
//! `LALSimInspiral.c`; unchanged in LALSimulation master as of 2025).
//!
//! | function | LAL | for |
//! |----------|-----|-----|
//! | [`simulate_td`] | `SimInspiralTD` | any [`Approximant`](crate::Approximant) and [`Binary`](crate::Binary) |
//! | [`simulate_fd`] | `SimInspiralFD` | any [`Approximant`](crate::Approximant) and [`Binary`](crate::Binary) |
//! | [`td_from_fd`] | `XLALSimInspiralTDFromFD` (= `generate_conditioned_td_waveform_from_fd`) | an [`FdGenerator`] |
//! | [`fd_from_td`] | `SimInspiralFD` (`generate_conditioned_fd_waveform_from_td`) | a [`TdGenerator`] |
//! | [`conditioned_fd`] | `SimInspiralFD` (`generate_conditioned_fd_waveform_from_fd`) | an [`FdGenerator`] |
//! | [`conditioned_td`] | `SimInspiralTD` (`generate_conditioned_td_waveform_from_td` and its `_fallback`) | a [`TdGenerator`] |
//! | [`FdWaveform::time_domain`], [`FdWaveform::conditioned`] | `SimInspiralTD`, `SimInspiralFD` | the frequency-domain models of this crate |
//! | [`nrsur7dq4_fd`], [`nrsur7dq4_td`] | `SimInspiralFD`, `SimInspiralTD` | NRSur7dq4 |
//!
//! [`td_from_fd`], [`fd_from_td`], [`simulate_td`] and [`simulate_fd`] need
//! the Fourier transforms of `crate::fft` (feature `fft`, on by default).
//!
//! A model enters the generic functions through the small traits
//! [`FdGenerator`] and [`TdGenerator`]: the masses and aligned spins that
//! LAL's duration bounds use, a few flags that LAL keys on the approximant,
//! and a generator that reproduces `SimInspiralChooseFDWaveform` on the
//! uniform grid `k delta_f` (respectively `SimInspiralChooseTDWaveform`)
//! for the sampling, starting frequency and reference frequency the
//! conditioning asks for (unlike the crate's
//! [`FrequencyDomainModel`](crate::FrequencyDomainModel), whose parameters
//! are fixed, since the conditioning changes `f_min` and, through it, the
//! default reference frequency). Closures can be used through
//! [`FdGeneratorFn`] and [`TdGeneratorFn`]; [`FdWaveform`] (an
//! [`Approximant`](crate::Approximant), a [`Binary`](crate::Binary) and
//! [`WaveformOptions`](crate::WaveformOptions)) adapts the
//! frequency-domain models of this crate (IMRPhenomD, IMRPhenomXAS/XHM/XP/
//! XPHM, TaylorF2), and [`TdModel`] any
//! [`TimeDomainModel`](crate::TimeDomainModel) with its binary (NRSur7dq4,
//! including boxed models from [`crate::choose_td_model`]).
//!
//! # From the frequency to the time domain
//!
//! For a sampling interval `delta_t`, starting frequency `f_min` and
//! reference frequency `f_ref`, [`td_from_fd`]
//!
//! 1. replaces `f_ref = 0` by `f_min` for the models that interpret it so
//!    (`FIX_REFERENCE_FREQUENCY`: IMRPhenomXP/XPHM here) and lowers `f_min`
//!    to the lowest Kerr ISCO frequency if it is above it;
//! 2. bounds the chirp time `t_chirp` from `f_min`, the merger and ringdown
//!    time `t_merge` and sets `t_extra = 3 / f_min`;
//! 3. generates the conditioned frequency-domain waveform ([`conditioned_fd`]
//!    with `delta_f = 0` and `f_max = 0.5 / delta_t`): the model starts at
//!    `f_start`, where the Newtonian chirp lasts `1.1 t_chirp`, on a grid
//!    with `delta_f = 1 / (N delta_t)` and `N` the power of two above
//!    `(t_chirp(f_start) + t_merge + 2 t_extra) / delta_t` (twice that if
//!    it is a power of two already); it is tapered with a Hann window from
//!    `f_start` to `f_min`, zeroed at Nyquist, and shifted by `t_merge`
//!    (rounded to samples; zero for inspiral-only models) backwards in time
//!    with the epoch moved accordingly;
//! 4. shifts it by a further `t_extra` (rounded to samples), transforms to
//!    the time domain (`XLALREAL8FreqTimeFFT`), high-pass filters at the
//!    original `f_min` (8th-order Butterworth, attenuation 0.99, forwards
//!    and backwards), keeps the `round((t_chirp(f_start') + t_merge) /
//!    delta_t)` samples ending `t_extra` before the end (`f_start'` from
//!    `1.1 t_chirp` again), and tapers one cycle at `f_min` at the start
//!    and one cycle at the Schwarzschild ISCO frequency at the end
//!    (`XLALSimInspiralTDConditionStage2`).
//!
//! The epoch is that of the model's frequency series (`-1 / delta_f` for
//! all the models here) plus the two shifts and the samples cut at the
//! start, each rounded to nanoseconds as `XLALGPSAdd` does.
//!
//! This is `SimInspiralTD` only for frequency-domain approximants without
//! a time-domain generator (TaylorF2 here). In LALSimulation 6.2.0,
//! IMRPhenomD, IMRPhenomXAS, IMRPhenomXHM, IMRPhenomXP, IMRPhenomXPHM and
//! IMRPhenomXPNR do have one: `SimInspiralChooseTDWaveform` is `XLALSimInspiralTDFromFD` and,
//! for IMRPhenomD/XAS/XHM, moves the epoch so that the peak of `h_+^2 +
//! h_x^2` is at `t = 0` (IMRPhenomD is transformed face-on and the
//! inclination applied in the time domain). `SimInspiralTD` then conditions
//! that as a time-domain waveform ([`conditioned_td`]): it is generated
//! from a frequency at which the chirp is longer by `0.1 t_chirp + t_merge
//! + t_extra`, tapered over the extra time, high-pass filtered at `f_min`
//! ([`crate::conditioning::td_condition_stage_1`]) and tapered at both ends
//! again. [`FdWaveform::time_domain`] does all this.
//!
//! # From the time to the frequency domain (`SimInspiralFD`)
//!
//! For a frequency spacing `delta_f` (zero for automatic), `f_min`, `f_max`
//! and `f_ref`, [`fd_from_td`] rounds `f_max / delta_f` up to a power of
//! two if it is not one (`f_nyquist`), generates the conditioned
//! time-domain waveform ([`conditioned_td`]) with `delta_t = 0.5 /
//! f_nyquist`, zero pads (or truncates) it at the start to `N = 2
//! f_nyquist / delta_f` samples, or to the power of two above its length
//! for `delta_f = 0`, and Fourier transforms it (`XLALREAL8TimeFreqFFT`:
//! scaled by `delta_t`, `delta_f = 1 / (N delta_t)`, the epoch of the
//! padded time series). Models that start at the requested frequency
//! (NRSur7dq4, which allows `f_min = 0`) are only tapered at the start with
//! [`crate::conditioning::wave_taper`]; others are generated from a lower
//! frequency and conditioned with
//! [`crate::conditioning::td_condition_stage_1`] and
//! [`crate::conditioning::td_condition_stage_2`].
//!
//! Note that `SimInspiralChooseFDWaveform` itself does not convert
//! time-domain models in LALSimulation 6.2.0 (it fails with "generator
//! does not provide a method to generate frequency-domain waveforms"); only
//! `SimInspiralFD` does.
//!
//! # Accuracy
//!
//! The bounds, tapers, filter, epochs and lengths reproduce LAL exactly
//! (the Butterworth filter to rounding for identical input). The Fourier
//! transforms are not FFTW's and differ from it by a few `1e-16` of the
//! largest value. Transforms to the frequency domain (NRSur7dq4) therefore
//! agree with LAL to `~3e-16`. Transforms to the time domain are high-pass
//! filtered afterwards, and the direct-form IIR sections of LAL's filter,
//! with poles close to `z = 1`, amplify last-bit differences of their input
//! by `~1e4`: these agree with LAL to `~1e-12` of the peak (and the
//! frequency-domain models contribute their own rounding differences, up
//! to `~3e-11` for the long TaylorF2 inspirals of binary neutron stars).
#![cfg_attr(
    feature = "fft",
    doc = "\n",
    doc = include_str!("examples.md")
)]

#[cfg(feature = "fft")]
use crate::conditioning::high_pass_pair;
use crate::conditioning::{
    chirp_start_frequency_bound, chirp_time_bound, final_black_hole_spin_bound,
    kerr_isco_frequency, merge_time_bound, ringdown_time_bound_kg, schwarzschild_isco_frequency,
    taper_frequency_band, td_condition_stage_1, td_condition_stage_2, time_shift_pair, wave_taper,
    WaveTaper,
};
use crate::constants::MSUN_SI;
#[cfg(feature = "fft")]
use crate::fft::{frequency_to_time, time_to_frequency};
pub use crate::series::{
    FrequencyDomainPolarizations, FrequencySeries, GpsTime, TimeDomainPolarizations, TimeSeries,
};
use crate::spintaylor::fp::fma;

mod models;

#[cfg(feature = "fft")]
pub use models::{nrsur7dq4_fd, simulate_fd, simulate_td};
pub use models::{nrsur7dq4_td, FdWaveform, TdModel};

/// Fraction of the chirp time added before `f_min` for tapering.
const EXTRA_TIME_FRACTION: f64 = 0.1;
/// Cycles at `f_min` added as padding.
const EXTRA_CYCLES: f64 = 3.0;

/// The parameters of a binary that LAL's conditioning uses: the component
/// masses (solar masses) and the spin components along the orbital angular
/// momentum.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ConditioningParameters {
    pub mass_1: f64,
    pub mass_2: f64,
    pub spin_1z: f64,
    pub spin_2z: f64,
}

impl ConditioningParameters {
    /// `t_merge`: [`merge_time_bound`] plus [`ringdown_time_bound_kg`] for
    /// the [`final_black_hole_spin_bound`].
    fn merge_and_ringdown_time(&self) -> f64 {
        let s = final_black_hole_spin_bound(self.spin_1z, self.spin_2z);
        merge_time_bound(self.mass_1, self.mass_2)
            + ringdown_time_bound_kg(self.mass_1 * MSUN_SI + self.mass_2 * MSUN_SI, s)
    }

    fn chirp_time(&self, f_start: f64) -> f64 {
        chirp_time_bound(
            f_start,
            self.mass_1,
            self.mass_2,
            self.spin_1z,
            self.spin_2z,
        )
    }

    fn chirp_start_frequency(&self, t_chirp: f64) -> f64 {
        chirp_start_frequency_bound(t_chirp, self.mass_1, self.mass_2)
    }
}

/// What the conditioning asks a frequency-domain model for: the arguments
/// of `SimInspiralChooseFDWaveform` that it changes.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FdRequest {
    /// Frequency spacing (Hz) of the grid `k delta_f`, `k = 0, 1, ...`.
    pub delta_f: f64,
    /// Starting frequency (Hz): the model is zero below it.
    pub f_min: f64,
    /// Maximum frequency (Hz), LAL's `f_max`.
    pub f_max: f64,
    /// Reference frequency (Hz), after `FIX_REFERENCE_FREQUENCY`.
    pub f_ref: f64,
}

/// What the conditioning asks a time-domain model for: the arguments of
/// `SimInspiralChooseTDWaveform` that it changes.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TdRequest {
    /// Sampling interval (s).
    pub delta_t: f64,
    /// Starting frequency (Hz).
    pub f_min: f64,
    /// Reference frequency (Hz), after `FIX_REFERENCE_FREQUENCY`.
    pub f_ref: f64,
}

/// A frequency-domain model as LAL's conditioning sees it.
pub trait FdGenerator {
    /// The error of the model.
    type Error;

    /// The masses and aligned spins used by the duration bounds.
    fn conditioning_parameters(&self) -> ConditioningParameters;

    /// Whether the model interprets `f_ref = 0` as `f_min` (the
    /// approximants of LAL's `FIX_REFERENCE_FREQUENCY`, e.g. IMRPhenomXP
    /// and IMRPhenomXPHM): the conditioning then passes the requested
    /// `f_min` rather than its lower starting frequency.
    fn zero_reference_frequency_is_f_min(&self) -> bool {
        false
    }

    /// Whether the model is inspiral only (TaylorF2): no merger-ringdown
    /// time is allowed for in the frequency-domain conditioning.
    fn is_inspiral_only(&self) -> bool {
        false
    }

    /// `(h_+, h_x)` on the grid `k delta_f`, `k = 0 .. n - 1`, as
    /// `SimInspiralChooseFDWaveform(..., delta_f, f_min, f_max, f_ref, ...)`
    /// returns them, including the length `n` and the epoch chosen by the
    /// model.
    fn generate_fd(&self, request: &FdRequest)
        -> Result<FrequencyDomainPolarizations, Self::Error>;
}

/// A time-domain model as LAL's conditioning sees it.
pub trait TdGenerator {
    /// The error of the model.
    type Error;

    /// The masses and aligned spins used by the duration bounds.
    fn conditioning_parameters(&self) -> ConditioningParameters;

    /// Whether the model interprets `f_ref = 0` as `f_min`
    /// (`FIX_REFERENCE_FREQUENCY`).
    fn zero_reference_frequency_is_f_min(&self) -> bool {
        false
    }

    /// Whether the waveform must start at the requested frequency: the
    /// models whose spins are defined at `f_min`, or which allow `f_min =
    /// 0` (NRSur7dq4). These are only tapered at the start
    /// (`generate_conditioned_td_waveform_from_td_fallback`) instead of
    /// being generated from a lower frequency.
    fn fixed_start_frequency(&self) -> bool {
        false
    }

    /// `(h_+, h_x)` as `SimInspiralChooseTDWaveform(..., delta_t, f_min,
    /// f_ref, ...)` returns them.
    fn generate_td(&self, request: &TdRequest) -> Result<TimeDomainPolarizations, Self::Error>;
}

/// Errors of the conversions.
#[derive(Clone, Debug, PartialEq)]
pub enum ConversionError<E> {
    /// The model failed.
    Model(E),
    /// The sampling arguments are invalid.
    InvalidArgument(&'static str),
}

impl<E: std::fmt::Display> std::fmt::Display for ConversionError<E> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Model(err) => write!(f, "{err}"),
            Self::InvalidArgument(why) => write!(f, "{why}"),
        }
    }
}

impl<E: std::error::Error + 'static> std::error::Error for ConversionError<E> {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Model(err) => Some(err),
            Self::InvalidArgument(_) => None,
        }
    }
}

/// The errors of the models of this crate pass through; invalid sampling
/// arguments become [`crate::Error::InvalidRequest`].
impl From<ConversionError<crate::Error>> for crate::Error {
    fn from(err: ConversionError<crate::Error>) -> Self {
        match err {
            ConversionError::Model(err) => err,
            ConversionError::InvalidArgument(why) => crate::Error::InvalidRequest(why),
        }
    }
}

/// `FIX_REFERENCE_FREQUENCY`.
fn fix_reference_frequency(f_ref: f64, f_min: f64, zero_is_f_min: bool) -> f64 {
    if f_ref == 0.0 && zero_is_f_min {
        f_min
    } else {
        f_ref
    }
}

/// The exponent `e` of C's `frexp(x, &e)`, `x = m 2^e` with `m` in
/// `[0.5, 1)`, for non-negative finite `x` (zero for zero).
fn frexp_exponent(x: f64) -> i32 {
    if x == 0.0 {
        return 0;
    }
    let (mantissa_scale, x) = if x < f64::MIN_POSITIVE {
        (-54, x * 2f64.powi(54))
    } else {
        (0, x)
    };
    ((x.to_bits() >> 52) & 0x7ff) as i32 - 1022 + mantissa_scale
}

/// `ldexp(1.0, frexp_exponent(x))`: the smallest power of two above `x`
/// (twice `x` if `x` is a power of two), as LAL rounds lengths.
fn next_power_of_two_above(x: f64) -> f64 {
    2f64.powi(frexp_exponent(x))
}

/// The Nyquist frequency of `SimInspiralFD`: `f_max`, or `f_max / delta_f`
/// rounded up to a power of two times `delta_f` if it is not one.
fn nyquist_frequency(f_max: f64, delta_f: f64) -> f64 {
    if delta_f != 0.0 {
        let n = (f_max / delta_f).round() as u64;
        if n & n.wrapping_sub(1) != 0 {
            return next_power_of_two_above(n as f64) * delta_f;
        }
    }
    f_max
}

/// Whether `x` is positive (false for NaN).
pub(crate) fn positive(x: f64) -> bool {
    x > 0.0
}

/// Whether `x` is positive or zero (false for NaN).
pub(crate) fn non_negative(x: f64) -> bool {
    x >= 0.0
}

/// `(int)(a - b)` for C's `size_t` operands: the difference wrapped to 32
/// bits.
#[cfg(feature = "fft")]
fn c_int_difference(a: usize, b: usize) -> i64 {
    a.wrapping_sub(b) as u32 as i32 as i64
}

/// `SimInspiralFD` for a frequency-domain model
/// (`generate_conditioned_fd_waveform_from_fd`): the model on the grid `k
/// delta_f` from a lower starting frequency, tapered between that and
/// `f_min`, with zero at Nyquist and shifted backwards in time by the
/// merger-ringdown time bound so that it can be transformed to the time
/// domain without wrapping around.
///
/// `delta_f = 0` chooses the spacing for which the transform of the
/// waveform (with padding) fits into a power-of-two number of samples at
/// `delta_t = 0.5 / f_max`. A `delta_f` that is too large (a waveform
/// longer than `1 / delta_f`) is used anyway, as LAL does with a warning.
pub fn conditioned_fd<M: FdGenerator + ?Sized>(
    model: &M,
    delta_f: f64,
    f_min: f64,
    f_max: f64,
    f_ref: f64,
) -> Result<FrequencyDomainPolarizations, ConversionError<M::Error>> {
    if !non_negative(delta_f) || !positive(f_min) || !positive(f_max) {
        return Err(ConversionError::InvalidArgument(
            "delta_f must be non-negative, f_min and f_max positive",
        ));
    }
    let p = model.conditioning_parameters();
    let f_ref = fix_reference_frequency(f_ref, f_min, model.zero_reference_frequency_is_f_min());

    let f_nyquist = nyquist_frequency(f_max, delta_f);
    let delta_t = 0.5 / f_nyquist;

    let f_min = f_min.min(kerr_isco_frequency(p.mass_1, p.mass_2));
    let t_chirp = p.chirp_time(f_min);
    let t_merge = if model.is_inspiral_only() {
        0.0
    } else {
        p.merge_and_ringdown_time()
    };
    let t_extra = EXTRA_CYCLES / f_min;
    let f_start = p.chirp_start_frequency((1.0 + EXTRA_TIME_FRACTION) * t_chirp);
    let t_chirp = p.chirp_time(f_start);

    let chirplen = (fma(2.0, t_extra, t_chirp + t_merge) / delta_t).round();
    let chirplen = next_power_of_two_above(chirplen);
    let delta_f = if delta_f == 0.0 {
        1.0 / (chirplen * delta_t)
    } else {
        delta_f
    };

    let request = FdRequest {
        delta_f,
        f_min: f_start,
        f_max,
        f_ref,
    };
    let mut h = model
        .generate_fd(&request)
        .map_err(ConversionError::Model)?;
    taper_frequency_band(&mut h.hplus, f_start, f_min);
    taper_frequency_band(&mut h.hcross, f_start, f_min);
    let t_shift = (t_merge / delta_t).round() * delta_t;
    time_shift_pair(&mut h.hplus, &mut h.hcross, delta_f, t_shift);
    Ok(h)
}

/// `SimInspiralTD` for a frequency-domain model
/// (`generate_conditioned_td_waveform_from_fd`, see the [module
/// documentation](self)).
#[cfg(feature = "fft")]
pub fn td_from_fd<M: FdGenerator + ?Sized>(
    model: &M,
    delta_t: f64,
    f_min: f64,
    f_ref: f64,
) -> Result<TimeDomainPolarizations, ConversionError<M::Error>> {
    if !positive(delta_t) || !positive(f_min) {
        return Err(ConversionError::InvalidArgument(
            "delta_t and f_min must be positive",
        ));
    }
    let p = model.conditioning_parameters();
    let original_f_min = f_min;
    let f_max = 0.5 / delta_t;
    let f_ref = fix_reference_frequency(f_ref, f_min, model.zero_reference_frequency_is_f_min());

    let f_min = f_min.min(kerr_isco_frequency(p.mass_1, p.mass_2));
    let t_chirp = p.chirp_time(f_min);
    let t_merge = p.merge_and_ringdown_time();
    let t_extra = EXTRA_CYCLES / f_min;

    let mut fd = conditioned_fd(model, 0.0, original_f_min, f_max, f_ref)?;

    let t_shift = (t_extra / delta_t).round() * delta_t;
    let delta_f = fd.hplus.delta_f;
    time_shift_pair(&mut fd.hplus, &mut fd.hcross, delta_f, t_shift);
    let mut td = [fd.hplus, fd.hcross].map(|h| {
        let chirplen = 2 * (h.data.len() - 1);
        frequency_to_time(&h, chirplen)
    });
    {
        let [hplus, hcross] = &mut td;
        // LAL ignores failures of the filter here
        let _ = high_pass_pair(hplus, hcross, original_f_min, 0.99, 8);
    }

    let f_start = p.chirp_start_frequency((1.0 + EXTRA_TIME_FRACTION) * t_chirp);
    let t_chirp = p.chirp_time(f_start);
    let chirplen = ((t_chirp + t_merge) / delta_t).round() as usize;
    for series in &mut td {
        let end = (series.data.len() as f64 - (t_shift / delta_t).round()) as usize;
        series.resize(c_int_difference(end, chirplen), chirplen);
    }
    let [mut hplus, mut hcross] = td;
    let f_isco = schwarzschild_isco_frequency(p.mass_1, p.mass_2);
    td_condition_stage_2(&mut hplus, &mut hcross, f_min, f_isco);
    Ok(TimeDomainPolarizations { hplus, hcross })
}

/// `SimInspiralTD` for a time-domain model: for models with a
/// [fixed start frequency](TdGenerator::fixed_start_frequency) the
/// waveform at `f_min` tapered at the start with
/// [`wave_taper`] (`generate_conditioned_td_waveform_from_td_fallback`);
/// otherwise the waveform from a lower frequency (where the chirp is
/// longer by `0.1 t_chirp + t_merge + 3 / f_min`), conditioned with
/// [`td_condition_stage_1`] and [`td_condition_stage_2`]
/// (`generate_conditioned_td_waveform_from_td`).
pub fn conditioned_td<M: TdGenerator + ?Sized>(
    model: &M,
    delta_t: f64,
    f_min: f64,
    f_ref: f64,
) -> Result<TimeDomainPolarizations, ConversionError<M::Error>> {
    if !positive(delta_t) || !non_negative(f_min) {
        return Err(ConversionError::InvalidArgument(
            "delta_t must be positive and f_min non-negative",
        ));
    }
    if model.fixed_start_frequency() {
        let request = TdRequest {
            delta_t,
            f_min,
            f_ref,
        };
        let mut h = model
            .generate_td(&request)
            .map_err(ConversionError::Model)?;
        wave_taper(&mut h.hplus.data, WaveTaper::Start);
        wave_taper(&mut h.hcross.data, WaveTaper::Start);
        return Ok(h);
    }
    if f_min == 0.0 {
        return Err(ConversionError::InvalidArgument(
            "f_min must be positive for this model",
        ));
    }
    let p = model.conditioning_parameters();
    let original_f_min = f_min;
    let f_ref = fix_reference_frequency(f_ref, f_min, model.zero_reference_frequency_is_f_min());
    let f_min = f_min.min(kerr_isco_frequency(p.mass_1, p.mass_2));
    let t_chirp = p.chirp_time(f_min);
    let t_merge = p.merge_and_ringdown_time();
    let t_extra = EXTRA_CYCLES / f_min;
    let f_start =
        p.chirp_start_frequency(fma(1.0 + EXTRA_TIME_FRACTION, t_chirp, t_merge) + t_extra);
    let request = TdRequest {
        delta_t,
        f_min: f_start,
        f_ref,
    };
    let mut h = model
        .generate_td(&request)
        .map_err(ConversionError::Model)?;
    td_condition_stage_1(
        &mut h.hplus,
        &mut h.hcross,
        fma(EXTRA_TIME_FRACTION, t_chirp, t_extra),
        original_f_min,
    );
    let f_isco = schwarzschild_isco_frequency(p.mass_1, p.mass_2);
    td_condition_stage_2(&mut h.hplus, &mut h.hcross, f_min, f_isco);
    Ok(h)
}

/// `SimInspiralFD` for a time-domain model
/// (`generate_conditioned_fd_waveform_from_td`, see the [module
/// documentation](self)). `delta_f = 0` chooses the spacing from the
/// length of the waveform.
#[cfg(feature = "fft")]
pub fn fd_from_td<M: TdGenerator + ?Sized>(
    model: &M,
    delta_f: f64,
    f_min: f64,
    f_max: f64,
    f_ref: f64,
) -> Result<FrequencyDomainPolarizations, ConversionError<M::Error>> {
    if !non_negative(delta_f) || !positive(f_max) {
        return Err(ConversionError::InvalidArgument(
            "delta_f must be non-negative and f_max positive",
        ));
    }
    let f_ref = fix_reference_frequency(f_ref, f_min, model.zero_reference_frequency_is_f_min());
    let f_nyquist = nyquist_frequency(f_max, delta_f);
    let delta_t = 0.5 / f_nyquist;
    let h = conditioned_td(model, delta_t, f_min, f_ref)?;
    let length = h.hplus.data.len();
    let chirplen = if delta_f == 0.0 {
        next_power_of_two_above(length as f64)
    } else {
        2.0 * f_nyquist / delta_f
    } as usize;
    if chirplen == 0 {
        return Err(ConversionError::InvalidArgument(
            "the frequency series would be empty",
        ));
    }
    let [hplus, hcross] = [h.hplus, h.hcross].map(|mut series| {
        series.resize(c_int_difference(series.data.len(), chirplen), chirplen);
        time_to_frequency(&series)
    });
    Ok(FrequencyDomainPolarizations { hplus, hcross })
}

/// A [`FdGenerator`] from a closure `Fn(&FdRequest) ->
/// Result<FrequencyDomainPolarizations, E>`.
#[derive(Clone, Copy, Debug)]
pub struct FdGeneratorFn<F> {
    /// The masses and aligned spins for the duration bounds.
    pub parameters: ConditioningParameters,
    /// See [`FdGenerator::zero_reference_frequency_is_f_min`].
    pub zero_reference_frequency_is_f_min: bool,
    /// See [`FdGenerator::is_inspiral_only`].
    pub inspiral_only: bool,
    /// The generator.
    pub generate: F,
}

impl<F, E> FdGenerator for FdGeneratorFn<F>
where
    F: Fn(&FdRequest) -> Result<FrequencyDomainPolarizations, E>,
{
    type Error = E;

    fn conditioning_parameters(&self) -> ConditioningParameters {
        self.parameters
    }

    fn zero_reference_frequency_is_f_min(&self) -> bool {
        self.zero_reference_frequency_is_f_min
    }

    fn is_inspiral_only(&self) -> bool {
        self.inspiral_only
    }

    fn generate_fd(&self, request: &FdRequest) -> Result<FrequencyDomainPolarizations, E> {
        (self.generate)(request)
    }
}

/// A [`TdGenerator`] from a closure `Fn(&TdRequest) ->
/// Result<TimeDomainPolarizations, E>`.
#[derive(Clone, Copy, Debug)]
pub struct TdGeneratorFn<F> {
    /// The masses and aligned spins for the duration bounds.
    pub parameters: ConditioningParameters,
    /// See [`TdGenerator::zero_reference_frequency_is_f_min`].
    pub zero_reference_frequency_is_f_min: bool,
    /// See [`TdGenerator::fixed_start_frequency`].
    pub fixed_start_frequency: bool,
    /// The generator.
    pub generate: F,
}

impl<F, E> TdGenerator for TdGeneratorFn<F>
where
    F: Fn(&TdRequest) -> Result<TimeDomainPolarizations, E>,
{
    type Error = E;

    fn conditioning_parameters(&self) -> ConditioningParameters {
        self.parameters
    }

    fn zero_reference_frequency_is_f_min(&self) -> bool {
        self.zero_reference_frequency_is_f_min
    }

    fn fixed_start_frequency(&self) -> bool {
        self.fixed_start_frequency
    }

    fn generate_td(&self, request: &TdRequest) -> Result<TimeDomainPolarizations, E> {
        (self.generate)(request)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn frexp_like_c() {
        assert_eq!(frexp_exponent(1.0), 1);
        assert_eq!(frexp_exponent(0.5), 0);
        assert_eq!(frexp_exponent(3.0), 2);
        assert_eq!(frexp_exponent(4.0), 3);
        assert_eq!(frexp_exponent(4095.0), 12);
        assert_eq!(next_power_of_two_above(4096.0), 8192.0);
        assert_eq!(next_power_of_two_above(4097.0), 8192.0);
        assert_eq!(next_power_of_two_above(0.0), 1.0);
        assert_eq!(nyquist_frequency(2048.0, 0.25), 2048.0);
        assert_eq!(nyquist_frequency(1000.0, 0.25), 1024.0);
        assert_eq!(nyquist_frequency(1000.0, 0.0), 1000.0);
    }

    #[cfg(feature = "fft")]
    #[test]
    fn int_difference_wraps() {
        assert_eq!(c_int_difference(10, 3), 7);
        assert_eq!(c_int_difference(3, 10), -7);
    }
}

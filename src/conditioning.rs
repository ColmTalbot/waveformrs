//! LALSimulation's generic waveform conditioning.
//!
//! The bounds on the duration of a signal, the tapers, the Butterworth
//! high-pass filter and the two stages of time-domain conditioning used by
//! `XLALSimInspiralTD` and `XLALSimInspiralFD` (LALSimInspiral.c and
//! LALSimInspiralGeneratorConditioning.c of LALSimulation 6.2.0, unchanged
//! in later versions up to 2025), ported operation by operation:
//!
//! | function | LAL |
//! |----------|-----|
//! | [`chirp_time_bound`] | `XLALSimInspiralChirpTimeBound` |
//! | [`merge_time_bound`] | `XLALSimInspiralMergeTimeBound` |
//! | [`ringdown_time_bound`] | `XLALSimInspiralRingdownTimeBound` |
//! | [`final_black_hole_spin_bound`] | `XLALSimInspiralFinalBlackHoleSpinBound` |
//! | [`chirp_start_frequency_bound`] | `XLALSimInspiralChirpStartFrequencyBound` |
//! | [`kerr_isco_frequency`], [`schwarzschild_isco_frequency`] | the `fisco` of the conditioning routines |
//! | [`wave_taper`] | `XLALSimInspiralREAL8WaveTaper` (Planck-like taper) |
//! | [`hann_taper_start`], [`hann_taper_end`] | the cosine (Tukey) tapers of the conditioning stages |
//! | [`high_pass`], [`high_pass_pair`] | `XLALHighPassREAL8TimeSeries` (LAL's `BandPassTimeSeries`) |
//! | [`td_condition_stage_1`], [`td_condition_stage_2`] | `XLALSimInspiralTDConditionStage1/2` |
//! | [`taper_frequency_band`] | the frequency-domain taper of `XLALSimInspiralFD` |
//! | [`time_shift`], [`time_shift_pair`] | the phase winding `exp(2 pi i f t_shift)` of the FD conditioning |
//! | [`peak_index`] | the peak alignment of `XLALSimInspiralChooseTDWaveform` for IMRPhenomD/XAS/XHM |
//!
//! Masses are in solar masses; they are converted to kilograms as
//! `m * MSUN_SI` for each body (as callers of LAL do) before LAL's
//! expressions are evaluated, with the fused multiply-adds of clang's arm64
//! build (see [`crate::spintaylor::fp`]). Everything here except the filter
//! coefficients is bit for bit identical to LAL on Apple silicon; the
//! filter agrees to rounding.

use std::f64::consts::PI;

use num_complex::Complex;

use crate::cmath::LibmCbrt;
use crate::constants::{C_SI, G_SI, MRSUN_SI, MSUN_SI, MTSUN_SI};
use crate::series::{FrequencySeries, GpsTime, TimeSeries};
use crate::spintaylor::fp::{fma, fms};

/// `LAL_TWOPI`.
const TWO_PI: f64 = 2.0 * PI;

/// The masses in kilograms as LAL receives them.
#[inline]
fn kilograms(mass_1: f64, mass_2: f64) -> (f64, f64) {
    (mass_1 * MSUN_SI, mass_2 * MSUN_SI)
}

/// `G / c^3`, the conversion from kilograms to seconds of the PN
/// coefficient functions.
#[inline]
fn kg_to_seconds() -> f64 {
    G_SI / C_SI.powf(3.0)
}

/// An upper bound on the time from the gravitational-wave frequency
/// `f_start` (Hz) to coalescence of a point-particle inspiral: the 2PN
/// chirp time with the negative terms dropped and the 1.5PN spin term of
/// the larger spin (`XLALSimInspiralChirpTimeBound`).
pub fn chirp_time_bound(f_start: f64, mass_1: f64, mass_2: f64, spin_1: f64, spin_2: f64) -> f64 {
    let (m1, m2) = kilograms(mass_1, mass_2);
    let m = m1 + m2;
    let mu = m1 * m2 / m;
    let eta = mu / m;
    let chi = if spin_1.abs() > spin_2.abs() {
        spin_1
    } else {
        spin_2
    }
    .abs();
    // XLALSimInspiralTaylorT2Timing_0PNCoeff, _2PNCoeff, _4PNCoeff
    let c0 = (-5.0 * (m * kg_to_seconds()) / (256.0 * eta)).abs();
    let c2 = fma(11.0 / 3.0, eta, 7.43 / 2.52);
    let c3 = (226.0 / 15.0) * chi;
    let c4 = fma(
        61.7 / 7.2 * eta,
        eta,
        fma(54.29 / 5.04, eta, 30.58673 / 5.08032),
    );
    let v = (PI * G_SI * m * f_start).libm_cbrt() / C_SI;
    let inner = fma(fma(c4, v, c3), v, c2);
    c0 * v.powf(-8.0) * fma(inner * v, v, 1.0)
}

/// An upper bound on the plunge and merger time: one orbit at `r = 9 M`
/// with speed `c / 3` (`XLALSimInspiralMergeTimeBound`).
pub fn merge_time_bound(mass_1: f64, mass_2: f64) -> f64 {
    let (m1, m2) = kilograms(mass_1, mass_2);
    let m = m1 + m2;
    let r = 9.0 * m * MRSUN_SI / MSUN_SI;
    let v = C_SI / 3.0;
    1.0 * (TWO_PI * r / v)
}

/// An upper bound on the ringdown time: 11 e-folds of the fundamental
/// (2,2) quasinormal mode of a black hole of mass `total_mass` (solar
/// masses, the sum `m1 + m2` in kilograms is used, see
/// [`ringdown_time_bound_kg`]) and spin `final_spin`, from the fits of
/// Berti, Cardoso & Will (2006) (`XLALSimInspiralRingdownTimeBound`).
pub fn ringdown_time_bound(total_mass: f64, final_spin: f64) -> f64 {
    ringdown_time_bound_kg(total_mass * MSUN_SI, final_spin)
}

/// [`ringdown_time_bound`] for the mass in kilograms, as LAL takes it.
pub fn ringdown_time_bound_kg(total_mass_kg: f64, final_spin: f64) -> f64 {
    const NEFOLDS: f64 = 11.0;
    const F1: f64 = 1.5251;
    const F2: f64 = -1.1568;
    const F3: f64 = 0.1292;
    const Q1: f64 = 0.7000;
    const Q2: f64 = 1.4187;
    const Q3: f64 = -0.4990;
    let omega = fma(F2, (1.0 - final_spin).powf(F3), F1) / (total_mass_kg * MTSUN_SI / MSUN_SI);
    let q = fma(Q2, (1.0 - final_spin).powf(Q3), Q1);
    let tau = 2.0 * q / omega;
    NEFOLDS * tau
}

/// An upper bound on the dimensionless spin of the remnant black hole
/// (Tichy & Marronetti 2008, or the larger spin), at most 0.998
/// (`XLALSimInspiralFinalBlackHoleSpinBound`).
pub fn final_black_hole_spin_bound(spin_1z: f64, spin_2z: f64) -> f64 {
    const MAXIMUM_BLACK_HOLE_SPIN: f64 = 0.998;
    let mut s = fma(0.15, spin_1z + spin_2z, 0.686);
    if s < spin_1z.abs() {
        s = spin_1z.abs();
    }
    if s < spin_2z.abs() {
        s = spin_2z.abs();
    }
    if s > MAXIMUM_BLACK_HOLE_SPIN {
        s = MAXIMUM_BLACK_HOLE_SPIN;
    }
    s
}

/// A lower bound on the frequency at which the Newtonian chirp lasts
/// `t_chirp` seconds (`XLALSimInspiralChirpStartFrequencyBound`).
pub fn chirp_start_frequency_bound(t_chirp: f64, mass_1: f64, mass_2: f64) -> f64 {
    let (m1, m2) = kilograms(mass_1, mass_2);
    let m = m1 + m2;
    let mu = m1 * m2 / m;
    let eta = mu / m;
    // XLALSimInspiralTaylorT3Frequency_0PNCoeff
    let c0 = 1.0 / (8.0 * PI * (m * kg_to_seconds()));
    c0 * (5.0 * m * (MTSUN_SI / MSUN_SI) / (eta * t_chirp)).powf(3.0 / 8.0)
}

/// The lowest Kerr ISCO frequency, `1 / (9^(3/2) pi M)`, below which the
/// conditioning routines start the waveform.
pub fn kerr_isco_frequency(mass_1: f64, mass_2: f64) -> f64 {
    let (m1, m2) = kilograms(mass_1, mass_2);
    1.0 / (9.0f64.powf(1.5) * PI * (m1 + m2) * MTSUN_SI / MSUN_SI)
}

/// The Schwarzschild ISCO frequency, `1 / (6^(3/2) pi M)`, at which the
/// end of the waveform is tapered.
pub fn schwarzschild_isco_frequency(mass_1: f64, mass_2: f64) -> f64 {
    let (m1, m2) = kilograms(mass_1, mass_2);
    1.0 / (6.0f64.powf(1.5) * PI * (m1 + m2) * MTSUN_SI / MSUN_SI)
}

/// The ends of a waveform tapered by [`wave_taper`]
/// (`LALSimInspiralApplyTaper`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum WaveTaper {
    /// `LAL_SIM_INSPIRAL_TAPER_START`
    Start,
    /// `LAL_SIM_INSPIRAL_TAPER_END`
    End,
    /// `LAL_SIM_INSPIRAL_TAPER_STARTEND`
    StartEnd,
}

/// Local extrema closer than this to the end of the signal are ignored
/// (`LALSIMULATION_RINGING_EXTENT`).
const RINGING_EXTENT: i64 = 19;

/// `XLALSimInspiralREAL8WaveTaper`: a Planck-like taper
/// `1 / (exp(z) + 1)` from the first (last) non-zero sample to the second
/// local extremum of `|h|` more than 19 samples in (or to the middle of
/// the signal), which is how `XLALSimInspiralTD` and `XLALSimInspiralFD`
/// condition models that start at the requested frequency (NRSur7dq4).
///
/// Does nothing if there is no or too short (fewer than 3 samples) a
/// signal.
pub fn wave_taper(data: &mut [f64], taper: WaveTaper) {
    let length = data.len() as i64;
    let Some(start) = data.iter().position(|&x| x != 0.0) else {
        return;
    };
    let start = start as i64;
    let end = data.iter().rposition(|&x| x != 0.0).unwrap() as i64;
    if end - start <= 1 {
        return;
    }
    fn at(data: &[f64], i: i64) -> f64 {
        data[i as usize].abs()
    }
    let mid = (start + end) / 2;

    if taper != WaveTaper::End {
        let mut flag = 0;
        let mut i = start + 1;
        let mut n = 0;
        while flag < 2 && i != mid && i + 1 < length {
            if at(data, i) >= at(data, i - 1) && at(data, i) >= at(data, i + 1) {
                if at(data, i) == at(data, i + 1) {
                    i += 1;
                }
                if i - start > RINGING_EXTENT {
                    flag += 1;
                }
                n = i - start;
            }
            i += 1;
        }
        if flag < 2 {
            n = mid - start;
        }
        let real_n = n as f64;
        data[start as usize] = 0.0;
        let mut i = start + 1;
        while i < start + n - 1 {
            let real_i = (i - start) as f64;
            let z = (real_n - 1.0) / real_i + (real_n - 1.0) / (real_i - (real_n - 1.0));
            let sigma = 1.0 / (z.exp() + 1.0);
            data[i as usize] *= sigma;
            i += 1;
        }
    }

    if taper != WaveTaper::Start {
        let mut i = end - 1;
        let mut flag = 0;
        let mut n = 0;
        while flag < 2 && i != mid && i >= 1 {
            if at(data, i) >= at(data, i + 1) && at(data, i) >= at(data, i - 1) {
                if at(data, i) == at(data, i - 1) {
                    i -= 1;
                }
                if end - i > RINGING_EXTENT {
                    flag += 1;
                }
                n = end - i;
            }
            i -= 1;
        }
        if flag < 2 {
            n = end - mid;
        }
        let real_n = n as f64;
        data[end as usize] = 0.0;
        let mut i = end - 1;
        while i > end - n + 1 {
            let real_i = (end - i) as f64;
            let z = (real_n - 1.0) / real_i + (real_n - 1.0) / (real_i - (real_n - 1.0));
            let sigma = 1.0 / (z.exp() + 1.0);
            data[i as usize] *= sigma;
            i -= 1;
        }
    }
}

/// The Hann window `0.5 - 0.5 cos(pi j / n)` of the conditioning stages.
#[inline]
fn hann(j: usize, n: usize) -> f64 {
    0.5 - 0.5 * (j as f64 * PI / n as f64).cos()
}

/// Multiply the first `n` samples by `0.5 - 0.5 cos(pi j / n)`, `j = 0 ..
/// n - 1` (the start taper of the conditioning stages).
pub fn hann_taper_start(data: &mut [f64], n: usize) {
    for (j, x) in data.iter_mut().take(n).enumerate() {
        *x *= hann(j, n);
    }
}

/// Multiply the samples `len - j`, `j = 1 .. n - 1`, by `0.5 - 0.5 cos(pi
/// j / n)`: the end taper of [`td_condition_stage_2`], which makes the
/// waveform vanish one sample after its end.
pub fn hann_taper_end(data: &mut [f64], n: usize) {
    let length = data.len();
    for j in 1..n.min(length + 1) {
        data[length - j] *= hann(j, n);
    }
}

/// Errors of the conditioning and conversions.
#[derive(Clone, Debug, PartialEq)]
pub enum ConditioningError {
    /// The high-pass frequency is not within the Nyquist band
    /// (`XLALParsePassBandParamStruc`, `XLAL_EINVAL`).
    FilterFrequency { frequency: f64, delta_t: f64 },
    /// Invalid arguments.
    InvalidArgument(String),
}

impl std::fmt::Display for ConditioningError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::FilterFrequency { frequency, delta_t } => write!(
                f,
                "filter frequency {frequency} Hz is not within the Nyquist band of a series \
                 sampled at {delta_t} s"
            ),
            Self::InvalidArgument(why) => write!(f, "{why}"),
        }
    }
}

impl std::error::Error for ConditioningError {}

/// A real IIR filter (`REAL8IIRFilter`): `direct` and `recursive`
/// coefficients as LAL stores them (`recursive[0] = -1`).
#[derive(Clone, Debug)]
struct IirFilter {
    direct: Vec<f64>,
    recursive: Vec<f64>,
}

/// `XLALWToZCOMPLEX16ZPGFilter` for one zero or pole `w`, giving its image
/// in the `z` plane and its gain correction factor (`None` for `w = -i`,
/// which maps to `z = infinity`).
fn bilinear(w: Complex<f64>, is_pole: bool) -> (Option<Complex<f64>>, Complex<f64>) {
    let (ar, ai) = (w.re, w.im);
    if ar == 0.0 {
        if ai == -1.0 {
            let g = if is_pole {
                Complex::new(0.0, -0.5)
            } else {
                Complex::new(0.0, 2.0)
            };
            (None, g)
        } else {
            let z = Complex::new((1.0 - ai) / (1.0 + ai), 0.0);
            let g = if is_pole {
                Complex::new(0.0, 1.0 / (1.0 + ai))
            } else {
                Complex::new(0.0, -(1.0 + ai))
            };
            (Some(z), g)
        }
    } else if (1.0 + ai).abs() > ar.abs() {
        let ratio = -ar / (1.0 + ai);
        let denom = fma(-ratio, ar, 1.0 + ai);
        let z = Complex::new(
            fma(ratio, ar, 1.0 - ai) / denom,
            fma(-ratio, 1.0 - ai, ar) / denom,
        );
        let g = if is_pole {
            Complex::new(ratio / denom, 1.0 / denom)
        } else {
            Complex::new(-ar, -(1.0 + ai))
        };
        (Some(z), g)
    } else {
        let ratio = -(1.0 + ai) / ar;
        let denom = fma(ratio, 1.0 + ai, -ar);
        let z = Complex::new(
            fma(1.0 - ai, ratio, ar) / denom,
            (fms(ar, ratio, 1.0) + ai) / denom,
        );
        let g = if is_pole {
            Complex::new(ratio / denom, 1.0 / denom)
        } else {
            Complex::new(-ar, -(1.0 + ai))
        };
        (Some(z), g)
    }
}

/// The IIR filter of a section with the given `w`-plane zeros, poles and
/// gain: `XLALWToZCOMPLEX16ZPGFilter` followed by
/// `XLALCreateREAL8IIRFilter`.
fn iir_from_wplane(
    zeros: &[Complex<f64>],
    poles: &[Complex<f64>],
    gain: Complex<f64>,
) -> IirFilter {
    let num = zeros.len().max(poles.len());
    let mut z_zeros = Vec::new();
    let mut z_poles = Vec::new();
    let mut factors = Vec::new();
    for &w in zeros {
        let (z, g) = bilinear(w, false);
        z_zeros.extend(z);
        factors.push(g);
    }
    let finite_zeros = num - zeros.iter().filter(|w| w.re == 0.0 && w.im == -1.0).count();
    z_zeros.resize(finite_zeros, Complex::new(-1.0, 0.0));
    for &w in poles {
        let (z, g) = bilinear(w, true);
        z_poles.extend(z);
        factors.push(g);
    }
    let finite_poles = num - poles.iter().filter(|w| w.re == 0.0 && w.im == -1.0).count();
    z_poles.resize(finite_poles, Complex::new(-1.0, 0.0));

    // multiply alternately by small and large factors, as LAL does
    let mut order: Vec<usize> = (0..factors.len()).collect();
    order.sort_by(|&a, &b| factors[a].norm_sqr().total_cmp(&factors[b].norm_sqr()));
    let mut gain = gain;
    let (mut i, mut j) = (0, factors.len().wrapping_sub(1));
    while !factors.is_empty() && i < j {
        gain *= factors[order[i]] * factors[order[j]];
        i += 1;
        j -= 1;
    }
    if !factors.is_empty() && i == j {
        gain *= factors[order[i]];
    }

    // XLALCreateREAL8IIRFilter
    let mut num_direct = 1usize;
    for z in &z_zeros {
        if z.im == 0.0 && z.re == 0.0 {
            num_direct -= 1;
        }
    }
    let mut num_recurs = 1usize;
    for p in &z_poles {
        if p.im == 0.0 && p.re == 0.0 {
            num_recurs -= 1;
        }
    }
    let num = z_poles.len().min(z_zeros.len());
    num_direct += num;
    num_recurs += num;
    let mut recursive = vec![0.0; num_recurs];
    recursive[0] = -1.0;
    for p in &z_poles {
        let (x, y) = (p.re, p.im);
        let mut j = num_recurs - 1;
        if y == 0.0 {
            while j > 0 {
                recursive[j] = fma(-x, recursive[j - 1], recursive[j]);
                j -= 1;
            }
        } else if y > 0.0 {
            while j > 1 {
                let modulus = fma(x, x, y * y);
                let t = fma(2.0 * x, recursive[j - 1], -(modulus * recursive[j - 2]));
                recursive[j] -= t;
                j -= 1;
            }
            recursive[j] = fma(-(2.0 * x), recursive[j - 1], recursive[j]);
        }
    }
    let mut direct = vec![0.0; num_direct];
    let offset = num - z_zeros.len().min(num);
    direct[offset] = gain.re;
    for zero in &z_zeros {
        let (x, y) = (zero.re, zero.im);
        let mut j = num_direct - 1;
        if y == 0.0 {
            while j > offset {
                direct[j] = fma(-x, direct[j - 1], direct[j]);
                j -= 1;
            }
        } else if y > 0.0 {
            while j > offset + 1 {
                let modulus = fma(x, x, y * y);
                let t = fma(2.0 * x, direct[j - 1], -(modulus * direct[j - 2]));
                direct[j] -= t;
                j -= 1;
            }
            direct[j] = fma(-(2.0 * x), direct[j - 1], direct[j]);
        }
    }
    IirFilter { direct, recursive }
}

impl IirFilter {
    /// `XLALIIRFilterREAL8Vector` (from zero history) when `reverse` is
    /// false, `XLALIIRFilterReverseREAL8Vector` when it is true.
    fn apply(&self, data: &mut [f64], reverse: bool) {
        if self.is_second_order() {
            return self.apply_second_order([data], reverse);
        }
        let num_hist = self.recursive.len() + 1;
        let recurs_order = self.recursive.len();
        let mut temp = vec![0.0; num_hist];
        let mut step = |x: &mut f64| {
            let mut j = num_hist - 1;
            while j >= recurs_order {
                temp[j] = temp[j - 1];
                j -= 1;
            }
            let mut w = *x;
            while j > 0 {
                temp[j] = temp[j - 1];
                w = fma(self.recursive[j], temp[j], w);
                j -= 1;
            }
            temp[0] = w;
            let mut datum = self.direct[0] * w;
            for (d, t) in self.direct.iter().zip(&temp).skip(1) {
                datum = fma(*d, *t, datum);
            }
            *x = datum;
        };
        if reverse {
            data.iter_mut().rev().for_each(&mut step);
        } else {
            data.iter_mut().for_each(&mut step);
        }
    }

    /// Whether this is a section with three direct and three recursive
    /// coefficients (the second-order sections of the Butterworth filter).
    fn is_second_order(&self) -> bool {
        self.direct.len() == 3 && self.recursive.len() == 3
    }

    /// [`Self::apply`] for a second-order section, on `N` sequences of the
    /// same length at once: the same operations as LAL for each sequence,
    /// interleaved so that the independent recursions overlap.
    fn apply_second_order<const N: usize>(&self, data: [&mut [f64]; N], reverse: bool) {
        let length = data[0].len();
        assert!(data.iter().all(|d| d.len() == length));
        let (d0, d1, d2) = (self.direct[0], self.direct[1], self.direct[2]);
        let (r1, r2) = (self.recursive[1], self.recursive[2]);
        // temp[0], temp[1], temp[2] of LAL's filter
        let mut t0 = [0.0; N];
        let mut t1 = [0.0; N];
        let mut t2 = [0.0; N];
        for step in 0..length {
            let index = if reverse { length - 1 - step } else { step };
            for k in 0..N {
                t2[k] = t1[k];
                let w = fma(r2, t2[k], data[k][index]);
                t1[k] = t0[k];
                let w = fma(r1, t1[k], w);
                t0[k] = w;
                data[k][index] = fma(d2, t2[k], fma(d1, t1[k], d0 * w));
            }
        }
    }
}

/// The sections of `XLALButterworthREAL8TimeSeries` for a high-pass
/// filter (`XLALParsePassBandParamStruc` with `f1 = a1 = -1`).
fn high_pass_sections(
    delta_t: f64,
    frequency: f64,
    attenuation: f64,
    order: u32,
) -> Result<Vec<IirFilter>, ConditioningError> {
    let w2 = delta_t * frequency;
    let w2_given = w2 > 0.0 && w2 < 0.5;
    if !w2_given || order == 0 {
        return Err(ConditioningError::FilterFrequency { frequency, delta_t });
    }
    let n = order as usize;
    let wc = if attenuation > 0.0 && attenuation < 1.0 {
        let w_high = (PI * w2).tan();
        let a_high = attenuation.sqrt();
        w_high * (1.0 / a_high - 1.0).powf(0.5 / n as f64)
    } else {
        (PI * w2).tan()
    };
    let mut sections = Vec::with_capacity(n / 2 + 1);
    let (mut i, mut j) = (0usize, n - 1);
    while i < j {
        let theta = PI * (i as f64 + 0.5) / n as f64;
        let ar = wc * theta.cos();
        let ai = wc * theta.sin();
        sections.push(iir_from_wplane(
            &[Complex::new(0.0, 0.0); 2],
            &[Complex::new(ar, ai), Complex::new(-ar, ai)],
            Complex::new(1.0, 0.0),
        ));
        i += 1;
        j -= 1;
    }
    if i == j {
        sections.push(iir_from_wplane(
            &[Complex::new(0.0, 0.0)],
            &[Complex::new(0.0, wc)],
            Complex::new(1.0, 0.0),
        ));
    }
    Ok(sections)
}

/// `XLALHighPassREAL8TimeSeries(series, frequency, attenuation, order)`:
/// an order-`order` Butterworth high-pass filter with attenuation
/// `attenuation` at `frequency` (Hz), applied forwards and backwards
/// (zero phase) section by section.
///
/// The conditioning routines use `attenuation = 0.99` and `order = 8`. If
/// `attenuation` is not in (0, 1), `frequency` is the characteristic
/// frequency of the filter, as in LAL. The result is bit for bit that of
/// LAL on Apple silicon.
pub fn high_pass(
    series: &mut TimeSeries,
    frequency: f64,
    attenuation: f64,
    order: u32,
) -> Result<(), ConditioningError> {
    for section in high_pass_sections(series.delta_t, frequency, attenuation, order)? {
        section.apply(&mut series.data, false);
        section.apply(&mut series.data, true);
    }
    Ok(())
}

/// [`high_pass`] of two series (the polarizations), with the same result
/// as filtering them one after the other, but faster when they have the
/// same sampling and length.
pub fn high_pass_pair(
    hplus: &mut TimeSeries,
    hcross: &mut TimeSeries,
    frequency: f64,
    attenuation: f64,
    order: u32,
) -> Result<(), ConditioningError> {
    if hplus.delta_t != hcross.delta_t || hplus.data.len() != hcross.data.len() {
        high_pass(hplus, frequency, attenuation, order)?;
        return high_pass(hcross, frequency, attenuation, order);
    }
    for section in high_pass_sections(hplus.delta_t, frequency, attenuation, order)? {
        for reverse in [false, true] {
            if section.is_second_order() {
                section.apply_second_order([&mut hplus.data, &mut hcross.data], reverse);
            } else {
                section.apply(&mut hplus.data, reverse);
                section.apply(&mut hcross.data, reverse);
            }
        }
    }
    Ok(())
}

/// `XLALSimInspiralTDConditionStage1`: taper the first `t_extra` seconds
/// with a Hann window, high-pass filter at `f_min` and remove the zeros
/// padding the end of the waveform.
pub fn td_condition_stage_1(
    hplus: &mut TimeSeries,
    hcross: &mut TimeSeries,
    t_extra: f64,
    f_min: f64,
) {
    let length = hplus.data.len();
    let mut nzeros = 0;
    while nzeros < length
        && hplus.data[length - nzeros - 1] == 0.0
        && hcross.data[hcross.data.len() - nzeros - 1] == 0.0
    {
        nzeros += 1;
    }
    let ntaper = (t_extra / hplus.delta_t).round() as usize;
    hann_taper_start(&mut hplus.data, ntaper);
    hann_taper_start(&mut hcross.data, ntaper);
    // LAL ignores failures of the filter here
    let _ = high_pass_pair(hplus, hcross, f_min, 0.99, 8);
    if nzeros > 0 {
        hplus.data.truncate(length - nzeros);
        let n = hcross.data.len();
        hcross.data.truncate(n - nzeros);
    }
}

/// `XLALSimInspiralTDConditionStage2`: taper one cycle at `f_max` at the
/// end and one cycle at `f_min` at the start (at least 4 samples each);
/// nothing for waveforms shorter than 8 samples.
pub fn td_condition_stage_2(
    hplus: &mut TimeSeries,
    hcross: &mut TimeSeries,
    f_min: f64,
    f_max: f64,
) {
    const MIN_TAPER_SAMPLES: usize = 4;
    if hplus.data.len() < 2 * MIN_TAPER_SAMPLES {
        return;
    }
    let ntaper = ((1.0 / (f_max * hplus.delta_t)).round() as usize).max(MIN_TAPER_SAMPLES);
    hann_taper_end(&mut hplus.data, ntaper);
    hann_taper_end(&mut hcross.data, ntaper);
    let ntaper = ((1.0 / (f_min * hplus.delta_t)).round() as usize).max(MIN_TAPER_SAMPLES);
    hann_taper_start(&mut hplus.data, ntaper);
    hann_taper_start(&mut hcross.data, ntaper);
}

/// The frequency-domain conditioning of `XLALSimInspiralFD`: zero below
/// `k0 = round(f_start / delta_f)`, multiply `k0 <= k < k1 =
/// round(f_min / delta_f)` by `0.5 - 0.5 cos(pi (k - k0) / (k1 - k0))` and
/// zero the last (Nyquist) sample.
pub fn taper_frequency_band(series: &mut FrequencySeries, f_start: f64, f_min: f64) {
    let k0 = (f_start / series.delta_f).round() as usize;
    let k1 = (f_min / series.delta_f).round() as usize;
    let length = series.data.len();
    let mut k = 0;
    while k < k0.min(length) {
        series.data[k] = Complex::new(0.0, 0.0);
        k += 1;
    }
    while k < k1.min(length) {
        let w = 0.5 - 0.5 * (PI * (k - k0) as f64 / (k1 - k0) as f64).cos();
        series.data[k] *= w;
        k += 1;
    }
    if let Some(last) = series.data.last_mut() {
        *last = Complex::new(0.0, 0.0);
    }
}

/// Shift a frequency series by `t_shift` seconds backwards in time and
/// compensate in the epoch, as the conditioning routines do: multiply
/// sample `k` by `cexp(2 pi i k delta_f t_shift)` (with `delta_f` given
/// explicitly since LAL sometimes uses a value other than the series') and
/// add `t_shift` to the epoch.
pub fn time_shift(series: &mut FrequencySeries, delta_f: f64, t_shift: f64) {
    for (k, h) in series.data.iter_mut().enumerate() {
        *h *= phase_factor(TWO_PI * k as f64 * delta_f * t_shift);
    }
    series.epoch = series.epoch.add_seconds(t_shift);
}

/// [`time_shift`] of both polarizations, computing each phase factor once
/// as LAL does.
pub fn time_shift_pair(
    hplus: &mut FrequencySeries,
    hcross: &mut FrequencySeries,
    delta_f: f64,
    t_shift: f64,
) {
    let common = hplus.data.len().min(hcross.data.len());
    for (k, (p, c)) in hplus
        .data
        .iter_mut()
        .zip(hcross.data.iter_mut())
        .enumerate()
    {
        let factor = phase_factor(TWO_PI * k as f64 * delta_f * t_shift);
        *p *= factor;
        *c *= factor;
    }
    for (series, rest) in [(&mut *hplus, common), (&mut *hcross, common)] {
        for (k, h) in series.data.iter_mut().enumerate().skip(rest) {
            *h *= phase_factor(TWO_PI * k as f64 * delta_f * t_shift);
        }
    }
    hplus.epoch = hplus.epoch.add_seconds(t_shift);
    hcross.epoch = hcross.epoch.add_seconds(t_shift);
}

/// `cexp(I phase)` as the C library computes it.
#[inline]
fn phase_factor(phase: f64) -> Complex<f64> {
    crate::imrphenomx::libm::cexp_i(phase)
}

/// The index of the largest `h_+^2 + h_x^2`, the last one among equal
/// values, as `XLALSimInspiralChooseTDWaveform` finds it to put the peak of
/// IMRPhenomD, IMRPhenomXAS and IMRPhenomXHM at `t = 0` (the length minus
/// one if all samples vanish).
pub fn peak_index(hplus: &[f64], hcross: &[f64]) -> usize {
    let mut max_amp = 0.0;
    let mut max_index = hplus.len().saturating_sub(1);
    for (i, (p, c)) in hplus.iter().zip(hcross).enumerate().rev() {
        let amp = fma(*p, *p, c * c);
        if amp > max_amp {
            max_index = i;
            max_amp = amp;
        }
    }
    max_index
}

/// `XLALGPSAdd(&zero, -1 / delta_f)`: the epoch of the frequency-domain
/// models of LALSimulation on a uniform grid (coalescence at `t = 0`).
pub fn fd_model_epoch(delta_f: f64) -> GpsTime {
    GpsTime::ZERO.add_seconds(-1.0 / delta_f)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bounds_are_sensible() {
        // GW150914-like: the Newtonian chirp time from 20 Hz is 0.84 s
        let t = chirp_time_bound(20.0, 36.0, 29.0, 0.3, -0.4);
        assert!(t > 0.84 && t < 2.0, "{t}");
        let s = final_black_hole_spin_bound(0.3, -0.4);
        assert_eq!(s, 0.686 + 0.15 * (0.3 - 0.4));
        assert_eq!(final_black_hole_spin_bound(0.99, 0.99), 0.99);
        assert_eq!(final_black_hole_spin_bound(1.0, 0.9), 0.998);
        let tm = merge_time_bound(36.0, 29.0) + ringdown_time_bound(65.0, s);
        assert!(tm > 0.0 && tm < 0.1, "{tm}");
        let f = chirp_start_frequency_bound(1.1 * t, 36.0, 29.0);
        assert!(f < 20.0 && f > 15.0, "{f}");
        assert!(kerr_isco_frequency(36.0, 29.0) < schwarzschild_isco_frequency(36.0, 29.0));
    }

    #[test]
    fn high_pass_removes_low_frequencies() {
        let delta_t = 1.0 / 1024.0;
        let n = 8192;
        let low = |t: f64| (2.0 * PI * 2.0 * t).sin();
        let high = |t: f64| (2.0 * PI * 100.0 * t).sin();
        let mut series = TimeSeries {
            epoch: GpsTime::ZERO,
            delta_t,
            data: (0..n)
                .map(|j| low(j as f64 * delta_t) + high(j as f64 * delta_t))
                .collect(),
        };
        high_pass(&mut series, 20.0, 0.99, 8).unwrap();
        // away from the edges, only the 100 Hz sinusoid is left
        for j in n / 4..3 * n / 4 {
            let expected = high(j as f64 * delta_t);
            assert!((series.data[j] - expected).abs() < 1e-2, "{j}");
        }
    }

    #[test]
    fn wave_taper_starts_at_zero() {
        let mut data: Vec<f64> = (0..400).map(|j| (0.3 * j as f64).sin() + 2.0).collect();
        wave_taper(&mut data, WaveTaper::StartEnd);
        assert_eq!(data[0], 0.0);
        assert_eq!(data[399], 0.0);
        assert!(data[1] < 1e-6);
        assert_eq!(data[200], (0.3 * 200.0f64).sin() + 2.0);
    }
}

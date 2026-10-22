// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2024-2026 Colm Talbot. Parts of waveformrs are translated
// from LALSuite, GSL and fdlibm: see NOTICE.md for their authors.
//! Real Fourier transforms with LAL's conventions (feature `fft`).
//!
//! LAL transforms real sequences with FFTW's half-complex plans
//! (`XLALREAL8ForwardFFT`, `XLALREAL8ReverseFFT`) and wraps them for time
//! and frequency series (`XLALREAL8TimeFreqFFT`, `XLALREAL8FreqTimeFFT` in
//! `TimeFreqFFT.c`):
//!
//! - forward: `X[k] = sum_j x[j] exp(-2 pi i j k / N)` for `k = 0 ..
//!   N / 2`, unnormalised; [`time_to_frequency`] multiplies by `delta_t`,
//!   keeps the epoch and sets `delta_f = 1 / (delta_t N)`;
//! - reverse: `x[j] = sum_{k=0}^{N-1} X[k] exp(2 pi i j k / N)` with the
//!   negative frequencies the complex conjugates of the positive ones and
//!   the imaginary parts of the zero and (for even `N`) Nyquist samples
//!   ignored, unnormalised; [`frequency_to_time`] multiplies by `delta_f`,
//!   keeps the epoch and sets `delta_t = 1 / (delta_f N)`.
//!
//! The transforms use [`rustfft`] with the usual packing of a real sequence
//! of even length `N` into a complex one of length `N / 2` (odd lengths use
//! a complex transform of length `N`), so any length works. They are not
//! bit for bit identical to FFTW: expect differences of a few `1e-16`
//! relative to the largest value of the output.

use std::collections::HashMap;
use std::f64::consts::PI;
use std::sync::{Arc, Mutex};

use num_complex::Complex;
use rustfft::{Fft, FftPlanner};

use crate::series::{FrequencySeries, TimeSeries};

type C64 = Complex<f64>;

/// A forward and reverse real FFT of a given length, with its twiddle
/// factors (`REAL8FFTPlan`). Plans are cheap to share between threads.
pub struct RealFft {
    length: usize,
    kind: Kind,
}

enum Kind {
    /// Even length `N`: a complex FFT of length `N / 2` and the twiddles
    /// `exp(-2 pi i k / N)`, `k = 0 .. N / 2`.
    Even {
        forward: Arc<dyn Fft<f64>>,
        inverse: Arc<dyn Fft<f64>>,
        twiddles: Vec<C64>,
    },
    /// Odd length: a complex FFT of the full length.
    Odd {
        forward: Arc<dyn Fft<f64>>,
        inverse: Arc<dyn Fft<f64>>,
    },
}

/// `exp(-2 pi i k / n)`, using the symmetries of the circle so that the
/// argument of `sin` and `cos` stays within `pi / 4`.
fn twiddle(k: usize, n: usize) -> C64 {
    // reduce k / n to [0, 1/8] and reconstruct with exact symmetries
    let k = k % n;
    let (eighths, rest) = ((8 * k) / n, (8 * k) % n);
    let angle = |num: usize| 2.0 * PI * num as f64 / (8 * n) as f64;
    // angle within the octant, measured from the octant boundary
    let (s, c) = if eighths % 2 == 0 {
        angle(rest).sin_cos()
    } else {
        angle(n - rest).sin_cos()
    };
    // (cos, sin) of 2 pi k / n
    let (cos, sin) = match eighths {
        0 => (c, s),
        1 => (s, c),
        2 => (-s, c),
        3 => (-c, s),
        4 => (-c, -s),
        5 => (-s, -c),
        6 => (s, -c),
        _ => (c, -s),
    };
    C64::new(cos, -sin)
}

impl RealFft {
    /// Plan the transforms of real sequences of length `length`.
    ///
    /// # Panics
    ///
    /// If `length` is zero.
    pub fn new(length: usize) -> Self {
        assert!(length > 0, "FFT length must be positive");
        let mut planner = FftPlanner::new();
        let kind = if length % 2 == 0 {
            let half = length / 2;
            Kind::Even {
                forward: planner.plan_fft_forward(half),
                inverse: planner.plan_fft_inverse(half),
                twiddles: (0..=half).map(|k| twiddle(k, length)).collect(),
            }
        } else {
            Kind::Odd {
                forward: planner.plan_fft_forward(length),
                inverse: planner.plan_fft_inverse(length),
            }
        };
        Self { length, kind }
    }

    /// The length of the real sequences.
    pub fn len(&self) -> usize {
        self.length
    }

    /// Whether the length is zero (never: plans have positive length).
    pub fn is_empty(&self) -> bool {
        false
    }

    /// The number of complex frequency samples, `N / 2 + 1`.
    pub fn complex_len(&self) -> usize {
        self.length / 2 + 1
    }

    /// `XLALREAL8ForwardFFT`: the unnormalised DFT `X[k] = sum_j x[j]
    /// exp(-2 pi i j k / N)` of `input` (length `N`) for `k = 0 .. N / 2`.
    ///
    /// # Panics
    ///
    /// If the lengths do not match the plan.
    pub fn forward(&self, input: &[f64], output: &mut [C64]) {
        assert_eq!(input.len(), self.length, "input length must match the plan");
        assert_eq!(
            output.len(),
            self.complex_len(),
            "output length must be N / 2 + 1"
        );
        match &self.kind {
            Kind::Even {
                forward, twiddles, ..
            } => {
                let half = self.length / 2;
                let mut z: Vec<C64> = input
                    .chunks_exact(2)
                    .map(|pair| C64::new(pair[0], pair[1]))
                    .collect();
                forward.process(&mut z);
                output[0] = C64::new(z[0].re + z[0].im, 0.0);
                output[half] = C64::new(z[0].re - z[0].im, 0.0);
                for k in 1..half {
                    let a = z[k];
                    let b = z[half - k].conj();
                    let even = (a + b) * 0.5;
                    let odd = (a - b) * C64::new(0.0, -0.5);
                    output[k] = even + twiddles[k] * odd;
                }
            }
            Kind::Odd { forward, .. } => {
                let mut z: Vec<C64> = input.iter().map(|&x| C64::new(x, 0.0)).collect();
                forward.process(&mut z);
                output.copy_from_slice(&z[..self.complex_len()]);
                output[0].im = 0.0;
            }
        }
    }

    /// `XLALREAL8ReverseFFT`: the unnormalised inverse DFT of the
    /// Hermitian sequence whose non-negative frequencies are `input`
    /// (length `N / 2 + 1`), `x[j] = sum_{k=0}^{N-1} X[k] exp(2 pi i j k /
    /// N)`. The imaginary parts of `input[0]` and (for even `N`)
    /// `input[N / 2]` are ignored.
    ///
    /// # Panics
    ///
    /// If the lengths do not match the plan.
    pub fn inverse(&self, input: &[C64], output: &mut [f64]) {
        assert_eq!(
            input.len(),
            self.complex_len(),
            "input length must be N / 2 + 1"
        );
        assert_eq!(
            output.len(),
            self.length,
            "output length must match the plan"
        );
        match &self.kind {
            Kind::Even {
                inverse, twiddles, ..
            } => {
                let half = self.length / 2;
                let mut z = vec![C64::new(0.0, 0.0); half];
                let (x0, xn) = (input[0].re, input[half].re);
                z[0] = C64::new(x0 + xn, x0 - xn);
                for k in 1..half {
                    let a = input[k];
                    let b = input[half - k].conj();
                    let even = a + b;
                    let odd = (a - b) * twiddles[k].conj();
                    z[k] = even + C64::new(-odd.im, odd.re);
                }
                inverse.process(&mut z);
                for (pair, w) in output.chunks_exact_mut(2).zip(&z) {
                    pair[0] = w.re;
                    pair[1] = w.im;
                }
            }
            Kind::Odd { inverse, .. } => {
                let n = self.length;
                let mut z = vec![C64::new(0.0, 0.0); n];
                z[0] = C64::new(input[0].re, 0.0);
                for k in 1..self.complex_len() {
                    z[k] = input[k];
                    z[n - k] = input[k].conj();
                }
                inverse.process(&mut z);
                for (x, w) in output.iter_mut().zip(&z) {
                    *x = w.re;
                }
            }
        }
    }
}

/// A cache of [`RealFft`] plans by length, shared by the conversions.
#[derive(Default)]
pub struct RealFftPlanner {
    plans: Mutex<HashMap<usize, Arc<RealFft>>>,
}

impl RealFftPlanner {
    /// An empty cache.
    pub fn new() -> Self {
        Self::default()
    }

    /// The plan for real sequences of length `length`, created on first
    /// use.
    pub fn plan(&self, length: usize) -> Arc<RealFft> {
        let mut plans = self.plans.lock().unwrap_or_else(|e| e.into_inner());
        plans
            .entry(length)
            .or_insert_with(|| Arc::new(RealFft::new(length)))
            .clone()
    }
}

/// The process-wide plan cache used by [`time_to_frequency`] and
/// [`frequency_to_time`].
pub fn global_planner() -> &'static RealFftPlanner {
    static PLANNER: once_cell::sync::Lazy<RealFftPlanner> =
        once_cell::sync::Lazy::new(RealFftPlanner::new);
    &PLANNER
}

/// `XLALREAL8TimeFreqFFT`: the Fourier transform `delta_t X[k]` of a time
/// series of length `N`, `k = 0 .. N / 2`, with the epoch of the time
/// series and `delta_f = 1 / (delta_t N)`.
///
/// # Panics
///
/// If the series is empty.
pub fn time_to_frequency(series: &TimeSeries) -> FrequencySeries {
    let plan = global_planner().plan(series.data.len());
    let mut data = vec![C64::new(0.0, 0.0); plan.complex_len()];
    plan.forward(&series.data, &mut data);
    for x in &mut data {
        *x *= series.delta_t;
    }
    FrequencySeries {
        epoch: series.epoch,
        delta_f: 1.0 / (series.delta_t * series.data.len() as f64),
        data,
    }
}

/// `XLALREAL8FreqTimeFFT`: the inverse Fourier transform `delta_f x[j]`
/// of a frequency series into a time series of length `length` (which must
/// satisfy `length / 2 + 1 == series.len()`), with the epoch of the
/// frequency series and `delta_t = 1 / (delta_f length)`.
///
/// # Panics
///
/// If the lengths are inconsistent.
pub fn frequency_to_time(series: &FrequencySeries, length: usize) -> TimeSeries {
    let plan = global_planner().plan(length);
    let mut data = vec![0.0; length];
    plan.inverse(&series.data, &mut data);
    for x in &mut data {
        *x *= series.delta_f;
    }
    TimeSeries {
        epoch: series.epoch,
        delta_t: 1.0 / (series.delta_f * length as f64),
        data,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn naive_dft(x: &[f64]) -> Vec<C64> {
        let n = x.len();
        (0..=n / 2)
            .map(|k| {
                x.iter()
                    .enumerate()
                    .map(|(j, &v)| {
                        let angle = -2.0 * PI * ((j * k) % n) as f64 / n as f64;
                        C64::new(v * angle.cos(), v * angle.sin())
                    })
                    .sum()
            })
            .collect()
    }

    #[test]
    fn twiddles_are_accurate() {
        for n in [1usize, 2, 3, 8, 12, 1000, 4096, 6000] {
            for k in 0..n {
                let angle = -2.0 * PI * k as f64 / n as f64;
                let expected = C64::new(angle.cos(), angle.sin());
                // the reference has the rounding error of its argument
                assert!((twiddle(k, n) - expected).norm() < 2e-15, "{k} {n}");
            }
        }
    }

    #[test]
    fn matches_naive_dft_and_round_trips() {
        for n in [1usize, 2, 3, 5, 6, 8, 30, 64, 98, 101, 1000] {
            let x: Vec<f64> = (0..n)
                .map(|j| ((j * j) as f64 * 0.37).sin() + 0.1)
                .collect();
            let plan = RealFft::new(n);
            let mut spectrum = vec![C64::new(0.0, 0.0); plan.complex_len()];
            plan.forward(&x, &mut spectrum);
            let expected = naive_dft(&x);
            for (a, b) in spectrum.iter().zip(&expected) {
                assert!((a - b).norm() < 1e-12 * n as f64, "n = {n}");
            }
            let mut back = vec![0.0; n];
            plan.inverse(&spectrum, &mut back);
            for (a, b) in back.iter().zip(&x) {
                assert!((a / n as f64 - b).abs() < 1e-14 * n as f64, "n = {n}");
            }
        }
    }

    #[test]
    fn series_conventions() {
        let series = TimeSeries {
            epoch: crate::series::GpsTime::from_seconds(-2.5),
            delta_t: 0.25,
            data: vec![1.0, 2.0, 0.0, -1.0, 0.5, 0.0, 0.0, 3.0],
        };
        let fd = time_to_frequency(&series);
        assert_eq!(fd.delta_f, 0.5);
        assert_eq!(fd.epoch, series.epoch);
        assert_eq!(fd.data.len(), 5);
        assert!((fd.data[0].re - 0.25 * 5.5).abs() < 1e-15);
        let td = frequency_to_time(&fd, 8);
        assert_eq!(td.delta_t, 0.25);
        for (a, b) in td.data.iter().zip(&series.data) {
            assert!((a - b).abs() < 1e-15);
        }
    }
}

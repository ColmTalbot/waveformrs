//! Uniformly sampled time and frequency series with LAL's epoch arithmetic.
//!
//! These mirror `REAL8TimeSeries` and `COMPLEX16FrequencySeries` (with
//! `f0 = 0`) closely enough that the conditioning and Fourier transforms of
//! [`crate::conditioning`], [`crate::fft`] and [`crate::conversions`] can
//! follow LALSimulation operation by operation. In particular the epoch is
//! a [`GpsTime`], an integer number of nanoseconds like LAL's
//! `LIGOTimeGPS`, so that shifting it (`XLALGPSAdd`) rounds exactly as LAL
//! does and epochs can be compared exactly.

use num_complex::Complex;

const BILLION: i64 = 1_000_000_000;

/// A GPS time with nanosecond resolution, LAL's `LIGOTimeGPS`.
///
/// Stored as a signed count of nanoseconds (`XLALGPSToINT8NS`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct GpsTime {
    nanoseconds: i64,
}

impl GpsTime {
    /// GPS time zero (`LIGOTIMEGPSZERO`).
    pub const ZERO: Self = Self { nanoseconds: 0 };

    /// The time `nanoseconds` after GPS time zero (`XLALINT8NSToGPS`).
    pub const fn from_nanoseconds(nanoseconds: i64) -> Self {
        Self { nanoseconds }
    }

    /// The time from LAL's `gpsSeconds` and `gpsNanoSeconds` fields
    /// (`XLALGPSSet`).
    pub const fn from_parts(seconds: i32, nanoseconds: i64) -> Self {
        Self {
            nanoseconds: BILLION * seconds as i64 + nanoseconds,
        }
    }

    /// `t` seconds rounded to the nearest nanosecond as `XLALGPSSetREAL8`
    /// does: `floor(t)` seconds plus `nearbyint((t - floor(t)) 1e9)`
    /// nanoseconds.
    ///
    /// # Panics
    ///
    /// If `t` is NaN or outside the range of a 32-bit count of seconds,
    /// where LAL fails.
    pub fn from_seconds(t: f64) -> Self {
        assert!(
            !t.is_nan() && t.abs() <= i32::MAX as f64,
            "GPS time {t} out of range"
        );
        let seconds = t.floor() as i32;
        let nanoseconds = ((t - seconds as f64) * 1e9).round_ties_even() as i32;
        Self::from_parts(seconds, nanoseconds as i64)
    }

    /// The number of nanoseconds since GPS time zero (`XLALGPSToINT8NS`).
    pub const fn nanoseconds(&self) -> i64 {
        self.nanoseconds
    }

    /// LAL's `gpsSeconds`: the seconds, truncated towards zero.
    pub const fn gps_seconds(&self) -> i64 {
        self.nanoseconds / BILLION
    }

    /// LAL's `gpsNanoSeconds`: the remaining nanoseconds, with the sign of
    /// the time.
    pub const fn gps_nanoseconds(&self) -> i64 {
        self.nanoseconds % BILLION
    }

    /// The time in seconds as `XLALGPSGetREAL8` (and `float()` of a
    /// `LIGOTimeGPS` in Python) computes it: `gpsSeconds + gpsNanoSeconds /
    /// 1e9`.
    pub fn seconds(&self) -> f64 {
        self.gps_seconds() as f64 + self.gps_nanoseconds() as f64 / 1e9
    }

    /// `XLALGPSAdd`: `dt` seconds, rounded to nanoseconds with
    /// [`GpsTime::from_seconds`], added to this time.
    pub fn add_seconds(self, dt: f64) -> Self {
        Self {
            nanoseconds: self.nanoseconds + Self::from_seconds(dt).nanoseconds,
        }
    }
}

/// A real time series (`REAL8TimeSeries`): `data[j]` is the value at time
/// `epoch + j delta_t`.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct TimeSeries {
    /// Time of the first sample.
    pub epoch: GpsTime,
    /// Sampling interval in seconds.
    pub delta_t: f64,
    /// The samples.
    pub data: Vec<f64>,
}

impl TimeSeries {
    /// A series of `length` zeros.
    pub fn zeros(epoch: GpsTime, delta_t: f64, length: usize) -> Self {
        Self {
            epoch,
            delta_t,
            data: vec![0.0; length],
        }
    }

    /// The number of samples.
    pub fn len(&self) -> usize {
        self.data.len()
    }

    /// Whether the series has no samples.
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    /// The sample times relative to GPS time zero, `epoch + j delta_t`.
    pub fn times(&self) -> impl Iterator<Item = f64> + '_ {
        let t0 = self.epoch.seconds();
        (0..self.data.len()).map(move |j| t0 + j as f64 * self.delta_t)
    }

    /// `XLALResizeREAL8TimeSeries(series, first, length)`: keep `length`
    /// samples starting at index `first` (which may be negative), padding
    /// with zeros where the old series has no samples, and move the epoch
    /// by `first delta_t` (rounded to nanoseconds).
    pub fn resize(&mut self, first: i64, length: usize) {
        resize_sequence(&mut self.data, first, length, 0.0);
        self.epoch = self.epoch.add_seconds(first as f64 * self.delta_t);
    }
}

/// A complex frequency series (`COMPLEX16FrequencySeries` with `f0 = 0`):
/// `data[k]` is the value at frequency `k delta_f`.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct FrequencySeries {
    /// The epoch of the time series this is the Fourier transform of.
    pub epoch: GpsTime,
    /// Frequency spacing in Hz.
    pub delta_f: f64,
    /// The samples, starting at zero frequency.
    pub data: Vec<Complex<f64>>,
}

impl FrequencySeries {
    /// A series of `length` zeros.
    pub fn zeros(epoch: GpsTime, delta_f: f64, length: usize) -> Self {
        Self {
            epoch,
            delta_f,
            data: vec![Complex::new(0.0, 0.0); length],
        }
    }

    /// The number of samples.
    pub fn len(&self) -> usize {
        self.data.len()
    }

    /// Whether the series has no samples.
    pub fn is_empty(&self) -> bool {
        self.data.is_empty()
    }

    /// The frequencies `k delta_f`.
    pub fn frequencies(&self) -> impl Iterator<Item = f64> + '_ {
        (0..self.data.len()).map(move |k| k as f64 * self.delta_f)
    }
}

/// Plus and cross polarizations in the time domain.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct TimeDomainPolarizations {
    pub hplus: TimeSeries,
    pub hcross: TimeSeries,
}

/// Plus and cross polarizations in the frequency domain.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct FrequencyDomainPolarizations {
    pub hplus: FrequencySeries,
    pub hcross: FrequencySeries,
}

/// `XLALResize*Sequence`: `new[j] = old[first + j]` where that exists,
/// `fill` elsewhere.
pub(crate) fn resize_sequence<T: Copy>(data: &mut Vec<T>, first: i64, length: usize, fill: T) {
    let old_length = data.len() as i64;
    if first >= 0 {
        let start = (first.min(old_length)) as usize;
        data.drain(..start);
        data.resize(length, fill);
    } else {
        let pad = (-first) as usize;
        let keep = old_length.min((length as i64 + first).max(0)) as usize;
        data.truncate(keep);
        let mut new = Vec::with_capacity(length);
        new.resize(pad.min(length), fill);
        new.extend_from_slice(data);
        new.resize(length, fill);
        *data = new;
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gps_rounding_matches_lal() {
        // XLALGPSSetREAL8: floor seconds, nearbyint nanoseconds
        let t = GpsTime::from_seconds(-0.9586624661);
        assert_eq!(t.nanoseconds(), -958_662_466);
        assert_eq!(t.gps_seconds(), 0);
        assert_eq!(t.gps_nanoseconds(), -958_662_466);
        assert_eq!(t.seconds(), -0.958662466);
        assert_eq!(GpsTime::from_seconds(-2.0000000004).seconds(), -2.0);
        // ties to even
        assert_eq!(GpsTime::from_seconds(2.5e-9).nanoseconds(), 2);
        let t = GpsTime::from_parts(-3, -218_017_578);
        assert_eq!(t.seconds(), -3.218017578);
        assert_eq!(
            GpsTime::ZERO.add_seconds(-8.0).add_seconds(0.25),
            GpsTime::from_parts(-7, -750_000_000)
        );
    }

    #[test]
    fn resize_pads_and_truncates() {
        let mut data = vec![1, 2, 3, 4];
        resize_sequence(&mut data, 1, 2, 0);
        assert_eq!(data, [2, 3]);
        let mut data = vec![1, 2, 3, 4];
        resize_sequence(&mut data, -2, 5, 0);
        assert_eq!(data, [0, 0, 1, 2, 3]);
        let mut data = vec![1, 2, 3, 4];
        resize_sequence(&mut data, -2, 8, 0);
        assert_eq!(data, [0, 0, 1, 2, 3, 4, 0, 0]);
        let mut data = vec![1, 2, 3, 4];
        resize_sequence(&mut data, 3, 3, 0);
        assert_eq!(data, [4, 0, 0]);
        let mut data = vec![1, 2, 3, 4];
        resize_sequence(&mut data, -6, 3, 0);
        assert_eq!(data, [0, 0, 0]);
        let mut data = vec![1, 2, 3, 4];
        resize_sequence(&mut data, 6, 2, 0);
        assert_eq!(data, [0, 0]);
        let mut series = TimeSeries {
            epoch: GpsTime::ZERO,
            delta_t: 0.25,
            data: vec![1.0; 4],
        };
        series.resize(-4, 8);
        assert_eq!(series.epoch, GpsTime::from_seconds(-1.0));
        assert_eq!(series.data, [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0]);
    }
}

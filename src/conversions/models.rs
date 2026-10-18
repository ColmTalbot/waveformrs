//! The models of this crate as [`FdGenerator`]s and [`TdGenerator`]s:
//! `SimInspiralChooseFDWaveform` on a uniform grid for the frequency-domain
//! approximants (built on [`crate::choose_fd_model`] and the uniform-grid
//! methods of IMRPhenomXP/XPHM), `SimInspiralChooseTDWaveform` for them
//! (`XLALSimInspiralTDFromFD` with LAL's peak alignment) and for the
//! time-domain models ([`crate::TimeDomainModel`]), and `SimInspiralTD` /
//! `SimInspiralFD` for any [`Approximant`].

use std::f64::consts::PI;

use num_complex::Complex;

use super::{
    non_negative, positive, ConditioningParameters, FdGenerator, FdRequest,
    FrequencyDomainPolarizations, FrequencySeries, GpsTime, TdGenerator, TdRequest,
    TimeDomainPolarizations, TimeSeries,
};
use crate::conditioning::fd_model_epoch;
#[cfg(feature = "fft")]
use crate::conditioning::peak_index;
use crate::constants::{MSUN_SI, MTSUN_SI};
use crate::error::{Error, Result};
use crate::imrphenomx::{IMRPhenomXAS, IMRPhenomXP, IMRPhenomXPHM};
#[cfg(feature = "fft")]
use crate::model::choose_td_model;
use crate::model::{choose_fd_model, Approximant, Binary, TimeDomainModel, WaveformOptions};
use crate::nrsur::NRSur7dq4;
#[cfg(feature = "fft")]
use crate::spintaylor::fp::fma;

type C64 = Complex<f64>;

/// `NextPow2` of LALSimulation's IMRPhenom models,
/// `pow(2, ceil(log2(n)))` for an integer `n`.
fn next_pow2(n: usize) -> usize {
    2.0f64.powf((n as f64).log2().ceil()) as usize
}

/// The polarizations on `k delta_f`, `k < n_full`, with the band `i_start
/// .. i_start + hp.len()` filled and zeros elsewhere, and the epoch `-1 /
/// delta_f` of LAL's frequency-domain models.
fn embed(
    delta_f: f64,
    n_full: usize,
    i_start: usize,
    hp: Vec<C64>,
    hc: Vec<C64>,
) -> FrequencyDomainPolarizations {
    let epoch = fd_model_epoch(delta_f);
    let fill = |band: Vec<C64>| {
        let mut data = vec![C64::new(0.0, 0.0); n_full.max(i_start + band.len())];
        data[i_start..i_start + band.len()].copy_from_slice(&band);
        data.truncate(n_full);
        FrequencySeries {
            epoch,
            delta_f,
            data,
        }
    };
    FrequencyDomainPolarizations {
        hplus: fill(hp),
        hcross: fill(hc),
    }
}

/// The conditioning parameters of a binary as LAL takes them.
fn conditioning_parameters(binary: &Binary) -> ConditioningParameters {
    ConditioningParameters {
        mass_1: binary.mass_1,
        mass_2: binary.mass_2,
        spin_1z: binary.spin_1[2],
        spin_2z: binary.spin_2[2],
    }
}

/// A frequency-domain approximant of this crate for a binary, as
/// LALSimulation 6.2.0 generates it with `longAscNodes = eccentricity = 0`
/// and the settings of `options` (multibanding is off, see
/// [`WaveformOptions`]):
///
/// - [`FdWaveform::polarizations_uniform`] is `SimInspiralChooseFDWaveform`
///   on the grid `k delta_f` (the [`FdGenerator`]);
/// - [`FdWaveform::choose_td`] is `SimInspiralChooseTDWaveform` (the
///   [`TdGenerator`]);
/// - [`FdWaveform::time_domain`] is `SimInspiralTD` and
///   [`FdWaveform::conditioned`] is `SimInspiralFD`.
///
/// The orbital phase, inclination, distance and reference frequency are
/// those of the binary (the generators take the reference frequency from
/// the request, which the conditioning may change).
#[derive(Clone, Debug, PartialEq)]
pub struct FdWaveform {
    pub approximant: Approximant,
    pub binary: Binary,
    pub options: WaveformOptions,
}

impl FdWaveform {
    /// The approximant for the binary. Returns [`Error::Unsupported`] for a
    /// time-domain approximant.
    pub fn new(approximant: Approximant, binary: Binary, options: WaveformOptions) -> Result<Self> {
        if !approximant.is_frequency_domain() {
            return Err(Error::Unsupported {
                approximant,
                what: "frequency-domain waveforms (it is a time-domain model)",
            });
        }
        Ok(Self {
            approximant,
            binary,
            options,
        })
    }

    fn invalid(&self, message: &'static str) -> Error {
        Error::InvalidParameter {
            approximant: self.approximant,
            message,
        }
    }

    /// `(h_+, h_x)` on the grid `k delta_f`, reproducing
    /// `SimInspiralChooseFDWaveform(..., delta_f, f_min, f_max, f_ref,
    /// ...)`: zero below `f_min` and above the model's cutoff, with the
    /// length and the epoch (`-1 / delta_f`) of LAL's series.
    pub fn polarizations_uniform(
        &self,
        request: &FdRequest,
    ) -> Result<FrequencyDomainPolarizations> {
        let &FdRequest {
            delta_f,
            f_min,
            f_max,
            f_ref,
        } = request;
        if !positive(delta_f) || !positive(f_min) || !non_negative(f_max) || !non_negative(f_ref) {
            return Err(
                self.invalid("delta_f and f_min must be positive, f_max and f_ref non-negative")
            );
        }
        let b = &self.binary;
        let options = &self.options;
        match self.approximant {
            Approximant::IMRPhenomXP => {
                if options.mode_array.is_some() {
                    return Err(Error::Unsupported {
                        approximant: self.approximant,
                        what: "a mode array",
                    });
                }
                let binary = Binary { f_ref, ..*b };
                let (hp, hc) = IMRPhenomXP::with_version(&binary, options.precession_version)?
                    .with_version_fallback(options.version_fallback)
                    .polarizations_uniform(delta_f, f_min, f_max)?;
                let n = hp.len();
                Ok(embed(delta_f, n, 0, hp, hc))
            }
            Approximant::IMRPhenomXPHM => {
                let binary = Binary { f_ref, ..*b };
                let mut model = IMRPhenomXPHM::with_version(&binary, options.precession_version)?
                    .with_version_fallback(options.version_fallback);
                if let Some(modes) = &options.mode_array {
                    // as crate::choose_fd_model: both signs of m' are twisted up
                    let mut positive: Vec<(u32, u32)> = modes
                        .iter()
                        .map(|&(ell, emm)| (ell, emm.unsigned_abs()))
                        .collect();
                    positive.sort_unstable();
                    positive.dedup();
                    model = model.with_mode_array(&positive)?;
                }
                let (hp, hc) = model.polarizations_uniform(delta_f, f_min, f_max)?;
                let n = hp.len();
                Ok(embed(delta_f, n, 0, hp, hc))
            }
            Approximant::TaylorF2 => {
                // XLALSimInspiralTaylorF2: m = m1_SI / MSUN_SI + m2_SI / MSUN_SI
                let m1 = b.mass_1 * MSUN_SI / MSUN_SI;
                let m2 = b.mass_2 * MSUN_SI / MSUN_SI;
                let pi_m = PI * ((m1 + m2) * MTSUN_SI);
                let v_isco = 1.0 / 6.0f64.sqrt();
                let f_isco = v_isco * v_isco * v_isco / pi_m;
                let f_max = if f_max == 0.0 { f_isco } else { f_max };
                if f_max <= f_min {
                    return Err(self.invalid("f_max must be above f_min"));
                }
                let n = (f_max / delta_f + 1.0) as usize;
                let epoch = fd_model_epoch(delta_f);
                let i_start = (f_min / delta_f).ceil() as usize;
                let frequencies: Vec<f64> = (i_start..n).map(|i| i as f64 * delta_f).collect();
                let model = choose_fd_model(self.approximant, &Binary { f_ref, ..*b }, options)?;
                let (mut hp, mut hc) = model.polarizations(&frequencies)?;
                // the phase term shift * f of XLALSimInspiralTaylorF2Core,
                // shift = 2 pi t_c for the epoch t_c = -1 / delta_f
                let shift =
                    2.0 * PI * (epoch.gps_seconds() as f64 + 1e-9 * epoch.gps_nanoseconds() as f64);
                for ((p, c), f) in hp.iter_mut().zip(hc.iter_mut()).zip(&frequencies) {
                    let (sin, cos) = (shift * f).sin_cos();
                    let factor = C64::new(cos, -sin);
                    *p *= factor;
                    *c *= factor;
                }
                Ok(embed(delta_f, n, i_start, hp, hc))
            }
            Approximant::IMRPhenomD => {
                // XLALSimIMRPhenomDGenerateFD
                let binary = Binary {
                    f_ref: if f_ref == 0.0 { f_min } else { f_ref },
                    ..*b
                };
                let model = choose_fd_model(self.approximant, &binary, options)?;
                let m_sec =
                    (b.mass_1 * MSUN_SI / MSUN_SI + b.mass_2 * MSUN_SI / MSUN_SI) * MTSUN_SI;
                let f_cut = 0.2 / m_sec;
                if f_cut <= f_min {
                    return Err(self.invalid("the cutoff frequency Mf = 0.2 is below f_min"));
                }
                let f_max_prime = if f_max != 0.0 { f_max } else { f_cut }.min(f_cut);
                if f_max_prime <= f_min {
                    return Err(self.invalid("f_max must be above f_min"));
                }
                let mut n = next_pow2((f_max_prime / delta_f) as usize) + 1;
                let i_start = (f_min / delta_f) as usize;
                let i_stop = (f_max_prime / delta_f) as usize;
                let frequencies: Vec<f64> = (i_start..i_stop).map(|i| i as f64 * delta_f).collect();
                let (hp, hc) = model.polarizations(&frequencies)?;
                if f_max_prime < f_max {
                    n = next_pow2((f_max / delta_f) as usize) + 1;
                }
                Ok(embed(delta_f, n, i_start, hp, hc))
            }
            Approximant::IMRPhenomXAS | Approximant::IMRPhenomXHM => {
                // XLALSimIMRPhenomXASGenerateFD, XLALSimIMRPhenomXHM2
                let binary = Binary {
                    f_ref: if f_ref == 0.0 { f_min } else { f_ref },
                    ..*b
                };
                let model = choose_fd_model(self.approximant, &binary, options)?;
                // the cutoff of IMRPhenomXSetWaveformVariables
                let (m1, m2, chi_1, chi_2) = (b.mass_1, b.mass_2, b.spin_1[2], b.spin_2[2]);
                let f_max_prime =
                    IMRPhenomXAS::from_components(m1, m2, chi_1, chi_2, b.luminosity_distance)?
                        .wf
                        .f_max_prime(f_max);
                let i_start = (f_min / delta_f) as usize;
                let i_stop = (f_max_prime / delta_f) as usize + 1;
                let frequencies: Vec<f64> = (i_start..i_stop).map(|i| i as f64 * delta_f).collect();
                let (hp, hc) = model.polarizations(&frequencies)?;
                let last_frequency = if f_max_prime < f_max {
                    f_max
                } else {
                    f_max_prime
                };
                let n = next_pow2((last_frequency / delta_f) as usize) + 1;
                Ok(embed(delta_f, n, i_start, hp, hc))
            }
            approximant => Err(Error::Unsupported {
                approximant,
                what: "frequency-domain waveforms",
            }),
        }
    }

    /// `SimInspiralChooseTDWaveform(..., delta_t, f_min, f_ref, ...)`, the
    /// unconditioned time-domain waveform, which TaylorF2 does not have.
    ///
    /// In LALSimulation 6.2.0 it is `XLALSimInspiralTDFromFD`
    /// ([`super::td_from_fd`]) and, for IMRPhenomD, IMRPhenomXAS and
    /// IMRPhenomXHM, the epoch is then set so that the last maximum of
    /// `h_+^2 + h_x^2` is at `t = 0` (IMRPhenomD is transformed face-on
    /// and the inclination applied in the time domain).
    #[cfg(feature = "fft")]
    pub fn choose_td(
        &self,
        delta_t: f64,
        f_min: f64,
        f_ref: f64,
    ) -> Result<TimeDomainPolarizations> {
        let f_ref = if f_ref == 0.0 && FdGenerator::zero_reference_frequency_is_f_min(self) {
            f_min
        } else {
            f_ref
        };
        let align = |h: &mut TimeDomainPolarizations| {
            let peak = peak_index(&h.hplus.data, &h.hcross.data);
            let epoch = GpsTime::from_seconds(-1.0 * delta_t * peak as f64);
            h.hplus.epoch = epoch;
            h.hcross.epoch = epoch;
        };
        match self.approximant {
            Approximant::TaylorF2 => Err(Error::Unsupported {
                approximant: self.approximant,
                what: "SimInspiralChooseTDWaveform (no time-domain generator in LALSimulation)",
            }),
            Approximant::IMRPhenomD => {
                let face_on = Self {
                    binary: Binary {
                        inclination: 0.0,
                        ..self.binary
                    },
                    ..self.clone()
                };
                let mut h = super::td_from_fd(&face_on, delta_t, f_min, f_ref)?;
                let peak = peak_index(&h.hplus.data, &h.hcross.data);
                let cfac = self.binary.inclination.cos();
                let pfac = 0.5 * fma(cfac, cfac, 1.0);
                h.hplus.data.iter_mut().for_each(|x| *x *= pfac);
                h.hcross.data.iter_mut().for_each(|x| *x *= cfac);
                let epoch = GpsTime::from_seconds(-1.0 * delta_t * peak as f64);
                h.hplus.epoch = epoch;
                h.hcross.epoch = epoch;
                Ok(h)
            }
            Approximant::IMRPhenomXAS | Approximant::IMRPhenomXHM => {
                let mut h = super::td_from_fd(self, delta_t, f_min, f_ref)?;
                align(&mut h);
                Ok(h)
            }
            _ => Ok(super::td_from_fd(self, delta_t, f_min, f_ref)?),
        }
    }

    /// `SimInspiralTD(..., delta_t, f_min, binary.f_ref, ...)`: the
    /// conditioned time-domain waveform.
    ///
    /// For all approximants but TaylorF2, LALSimulation 6.2.0 conditions
    /// the waveform of [`Self::choose_td`] as a time-domain waveform
    /// ([`super::conditioned_td`]): it is generated from a lower frequency,
    /// tapered over the extra time, high-pass filtered at `f_min` and
    /// tapered at both ends. For TaylorF2, which has no time-domain
    /// generator, `SimInspiralTD` is [`super::td_from_fd`].
    #[cfg(feature = "fft")]
    pub fn time_domain(&self, delta_t: f64, f_min: f64) -> Result<TimeDomainPolarizations> {
        if self.approximant == Approximant::TaylorF2 {
            Ok(super::td_from_fd(self, delta_t, f_min, self.binary.f_ref)?)
        } else {
            Ok(super::conditioned_td(
                self,
                delta_t,
                f_min,
                self.binary.f_ref,
            )?)
        }
    }

    /// `SimInspiralFD(..., delta_f, f_min, f_max, binary.f_ref, ...)`: the
    /// conditioned frequency-domain waveform (see [`super::conditioned_fd`]).
    pub fn conditioned(
        &self,
        delta_f: f64,
        f_min: f64,
        f_max: f64,
    ) -> Result<FrequencyDomainPolarizations> {
        Ok(super::conditioned_fd(
            self,
            delta_f,
            f_min,
            f_max,
            self.binary.f_ref,
        )?)
    }
}

impl FdGenerator for FdWaveform {
    type Error = Error;

    fn conditioning_parameters(&self) -> ConditioningParameters {
        conditioning_parameters(&self.binary)
    }

    fn zero_reference_frequency_is_f_min(&self) -> bool {
        matches!(
            self.approximant,
            Approximant::IMRPhenomXP | Approximant::IMRPhenomXPHM
        )
    }

    fn is_inspiral_only(&self) -> bool {
        self.approximant == Approximant::TaylorF2
    }

    fn generate_fd(&self, request: &FdRequest) -> Result<FrequencyDomainPolarizations> {
        self.polarizations_uniform(request)
    }
}

/// The frequency-domain approximants as `SimInspiralChooseTDWaveform`
/// generates them ([`FdWaveform::choose_td`]), which is what
/// `SimInspiralTD` conditions for all of them but TaylorF2.
#[cfg(feature = "fft")]
impl TdGenerator for FdWaveform {
    type Error = Error;

    fn conditioning_parameters(&self) -> ConditioningParameters {
        conditioning_parameters(&self.binary)
    }

    fn zero_reference_frequency_is_f_min(&self) -> bool {
        FdGenerator::zero_reference_frequency_is_f_min(self)
    }

    fn generate_td(&self, request: &TdRequest) -> Result<TimeDomainPolarizations> {
        self.choose_td(request.delta_t, request.f_min, request.f_ref)
    }
}

/// A time-domain model of this crate ([`crate::TimeDomainModel`], e.g.
/// [`NRSur7dq4`] or a model from [`crate::choose_td_model`]) with the binary
/// it was built for, as a [`TdGenerator`].
///
/// The reference frequency is fixed in the model: a request with another
/// `f_ref` than `binary.f_ref` is an error (LAL's conditioning never
/// changes it for NRSur7dq4, which allows `f_min = 0` and is only tapered).
#[derive(Clone, Copy, Debug)]
pub struct TdModel<'a, M: ?Sized> {
    pub model: &'a M,
    pub binary: Binary,
}

impl<M: TimeDomainModel + ?Sized> TdGenerator for TdModel<'_, M> {
    type Error = Error;

    fn conditioning_parameters(&self) -> ConditioningParameters {
        conditioning_parameters(&self.binary)
    }

    /// NRSur7dq4 allows `f_min = 0` (`LAL_SIM_INSPIRAL_ALLOW_ZERO_FMIN`),
    /// so LAL only tapers it.
    fn fixed_start_frequency(&self) -> bool {
        self.model.approximant() == Approximant::NRSur7dq4
    }

    fn generate_td(&self, request: &TdRequest) -> Result<TimeDomainPolarizations> {
        if request.f_ref != self.binary.f_ref {
            return Err(Error::InvalidRequest(
                "the reference frequency of a time-domain model is fixed by its binary",
            ));
        }
        let h = self.model.polarizations(request.delta_t, request.f_min)?;
        let epoch = GpsTime::from_seconds(h.epoch);
        let series = |data: Vec<f64>| TimeSeries {
            epoch,
            delta_t: h.delta_t,
            data,
        };
        Ok(TimeDomainPolarizations {
            hplus: series(h.hplus),
            hcross: series(h.hcross),
        })
    }
}

/// `SimInspiralFD(..., delta_f, f_min, f_max, f_ref, NRSur7dq4)` with the
/// binary (and reference frequency) of the model (see
/// [`super::fd_from_td`]).
#[cfg(feature = "fft")]
pub fn nrsur7dq4_fd(
    model: &NRSur7dq4<'_>,
    delta_f: f64,
    f_min: f64,
    f_max: f64,
) -> Result<FrequencyDomainPolarizations> {
    let binary = *model.binary();
    let generator = TdModel { model, binary };
    Ok(super::fd_from_td(
        &generator,
        delta_f,
        f_min,
        f_max,
        binary.f_ref,
    )?)
}

/// `SimInspiralTD(..., delta_t, f_min, f_ref, NRSur7dq4)`: the waveform of
/// `SimInspiralChooseTDWaveform` with [`crate::conditioning::wave_taper`]
/// at the start.
pub fn nrsur7dq4_td(
    model: &NRSur7dq4<'_>,
    delta_t: f64,
    f_min: f64,
) -> Result<TimeDomainPolarizations> {
    let binary = *model.binary();
    let generator = TdModel { model, binary };
    Ok(super::conditioned_td(
        &generator,
        delta_t,
        f_min,
        binary.f_ref,
    )?)
}

/// `SimInspiralTD(m1, m2, *spin_1, *spin_2, distance, inclination, phi_ref,
/// 0, 0, 0, delta_t, f_min, f_ref, params, approximant)` of LALSimulation
/// 6.2.0 for any approximant of the crate, with the `LALDict` settings in
/// `options` (see [`WaveformOptions`]; multibanding is off): the
/// conditioned time-domain waveform, transformed from the frequency domain
/// for the frequency-domain approximants ([`FdWaveform::time_domain`]).
#[cfg(feature = "fft")]
pub fn simulate_td(
    approximant: Approximant,
    binary: &Binary,
    options: &WaveformOptions,
    delta_t: f64,
    f_min: f64,
) -> Result<TimeDomainPolarizations> {
    if approximant.is_frequency_domain() {
        FdWaveform::new(approximant, *binary, options.clone())?.time_domain(delta_t, f_min)
    } else {
        let model = choose_td_model(approximant, binary, options)?;
        let generator = TdModel {
            model: &*model,
            binary: *binary,
        };
        Ok(super::conditioned_td(
            &generator,
            delta_t,
            f_min,
            binary.f_ref,
        )?)
    }
}

/// `SimInspiralFD(m1, m2, *spin_1, *spin_2, distance, inclination, phi_ref,
/// 0, 0, 0, delta_f, f_min, f_max, f_ref, params, approximant)` of
/// LALSimulation 6.2.0 for any approximant of the crate (see
/// [`simulate_td`]): the conditioned frequency-domain waveform, transformed
/// from the time domain for NRSur7dq4 ([`super::fd_from_td`]).
#[cfg(feature = "fft")]
pub fn simulate_fd(
    approximant: Approximant,
    binary: &Binary,
    options: &WaveformOptions,
    delta_f: f64,
    f_min: f64,
    f_max: f64,
) -> Result<FrequencyDomainPolarizations> {
    if approximant.is_frequency_domain() {
        FdWaveform::new(approximant, *binary, options.clone())?.conditioned(delta_f, f_min, f_max)
    } else {
        let model = choose_td_model(approximant, binary, options)?;
        let generator = TdModel {
            model: &*model,
            binary: *binary,
        };
        Ok(super::fd_from_td(
            &generator,
            delta_f,
            f_min,
            f_max,
            binary.f_ref,
        )?)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rejects_time_domain_approximants() {
        let binary = Binary {
            mass_1: 30.0,
            mass_2: 20.0,
            spin_1: [0.0; 3],
            spin_2: [0.0; 3],
            luminosity_distance: 100.0,
            inclination: 0.3,
            phi_ref: 0.0,
            f_ref: 0.0,
        };
        assert!(matches!(
            FdWaveform::new(Approximant::NRSur7dq4, binary, WaveformOptions::default()),
            Err(Error::Unsupported { .. })
        ));
        let in_plane = Binary {
            spin_1: [0.1, 0.0, 0.0],
            ..binary
        };
        let model = FdWaveform::new(
            Approximant::IMRPhenomD,
            in_plane,
            WaveformOptions::default(),
        )
        .unwrap();
        assert!(matches!(
            model.conditioned(0.0, 20.0, 1024.0),
            Err(Error::InPlaneSpins(Approximant::IMRPhenomD))
        ));
    }
}

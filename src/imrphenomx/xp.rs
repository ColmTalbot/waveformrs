//! The precessing (2,2)-mode model IMRPhenomXP, following
//! `XLALSimIMRPhenomXPFrequencySequence`, `XLALSimIMRPhenomXPGenerateFD` and
//! `IMRPhenomXPGenerateFD` in `LALSimIMRPhenomX.c`.

use std::f64::consts::PI;

use num_complex::Complex;

use crate::constants::MSUN_SI;
use crate::grid::{FrequencyGrid, ScalePowers};
use crate::imrphenomx::internals::{hz_to_mf, IMRPhenomXError, UsefulPowers, WaveformStruct};
use crate::imrphenomx::phase::phase_22;
use crate::imrphenomx::precessing_grid::{
    add_twisted_mode, coprecessing_mode, rotate_polarizations, AngleArrays, TwistArrays,
};
use crate::imrphenomx::precession::{PrecessingBinary, PrecessionStruct, PrecessionVersion};
use crate::imrphenomx::spintaylor_angles::{
    inspiral_angles_spintaylor, EulerAngles, SpinTaylorAngles,
};
use crate::imrphenomx::version_fallback::{check_version, forbid_spintaylor_fallback};
use crate::imrphenomx::wigner::twist_up_22;
use crate::imrphenomx::xas::IMRPhenomXAS;
use crate::spintaylor::fp::fma;

mod jframe;

/// `pWF->phifRef` as `IMRPhenomXPGenerateFD` evaluates it,
/// `-(phi22(MfRef) / eta + linb MfRef + lina) + 2 phi0 + pi / 4`, which LAL's
/// arm64 build contracts differently from `IMRPhenomXASGenerateFD`.
fn phif_ref_xp(xas: &IMRPhenomXAS, wf: &WaveformStruct) -> f64 {
    let inveta = 1.0 / wf.eta;
    let powers = UsefulPowers::new(wf.mf_ref);
    let phi = phase_22(wf.mf_ref, &powers, &xas.phase22);
    // the tidal phase phiTfRef vanishes for black holes
    let shift = fma(inveta, phi, xas.linb * wf.mf_ref) + xas.lina + 0.0;
    fma(2.0, wf.phi0, -shift) + std::f64::consts::FRAC_PI_4
}

/// `(h_+, h_x)` on a frequency grid.
pub type Polarizations = (Vec<Complex<f64>>, Vec<Complex<f64>>);

/// IMRPhenomXP (Pratten et al. 2021, arXiv:2004.06503): the aligned-spin
/// (2,2) mode of IMRPhenomXAS twisted up with precession Euler angles, with
/// the default LALSimulation settings (`PhenomXPConvention` 1,
/// `PhenomXPFinalSpinMod` 4) and the SpinTaylor angles
/// (`PhenomXPrecVersion` 320, Colleoni et al. 2024, arXiv:2412.16721), the
/// MSA angles (223, Chatziioannou et al. 2017, arXiv:1703.03967) or the
/// NNLO angles (102), see [`PrecessionVersion`].
///
/// # Conventions
///
/// [`IMRPhenomXP::polarizations`] reproduces the `hplus`, `hcross` of
/// `SimInspiralChooseFDWaveformSequence(phi_ref, m1, m2, *spin_1, *spin_2,
/// f_ref, distance, inclination, params, IMRPhenomXP, frequencies)` where
/// `params` sets `PhenomXPrecVersion`; [`IMRPhenomXP::polarizations_uniform`]
/// reproduces `SimInspiralChooseFDWaveform` with `deltaF > 0`. Masses are
/// in solar masses, the distance in Mpc and the spins are the dimensionless
/// spin vectors in the L frame at the reference frequency (see
/// [`PrecessingBinary`]).
///
/// # Precomputation
///
/// Everything that depends only on the parameters (the J-frame
/// construction, the precessing final spin, the NNLO or MSA angle
/// coefficients and the IMRPhenomXAS coefficients) is computed by
/// [`IMRPhenomXP::new`] when a reference frequency is given. The SpinTaylor angles depend on the
/// frequency grid in LALSimulation (they are integrated from just below the
/// first frequency and the merger-ringdown part of `gamma` is integrated
/// along the grid), as does the remnant spin they set, so they are computed
/// by [`IMRPhenomXP::polarizations`]. With `f_ref = 0` the reference
/// frequency is the first frequency of the grid and everything is computed
/// per grid.
///
/// # Fallbacks
///
/// Like LALSimulation, the SpinTaylor angles fall back to the NNLO angles
/// for in-plane spins below `1e-7` and, when the PN integration does not
/// provide the angles (it stops before reaching `f_min`, or the PN solution
/// starts above it), silently to the MSA angles (`PhenomXPrecVersion` 223).
/// The latter depends on the frequency grid. The MSA angles in turn fall
/// back to the NNLO angles when the MSA system cannot be initialised.
#[derive(Clone, Debug)]
pub struct IMRPhenomXP {
    binary: PrecessingBinary,
    version: PrecessionVersion,
    /// Whether the prescription may change as in LAL (see
    /// [`IMRPhenomXP::with_version_fallback`]).
    version_fallback: bool,
    /// Parameter-only set-up for a non-uniform grid (`f_ref > 0`).
    setup: Option<Setup>,
}

/// Waveform and precession structs (and the XAS coefficients when they do
/// not depend on the grid).
#[derive(Clone, Debug)]
struct Setup {
    wf: WaveformStruct,
    prec: PrecessionStruct,
    xas: Option<IMRPhenomXAS>,
}

impl Setup {
    /// `IMRPhenomXSetWaveformVariables` and
    /// `IMRPhenomXGetAndSetPrecessionVariables` as called by
    /// `XLALSimIMRPhenomXPFrequencySequence` / `XLALSimIMRPhenomXPGenerateFD`.
    fn new(
        binary: &PrecessingBinary,
        version: PrecessionVersion,
        delta_f: f64,
        f_min: f64,
        f_max: f64,
    ) -> Result<Self, IMRPhenomXError> {
        let (m1, m2, chi1, chi2) = binary.ordered();
        // The aligned-spin struct receives the masses and spins in the
        // original order (and orders them with its own convention).
        let mut wf = WaveformStruct::new(
            binary.mass_1,
            binary.mass_2,
            binary.spin_1[2],
            binary.spin_2[2],
            delta_f,
            binary.f_ref,
            binary.phi_ref,
            f_min,
            f_max,
            binary.luminosity_distance,
            binary.inclination,
        )?;
        // IMRPhenomXP only uses the (2,2) mode: M_MIN = M_MAX = 2.
        let prec = PrecessionStruct::new(
            &mut wf,
            m1 * MSUN_SI,
            m2 * MSUN_SI,
            chi1,
            chi2,
            version,
            2,
            2,
        )?;
        let xas = (prec.precessing_tag != 3).then(|| IMRPhenomXAS::from_waveform_struct(wf, true));
        Ok(Self { wf, prec, xas })
    }

    /// Update the frequency bounds of the waveform struct for a
    /// non-uniform grid, as `IMRPhenomXSetWaveformVariables` sets them.
    fn with_bounds(&self, f_min: f64, f_max: f64) -> Self {
        let mut setup = self.clone();
        let wf = &mut setup.wf;
        wf.f_min = f_min;
        wf.f_max = f_max;
        wf.mf_max = hz_to_mf(f_max, wf.m_tot);
        wf.f_max_prime = wf.f_max_prime(f_max);
        setup
    }

    /// The SpinTaylor angles on the grid `frequencies` and the XAS model
    /// with the final spin they set. Failures of the PN integration are
    /// returned as [`IMRPhenomXError::SpinTaylorFallbackToMsa`] and handled
    /// by [`fallback_to_msa`] (`IMRPhenomXGetAndSetPrecessionVariables`
    /// integration step and `IMRPhenomXPSpinTaylorAnglesIMR`). `fast`
    /// selects the vectorised merger-ringdown angles of the grid path.
    fn spintaylor_angles(
        &mut self,
        frequencies: &[f64],
        fast: bool,
    ) -> Result<(EulerAngles, f64), IMRPhenomXError> {
        let prec = &self.prec;
        let f_low = (self.wf.f_min - prec.integration_buffer) * 2. / prec.m_max as f64;
        if f_low.is_nan() || f_low <= 0. {
            return Err(IMRPhenomXError::InvalidFrequencies(
                "starting frequency for SpinTaylor angles must be positive",
            ));
        }
        let arrays = inspiral_angles_spintaylor(&self.wf, prec, f_low)?;
        let mut angles = SpinTaylorAngles::new(&mut self.wf, prec, &arrays)?;
        let mf_cut = self.wf.f_max_prime * self.wf.m_sec;
        let euler = angles.angles_imr(&self.wf, prec, frequencies, mf_cut, fast)?;
        self.xas = Some(IMRPhenomXAS::from_waveform_struct(self.wf, true));
        Ok((euler, angles.gamma_ref))
    }

    /// `IMRPhenomXPGenerateFD` on the grid `frequencies`, followed by the
    /// rotation of the polarisations by `2 zeta`.
    fn generate(&mut self, frequencies: &[f64]) -> Result<Polarizations, IMRPhenomXError> {
        let n = frequencies.len();
        let mut hp = vec![Complex::new(0.0, 0.0); n];
        let mut hc = vec![Complex::new(0.0, 0.0); n];
        let spintaylor = if self.prec.precessing_tag == 3 {
            Some(self.spintaylor_angles(frequencies, false)?)
        } else {
            None
        };
        let wf = &self.wf;
        let prec = &self.prec;
        let xas = self.xas.as_ref().expect("XAS coefficients");
        // Convention 1 sets phi0 = 0.
        let phif_ref = phif_ref_xp(xas, wf);
        let mf_cut = wf.f_max_prime * wf.m_sec;
        let y2 = &prec.ylm.y2;
        let eta = wf.eta;

        match &spintaylor {
            Some((angles, gamma_ref)) => {
                for (i, &f) in frequencies.iter().enumerate() {
                    let mf = wf.m_sec * f;
                    if mf <= mf_cut {
                        // LAL's powers of Mf (IMRPhenomX_Initialize_Powers)
                        let h22 = xas.h22_with_powers(mf, &UsefulPowers::new(mf), phif_ref);
                        let epsilon = -(angles.gamma[i] - gamma_ref) + prec.epsilon0;
                        (hp[i], hc[i]) =
                            twist_up_22(h22, angles.alpha[i], angles.cos_beta[i], epsilon, y2);
                    }
                }
            }
            None if prec.msa.is_some() => {
                for (i, &f) in frequencies.iter().enumerate() {
                    let mf = wf.m_sec * f;
                    if mf <= mf_cut {
                        // LAL's powers of Mf (IMRPhenomX_Initialize_Powers)
                        let h22 = xas.h22_with_powers(mf, &UsefulPowers::new(mf), phif_ref);
                        let (alpha, epsilon, cos_beta) =
                            prec.msa_angles((std::f64::consts::PI * mf).cbrt());
                        (hp[i], hc[i]) = twist_up_22(h22, alpha, cos_beta, epsilon, y2);
                    }
                }
            }
            None => {
                for (i, &f) in frequencies.iter().enumerate() {
                    let mf = wf.m_sec * f;
                    if mf <= mf_cut {
                        // LAL's powers of Mf (IMRPhenomX_Initialize_Powers)
                        let h22 = xas.h22_with_powers(mf, &UsefulPowers::new(mf), phif_ref);
                        let (alpha, epsilon, cos_beta) =
                            prec.nnlo_angles(std::f64::consts::PI * mf, eta);
                        (hp[i], hc[i]) = twist_up_22(h22, alpha, cos_beta, epsilon, y2);
                    }
                }
            }
        }

        let zeta = prec.zeta_polarization;
        if zeta.abs() > 0.0 {
            let cos_pol = (2.0 * zeta).cos();
            let sin_pol = (2.0 * zeta).sin();
            for (p, c) in hp.iter_mut().zip(hc.iter_mut()) {
                let (pp, pc) = (*p, *c);
                *p = pp * cos_pol + pc * sin_pol;
                *c = pc * cos_pol - pp * sin_pol;
            }
        }
        Ok((hp, hc))
    }

    /// [`Self::generate`] on a prepared grid, written into `hplus` and
    /// `hcross` (see [`crate::imrphenomx::precessing_grid`]).
    fn generate_grid(
        &mut self,
        grid: &FrequencyGrid,
        hplus: &mut [Complex<f64>],
        hcross: &mut [Complex<f64>],
    ) -> Result<(), IMRPhenomXError> {
        let frequencies = grid.frequencies();
        let spintaylor = if self.prec.precessing_tag == 3 {
            Some(self.spintaylor_angles(frequencies, true)?)
        } else {
            None
        };
        let wf = &self.wf;
        let prec = &self.prec;
        let xas = self.xas.as_ref().expect("XAS coefficients");
        let phif_ref = phif_ref_xp(xas, wf);
        let mf_cut = wf.f_max_prime * wf.m_sec;
        // the frequencies with Mf <= mf_cut
        let n = frequencies.partition_point(|&f| wf.m_sec * f <= mf_cut);

        let angles = match spintaylor {
            Some((mut euler, gamma_ref)) => {
                euler.alpha.truncate(n);
                euler.cos_beta.truncate(n);
                let epsilon = euler.gamma[..n]
                    .iter()
                    .map(|gamma| -(gamma - gamma_ref) + prec.epsilon0)
                    .collect();
                AngleArrays {
                    alpha: euler.alpha,
                    cos_beta: euler.cos_beta,
                    epsilon,
                }
            }
            None if prec.msa.is_some() => {
                // v = (pi Mf)^(1/3) as in the per-frequency path
                let omega: Vec<f64> = frequencies[..n]
                    .iter()
                    .map(|&f| PI * (wf.m_sec * f))
                    .collect();
                let mut angles = AngleArrays::with_capacity(n);
                prec.msa_angles_many(&omega, |alpha, epsilon, cos_beta| {
                    angles.push(alpha, cos_beta, epsilon)
                });
                angles
            }
            None => AngleArrays::nnlo(prec, wf.eta, grid, wf.m_sec, 2, n),
        };
        let twist = TwistArrays::new(&angles);

        let zero = Complex::new(0.0, 0.0);
        hplus.fill(zero);
        hcross.fill(zero);
        let mut h22 = vec![zero; n];
        let scale = ScalePowers::new(wf.m_sec);
        coprecessing_mode(&mut h22, 2, &angles.epsilon, |i| {
            let powers = grid.powers(&scale, i);
            let mf = powers.itself;
            let powers = UsefulPowers::from_grid(&powers);
            (
                xas.amplitude_22_full(mf, &powers),
                xas.phase_22_full(mf, &powers, phif_ref),
            )
        });
        add_twisted_mode(2, 2, &prec.ylm.y2, &twist, &h22, hplus, hcross);
        rotate_polarizations(prec.zeta_polarization, &mut hplus[..n], &mut hcross[..n]);
        Ok(())
    }
}

/// The first and last of the frequencies passed to the drivers, checking
/// that they are positive and strictly increasing (`None` if there are
/// none).
pub(crate) fn frequency_bounds(frequencies: &[f64]) -> Result<Option<(f64, f64)>, IMRPhenomXError> {
    let (Some(&f_min), Some(&f_max)) = (frequencies.first(), frequencies.last()) else {
        return Ok(None);
    };
    if f_min.is_nan() || f_min <= 0.0 {
        return Err(IMRPhenomXError::InvalidFrequencies(
            "frequencies must be positive",
        ));
    }
    if frequencies
        .windows(2)
        .any(|w| w[0] >= w[1] || w[0].is_nan() || w[1].is_nan())
    {
        return Err(IMRPhenomXError::InvalidFrequencies(
            "frequencies must be strictly increasing",
        ));
    }
    Ok(Some((f_min, f_max)))
}

/// `IMRPhenomXGetAndSetPrecessionVariables` falls back to the MSA angles
/// (version 223) when the SpinTaylor integration fails: redo the set-up
/// with [`PrecessionVersion::Msa`] (see [`PrecessionStruct::new`]) and
/// generate the polarisations on the same grid.
fn fallback_to_msa(
    result: Result<Polarizations, IMRPhenomXError>,
    binary: &PrecessingBinary,
    delta_f: f64,
    frequencies: &[f64],
    f_min: f64,
    f_max: f64,
) -> Result<Polarizations, IMRPhenomXError> {
    match result {
        Err(IMRPhenomXError::SpinTaylorFallbackToMsa(_)) => {
            Setup::new(binary, PrecessionVersion::Msa, delta_f, f_min, f_max)?.generate(frequencies)
        }
        result => result,
    }
}

impl IMRPhenomXP {
    /// IMRPhenomXP with the default SpinTaylor angles (`PhenomXPrecVersion`
    /// 320). Note that the default of LALSimulation 6.2.0 when no version is
    /// given is 223 ([`PrecessionVersion::Msa`]).
    ///
    /// Returns an error where LALSimulation would (masses, spins beyond the
    /// Kerr bound, negative distance or reference frequency, unphysical
    /// remnant fits).
    pub fn new(binary: &PrecessingBinary) -> Result<Self, IMRPhenomXError> {
        Self::with_version(binary, PrecessionVersion::default())
    }

    /// IMRPhenomXP with the given Euler angle prescription.
    pub fn with_version(
        binary: &PrecessingBinary,
        version: PrecessionVersion,
    ) -> Result<Self, IMRPhenomXError> {
        if binary.f_ref < 0.0 {
            return Err(IMRPhenomXError::NegativeReferenceFrequency);
        }
        // Check the parameters (and precompute) with placeholder bounds,
        // which do not enter the parameter-only quantities.
        let f_placeholder = if binary.f_ref > 0.0 {
            binary.f_ref
        } else {
            20.0
        };
        let setup = Setup::new(binary, version, 0.0, f_placeholder, 2.0 * f_placeholder)?;
        Ok(Self {
            binary: *binary,
            version,
            version_fallback: true,
            setup: (binary.f_ref > 0.0).then_some(setup),
        })
    }

    /// Whether the Euler-angle prescription may change automatically as in
    /// LALSimulation (the default, `true`): SpinTaylor to NNLO for in-plane
    /// spins below `1e-7`, SpinTaylor to MSA when the SpinTaylor angles
    /// cannot be computed, MSA to NNLO when the MSA system cannot be
    /// initialised. With `false` the requested prescription is used or
    /// the evaluation returns [`IMRPhenomXError::PrecessionVersionUnavailable`].
    pub fn with_version_fallback(mut self, allow: bool) -> Self {
        self.version_fallback = allow;
        self
    }

    /// Whether the Euler-angle prescription may change automatically (see
    /// [`IMRPhenomXP::with_version_fallback`]).
    pub fn version_fallback(&self) -> bool {
        self.version_fallback
    }

    /// Without version fallbacks, the failure of the SpinTaylor angles as an
    /// error (otherwise unchanged, for [`fallback_to_msa`]).
    fn strict<T>(&self, result: Result<T, IMRPhenomXError>) -> Result<T, IMRPhenomXError> {
        if self.version_fallback {
            result
        } else {
            forbid_spintaylor_fallback(result)
        }
    }

    /// Without version fallbacks, an error unless `setup` uses the requested
    /// prescription.
    fn check(&self, setup: &Setup) -> Result<(), IMRPhenomXError> {
        if self.version_fallback {
            Ok(())
        } else {
            check_version(&setup.prec, self.version)
        }
    }

    /// The binary parameters.
    pub fn binary(&self) -> &PrecessingBinary {
        &self.binary
    }

    /// The requested Euler angle prescription.
    pub fn version(&self) -> PrecessionVersion {
        self.version
    }

    /// The `PhenomXPrecVersion` LALSimulation uses for this binary before
    /// the frequency grid is known: 102 if NNLO angles were requested, the
    /// SpinTaylor angles were requested and the in-plane spins are below
    /// `1e-7`, or the MSA angles were requested and the MSA system cannot be
    /// initialised (only known here for `f_ref > 0`); 320 or 223 otherwise.
    /// With 320, the waveform may still fall back to 223 on a given grid.
    pub fn effective_version(&self) -> u32 {
        match &self.setup {
            Some(setup) => setup.prec.version,
            None => match self.version {
                PrecessionVersion::Nnlo => 102,
                PrecessionVersion::Msa => 223,
                PrecessionVersion::SpinTaylor => {
                    let s1 = self.binary.spin_1;
                    let s2 = self.binary.spin_2;
                    let chi_in_plane =
                        (s1[0] * s1[0] + s1[1] * s1[1] + s2[0] * s2[0] + s2[1] * s2[1]).sqrt();
                    if chi_in_plane >= 1e-7 {
                        320
                    } else {
                        102
                    }
                }
            },
        }
    }

    /// `(h_+, h_x)` at the given frequencies (Hz), reproducing
    /// `SimInspiralChooseFDWaveformSequence`: the first and last frequencies
    /// are LAL's `f_min` and `f_max`, frequencies above
    /// `min(f_max, f_cut)` with `f_cut = 0.3 / M` give zero. The frequencies
    /// must be positive and increasing.
    pub fn polarizations(&self, frequencies: &[f64]) -> Result<Polarizations, IMRPhenomXError> {
        let Some((f_min, f_max)) = frequency_bounds(frequencies)? else {
            return Ok((Vec::new(), Vec::new()));
        };
        let mut setup = self.setup(f_min, f_max)?;
        fallback_to_msa(
            self.strict(setup.generate(frequencies)),
            &self.binary,
            0.0,
            frequencies,
            f_min,
            f_max,
        )
    }

    /// [`Self::polarizations`] at the frequencies of a prepared
    /// [`FrequencyGrid`], written into `hplus` and `hcross` (each with the
    /// length of the grid).
    ///
    /// For evaluating many waveforms on the same frequencies: the Euler
    /// angles are those of [`Self::polarizations`], while the (2,2) mode
    /// uses the powers of `Mf` stored in the grid and the twisting up is
    /// done in vectorisable passes over the grid. The result agrees with
    /// [`Self::polarizations`] to rounding (~1e-13 relative), including the
    /// fallback of the SpinTaylor angles to the MSA angles.
    ///
    /// # Panics
    ///
    /// If `hplus` or `hcross` do not have the length of the grid.
    pub fn polarizations_on_grid(
        &self,
        grid: &FrequencyGrid,
        hplus: &mut [Complex<f64>],
        hcross: &mut [Complex<f64>],
    ) -> Result<(), IMRPhenomXError> {
        assert_eq!(hplus.len(), grid.len(), "output length must match the grid");
        assert_eq!(
            hcross.len(),
            grid.len(),
            "output length must match the grid"
        );
        let Some((f_min, f_max)) = frequency_bounds(grid.frequencies())? else {
            return Ok(());
        };
        match self.strict(self.setup(f_min, f_max)?.generate_grid(grid, hplus, hcross)) {
            Err(IMRPhenomXError::SpinTaylorFallbackToMsa(_)) => {
                Setup::new(&self.binary, PrecessionVersion::Msa, 0.0, f_min, f_max)?
                    .generate_grid(grid, hplus, hcross)
            }
            result => result,
        }
    }

    /// The set-up for a non-uniform grid from `f_min` to `f_max`.
    fn setup(&self, f_min: f64, f_max: f64) -> Result<Setup, IMRPhenomXError> {
        let setup = match &self.setup {
            Some(setup) => setup.with_bounds(f_min, f_max),
            None => Setup::new(&self.binary, self.version, 0.0, f_min, f_max)?,
        };
        self.check(&setup)?;
        Ok(setup)
    }

    /// `(h_+, h_x)` on the uniform grid `k delta_f`, `k = 0 .. n - 1`,
    /// reproducing `SimInspiralChooseFDWaveform(..., delta_f, f_min, f_max,
    /// f_ref, ...)` (`XLALSimIMRPhenomXPGenerateFD`): zero below `f_min` and
    /// above `min(f_max, f_cut)`, and `n - 1` the next power of two of
    /// `f_max / delta_f` (`f_max = 0` selects `f_cut = 0.3 / M`).
    pub fn polarizations_uniform(
        &self,
        delta_f: f64,
        f_min: f64,
        f_max: f64,
    ) -> Result<Polarizations, IMRPhenomXError> {
        if delta_f.is_nan() || delta_f <= 0.0 {
            return Err(IMRPhenomXError::InvalidFrequencies(
                "deltaF must be positive",
            ));
        }
        if f_min.is_nan() || f_min <= 0.0 {
            return Err(IMRPhenomXError::InvalidFrequencies(
                "f_min must be positive",
            ));
        }
        if f_max < 0.0 {
            return Err(IMRPhenomXError::InvalidFrequencies(
                "f_max must be non-negative",
            ));
        }
        if self.binary.luminosity_distance.is_nan() || self.binary.luminosity_distance <= 0.0 {
            return Err(IMRPhenomXError::NegativeDistance);
        }
        let mut setup = Setup::new(&self.binary, self.version, delta_f, f_min, f_max)?;
        self.check(&setup)?;
        let f_max_prime = setup.wf.f_max_prime;
        let i_start = (f_min / delta_f) as usize;
        let i_stop = (f_max_prime / delta_f) as usize + 1;
        let frequencies: Vec<f64> = (i_start..i_stop).map(|i| i as f64 * delta_f).collect();
        let (hp_band, hc_band) = fallback_to_msa(
            self.strict(setup.generate(&frequencies)),
            &self.binary,
            delta_f,
            &frequencies,
            f_min,
            f_max,
        )?;

        let last_freq = if f_max_prime < f_max {
            f_max
        } else {
            f_max_prime
        };
        let n_full = next_pow2((last_freq / delta_f) as usize) + 1;
        let mut hp = vec![Complex::new(0.0, 0.0); n_full.max(i_stop)];
        let mut hc = hp.clone();
        hp[i_start..i_stop].copy_from_slice(&hp_band);
        hc[i_start..i_stop].copy_from_slice(&hc_band);
        hp.truncate(n_full);
        hc.truncate(n_full);
        Ok((hp, hc))
    }

    /// The SpinTaylor Euler angles `(alpha, cos(beta), gamma)` on the grid
    /// `Mf_min + k delta_f M`, `k = 0 .. (f_max / delta_f) - (f_min /
    /// delta_f)`, reproducing `XLALSimIMRPhenomXPSpinTaylorAngles` (where
    /// `gamma` is offset to `-epsilon0` at the reference frequency). This is
    /// independent of the polarisations and meant for validation.
    pub fn spin_taylor_angles(
        binary: &PrecessingBinary,
        f_min: f64,
        f_max: f64,
        delta_f: f64,
    ) -> Result<EulerAngles, IMRPhenomXError> {
        if delta_f.is_nan() || delta_f <= 0.0 {
            return Err(IMRPhenomXError::InvalidFrequencies(
                "deltaF must be positive",
            ));
        }
        if f_min.is_nan() || f_max.is_nan() || f_min <= 0.0 || f_max <= 0.0 {
            return Err(IMRPhenomXError::InvalidFrequencies(
                "fmin and fmax must be positive",
            ));
        }
        if binary.f_ref > 0.0 && binary.f_ref < f_min {
            return Err(IMRPhenomXError::InvalidFrequencies(
                "fRef must be >= fmin or 0 to use fmin",
            ));
        }
        let i_start = (f_min / delta_f) as usize;
        let i_stop = (f_max / delta_f) as usize + 1;
        let output_length = i_stop - i_start;

        // Unlike the waveform drivers, the waveform struct is set up with the
        // ordered masses and spins, at 1 Mpc and zero inclination.
        let (m1, m2, chi1, chi2) = binary.ordered();
        let mut wf = WaveformStruct::new(
            m1,
            m2,
            chi1[2],
            chi2[2],
            delta_f,
            binary.f_ref,
            binary.phi_ref,
            f_min,
            f_max,
            1.0,
            0.0,
        )?;
        let prec = PrecessionStruct::new(
            &mut wf,
            m1 * MSUN_SI,
            m2 * MSUN_SI,
            chi1,
            chi2,
            PrecessionVersion::SpinTaylor,
            2,
            2,
        )?;
        if prec.version != 320 {
            return Err(IMRPhenomXError::SpinTaylorAngles(
                "in-plane spins too small for the SpinTaylor angles",
            ));
        }
        let f_low = (wf.f_min - prec.integration_buffer) * 2. / prec.m_max as f64;
        if f_low.is_nan() || f_low <= 0. {
            return Err(IMRPhenomXError::InvalidFrequencies(
                "starting frequency for SpinTaylor angles must be positive",
            ));
        }
        let arrays = inspiral_angles_spintaylor(&wf, &prec, f_low)?;
        let angles = SpinTaylorAngles::new(&mut wf, &prec, &arrays)?;

        let mf_min = hz_to_mf(f_min, wf.m_tot);
        if !angles.alpha_spline.contains(mf_min) {
            return Err(IMRPhenomXError::SpinTaylorAngles(
                "could not evaluate angles at fMin",
            ));
        }
        let alphamin = angles.alpha_spline.eval(mf_min);
        let cosbetamin = angles.cosbeta_spline.eval(mf_min);
        let alpha_ref = if wf.mf_ref < angles.ftrans_mrd {
            if !angles.alpha_spline.contains(wf.mf_ref) {
                return Err(IMRPhenomXError::SpinTaylorAngles(
                    "could not evaluate angles at fRef",
                ));
            }
            angles.alpha_spline.eval(wf.mf_ref)
        } else {
            angles.alpha_mrd.alpha(wf.mf_ref)
        };
        let alpha_offset = -alpha_ref + prec.alpha0;

        let mut out = EulerAngles {
            alpha: Vec::with_capacity(output_length),
            cos_beta: Vec::with_capacity(output_length),
            gamma: Vec::with_capacity(output_length),
        };
        let mut frequencies = Vec::with_capacity(output_length);
        frequencies.push(mf_min);
        out.alpha.push(alphamin - alpha_ref + prec.alpha0);
        out.cos_beta.push(cosbetamin);
        out.gamma.push(0.0);
        let mut hint = 0;
        let mut step_hint = 0;
        for i in 1..output_length {
            let mf = crate::spintaylor::fp::fma(i as f64, wf.delta_mf, mf_min);
            frequencies.push(mf);
            if mf < angles.ftrans_mrd {
                if !angles.alpha_spline.contains(mf) {
                    return Err(IMRPhenomXError::SpinTaylorAngles(
                        "interpolation of SpinTaylor angles failed",
                    ));
                }
                hint = angles.alpha_spline.locate(mf, hint);
                out.alpha
                    .push(angles.alpha_spline.eval_in(hint, mf) + alpha_offset);
                out.cos_beta.push(angles.cosbeta_spline.eval_in(hint, mf));
            } else {
                out.alpha.push(angles.alpha_mrd.alpha(mf) + alpha_offset);
                out.cos_beta.push(angles.beta_mrd.beta(mf).cos());
            }
            let deltagamma = angles.gamma_step(mf, wf.delta_mf, &mut step_hint)?;
            out.gamma.push(out.gamma[i - 1] + deltagamma);
        }
        let gamma_spline = crate::interpolate::CubicSpline::new(&frequencies, &out.gamma);
        let gamma_ref = gamma_spline.eval_unchecked(wf.mf_ref);
        for gamma in out.gamma.iter_mut() {
            *gamma = *gamma - gamma_ref - prec.epsilon0;
        }
        Ok(out)
    }
}

/// `NextPow2` of `LALSimIMRPhenomXUtilities.c`.
fn next_pow2(n: usize) -> usize {
    2.0f64.powf((n as f64).log2().ceil()) as usize
}

#[cfg(test)]
mod tests {
    use super::*;

    fn binary() -> PrecessingBinary {
        PrecessingBinary {
            mass_1: 36.0,
            mass_2: 29.0,
            spin_1: [0.4, 0.2, 0.3],
            spin_2: [-0.2, 0.3, -0.1],
            luminosity_distance: 500.0,
            inclination: 0.7,
            phi_ref: 1.2,
            f_ref: 0.0,
        }
    }

    fn grid() -> Vec<f64> {
        (0..400).map(|i| 20.0 * 1.012f64.powi(i)).collect()
    }

    #[test]
    fn default_reference_frequency_is_first_frequency() {
        let frequencies = grid();
        for version in [
            PrecessionVersion::SpinTaylor,
            PrecessionVersion::Nnlo,
            PrecessionVersion::Msa,
        ] {
            let default = IMRPhenomXP::with_version(&binary(), version)
                .unwrap()
                .polarizations(&frequencies)
                .unwrap();
            let explicit = IMRPhenomXP::with_version(
                &PrecessingBinary {
                    f_ref: frequencies[0],
                    ..binary()
                },
                version,
            )
            .unwrap()
            .polarizations(&frequencies)
            .unwrap();
            assert_eq!(default, explicit);
        }
    }

    #[test]
    fn body_labels_do_not_matter() {
        let frequencies = grid();
        let b = binary();
        let swapped = PrecessingBinary {
            mass_1: b.mass_2,
            mass_2: b.mass_1,
            spin_1: b.spin_2,
            spin_2: b.spin_1,
            ..b
        };
        for version in [
            PrecessionVersion::SpinTaylor,
            PrecessionVersion::Nnlo,
            PrecessionVersion::Msa,
        ] {
            let a = IMRPhenomXP::with_version(&b, version).unwrap();
            let s = IMRPhenomXP::with_version(&swapped, version).unwrap();
            assert_eq!(
                a.polarizations(&frequencies).unwrap(),
                s.polarizations(&frequencies).unwrap()
            );
        }
    }

    #[test]
    fn aligned_spins_fall_back_to_nnlo() {
        let aligned = PrecessingBinary {
            spin_1: [0.0, 0.0, 0.3],
            spin_2: [0.0, 0.0, -0.1],
            ..binary()
        };
        let frequencies = grid();
        let spintaylor = IMRPhenomXP::new(&aligned).unwrap();
        assert_eq!(spintaylor.effective_version(), 102);
        let nnlo = IMRPhenomXP::with_version(&aligned, PrecessionVersion::Nnlo).unwrap();
        assert_eq!(
            spintaylor.polarizations(&frequencies).unwrap(),
            nnlo.polarizations(&frequencies).unwrap()
        );
        assert_eq!(
            IMRPhenomXP::new(&binary()).unwrap().effective_version(),
            320
        );
    }

    /// For the NNLO and MSA angles nothing depends on the grid, so the
    /// uniform grid gives the same values as the non-uniform entry point.
    #[test]
    fn uniform_grid_matches_sequence_for_nnlo_and_msa() {
        for version in [PrecessionVersion::Nnlo, PrecessionVersion::Msa] {
            let model = IMRPhenomXP::with_version(
                &PrecessingBinary {
                    f_ref: 20.0,
                    ..binary()
                },
                version,
            )
            .unwrap();
            let delta_f = 0.25;
            let (hp, hc) = model.polarizations_uniform(delta_f, 20.0, 1024.0).unwrap();
            assert_eq!(hp.len(), 4097);
            let frequencies: Vec<f64> = (80..=4096).map(|i| i as f64 * delta_f).collect();
            let (hp_seq, hc_seq) = model.polarizations(&frequencies).unwrap();
            assert_eq!(&hp[80..], &hp_seq[..]);
            assert_eq!(&hc[80..], &hc_seq[..]);
            assert!(hp[..80].iter().all(|h| h.norm() == 0.0));
        }
    }

    /// Where the SpinTaylor integration does not reach the first frequency
    /// LALSimulation silently uses the MSA angles, and so do we.
    #[test]
    fn spintaylor_failure_falls_back_to_msa() {
        let heavy = PrecessingBinary {
            mass_1: 180.0,
            mass_2: 65.0,
            spin_1: [0.2, 0.1, -0.1],
            spin_2: [-0.1, -0.2, 0.1],
            luminosity_distance: 500.0,
            inclination: 0.4,
            phi_ref: 0.3,
            f_ref: 15.5,
        };
        let frequencies: Vec<f64> = (0..200).map(|i| 15.5 + i as f64).collect();
        let model = IMRPhenomXP::new(&heavy).unwrap();
        let msa = IMRPhenomXP::with_version(&heavy, PrecessionVersion::Msa).unwrap();
        assert_eq!(msa.effective_version(), 223);
        assert_eq!(
            model.polarizations(&frequencies).unwrap(),
            msa.polarizations(&frequencies).unwrap()
        );
        assert_eq!(
            model.polarizations_uniform(0.25, 15.5, 200.0).unwrap(),
            msa.polarizations_uniform(0.25, 15.5, 200.0).unwrap()
        );
        // the angles alone (as XLALSimIMRPhenomXPSpinTaylorAngles) fail
        assert!(matches!(
            IMRPhenomXP::spin_taylor_angles(&heavy, 15.5, 200.0, 0.25),
            Err(IMRPhenomXError::SpinTaylorFallbackToMsa(_))
        ));
    }

    /// A fallback through a failure that LAL reports with `XLAL_ERROR`
    /// (too few PN samples for the splines). The C function returns the MSA
    /// waveform (bit for bit equal to PrecVersion 223, checked with the C
    /// API), but it leaves `xlalErrno` set, so the SWIG Python wrapper of
    /// `SimInspiralChooseFDWaveformSequence` raises "Internal function call
    /// failed" for it.
    #[test]
    fn spintaylor_failure_with_xlal_error_falls_back_to_msa() {
        let heavy = PrecessingBinary {
            mass_1: 170.0,
            mass_2: 80.0,
            spin_1: [-0.6, 0.2, 0.5],
            spin_2: [0.3, 0.6, -0.2],
            luminosity_distance: 1000.0,
            inclination: 2.2,
            phi_ref: 4.4,
            f_ref: 35.0,
        };
        let frequencies: Vec<f64> = (0..64).map(|i| 20.0 * 1.0195f64.powi(i)).collect();
        let msa = IMRPhenomXP::with_version(&heavy, PrecessionVersion::Msa).unwrap();
        assert_eq!(
            IMRPhenomXP::new(&heavy)
                .unwrap()
                .polarizations(&frequencies)
                .unwrap(),
            msa.polarizations(&frequencies).unwrap()
        );
    }

    #[test]
    fn rejects_invalid_input() {
        let model = IMRPhenomXP::new(&binary()).unwrap();
        assert!(model.polarizations(&[30.0, 20.0, 40.0]).is_err());
        assert!(model.polarizations(&[0.0, 20.0]).is_err());
        assert_eq!(model.polarizations(&[]).unwrap(), (Vec::new(), Vec::new()));
        let spinning = PrecessingBinary {
            spin_1: [0.9, 0.5, 0.0],
            ..binary()
        };
        assert!(matches!(
            IMRPhenomXP::new(&spinning),
            Err(IMRPhenomXError::UnphysicalSpinMagnitude(_))
        ));
        let negative = PrecessingBinary {
            f_ref: -1.0,
            ..binary()
        };
        assert_eq!(
            IMRPhenomXP::new(&negative).unwrap_err(),
            IMRPhenomXError::NegativeReferenceFrequency
        );
    }
}

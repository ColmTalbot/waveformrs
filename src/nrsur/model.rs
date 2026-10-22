// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2024-2026 Colm Talbot. Parts of waveformrs are translated
// from LALSuite, GSL and fdlibm: see NOTICE.md for their authors.
//! The NRSur7dq4 waveform: `PrecessingNRSur_core` and the LAL drivers
//! `XLALSimInspiralPrecessingNRSurPolarizations`,
//! `XLALSimInspiralPrecessingNRSurModes` and `XLALPrecessingNRSurDynamics`.

use std::f64::consts::{FRAC_PI_2, PI};

use num_complex::Complex;

use crate::complex::cexp_i;
use crate::constants::{C_SI, MSUN_SI, MTSUN_SI, PC_SI};
use crate::imrphenomx::PrecessingBinary;
use crate::nrsur::data::{mode_index, DataPiece, NRSur7dq4Data, LMAX, N_MODES};
use crate::nrsur::dynamics::{Dynamics, State, N_STATE};
use crate::nrsur::fit::MassRatioParameters;
use crate::nrsur::frames::coorbital_to_inertial;
use crate::nrsur::spline::{cubic, FixedKnots};
use crate::nrsur::NRSurError;
use crate::spintaylor::fp::{fma, fms, fnms};

/// `h_+` and `h_x` sampled uniformly in time, as the `REAL8TimeSeries` of
/// `SimInspiralChooseTDWaveform`.
#[derive(Clone, Debug, Default)]
pub struct TimeSeriesPolarizations {
    /// GPS time of the first sample relative to the peak of the waveform
    /// (`(*hplus)->epoch`, rounded to nanoseconds as LAL's `LIGOTimeGPS`).
    pub epoch: f64,
    /// Sampling interval in seconds.
    pub delta_t: f64,
    pub hplus: Vec<f64>,
    pub hcross: Vec<f64>,
}

/// One inertial-frame mode `h_{ell m}`.
#[derive(Clone, Debug, Default)]
pub struct Mode {
    pub ell: u32,
    pub m: i32,
    pub data: Vec<Complex<f64>>,
}

/// The modes `h_{ell m}` sampled uniformly in time, as the
/// `SphHarmTimeSeries` of `SimInspiralChooseTDModes`.
#[derive(Clone, Debug, Default)]
pub struct TimeSeriesModes {
    /// GPS time of the first sample (as for [`TimeSeriesPolarizations`]).
    pub epoch: f64,
    /// Sampling interval in seconds.
    pub delta_t: f64,
    /// The modes, ordered by `ell` then `m`.
    pub modes: Vec<Mode>,
}

impl TimeSeriesModes {
    /// The mode `(ell, m)`, if present.
    pub fn mode(&self, ell: u32, m: i32) -> Option<&[Complex<f64>]> {
        self.modes
            .iter()
            .find(|mode| mode.ell == ell && mode.m == m)
            .map(|mode| mode.data.as_slice())
    }
}

/// The surrogate dynamics on the dynamics nodes (`XLALPrecessingNRSurDynamics`).
#[derive(Clone, Debug, Default)]
pub struct NRSurDynamics {
    /// The dynamics node times in units of the total mass.
    pub t: Vec<f64>,
    /// The coprecessing-frame quaternion.
    pub quaternion: [Vec<f64>; 4],
    /// The orbital phase in the coprecessing frame.
    pub orbital_phase: Vec<f64>,
    /// The spin of the heavier black hole in the coprecessing frame.
    pub chi_a: [Vec<f64>; 3],
    /// The spin of the lighter black hole in the coprecessing frame.
    pub chi_b: [Vec<f64>; 3],
}

/// Reusable buffers for evaluating NRSur7dq4 many times without
/// allocating. Create one per thread with [`NRSurWorkspace::new`] and pass it
/// to [`NRSur7dq4::polarizations_with`] or [`NRSur7dq4::modes_with`]; it is
/// sized on first use.
///
/// The workspace also keeps the inertial-frame modes of the last binary.
/// When the next call has the same intrinsic parameters (masses, spins,
/// `f_ref`, `f_min`, mode array and extrapolation setting; compared bit for
/// bit), only the spherical harmonics, the distance scaling and the
/// resampling onto the output times are redone, e.g. when a sampler varies
/// only the inclination, phase, distance or sampling interval.
#[derive(Clone, Debug, Default)]
pub struct NRSurWorkspace {
    n_coorb: usize,
    /// The intrinsic parameters whose inertial modes are in `h_re`, `h_im`,
    /// with the labelled binary and the start of the output.
    cached: Option<(IntrinsicKey, Labelled, f64)>,
    /// ODE state at the dynamics nodes.
    dynamics: Vec<State>,
    /// The orbital frequency fit at the nodes for the reference state.
    omega_cache: Vec<f64>,
    /// The dynamics on the coorbital times: quaternion (4), orbital phase,
    /// chiA (3) and chiB (3) in the coorbital frame, `[component][time]`.
    series: Vec<f64>,
    /// cos(m phi), sin(m phi) for m = 1..=4.
    cos_mphi: Vec<[f64; LMAX]>,
    sin_mphi: Vec<[f64; LMAX]>,
    /// Spline scratch (values, second derivatives, substitution).
    spline_y: Vec<f64>,
    spline_c: Vec<f64>,
    spline_z: Vec<f64>,
    /// Data piece evaluation: fits at the empirical nodes and the result.
    nodes: Vec<f64>,
    piece: Vec<f64>,
    /// Coorbital and inertial modes, `[mode][time]`.
    h_coorb_re: Vec<f64>,
    h_coorb_im: Vec<f64>,
    h_re: Vec<f64>,
    h_im: Vec<f64>,
    /// `h_+`, `h_x` on the coorbital times and their spline coefficients.
    hp_model: Vec<f64>,
    hc_model: Vec<f64>,
    c_plus: Vec<f64>,
    c_cross: Vec<f64>,
    /// The output times located in the coorbital times, and the spline
    /// coefficients of the intervals they fall in.
    grid: OutputGrid,
    coef_b: Vec<f64>,
    coef_d: Vec<f64>,
}

/// The length (in units of the total mass) of the taper at the end of the
/// surrogate with [`NRSur7dq4::with_continuous`].
pub const END_TAPER: f64 = 10.0;

/// What the inertial modes depend on.
#[derive(Clone, Copy, Debug, PartialEq)]
struct IntrinsicKey {
    data: usize,
    bits: [u64; 10],
    unlimited_extrapolation: bool,
    continuous: bool,
    active: [bool; N_MODES],
}

impl NRSurWorkspace {
    pub fn new() -> Self {
        Self::default()
    }

    fn prepare(&mut self, data: &NRSur7dq4Data) {
        let n_ds = data.t_ds.len();
        let n = data.t_coorb.len();
        if self.n_coorb == n && self.dynamics.len() == n_ds {
            return;
        }
        self.cached = None;
        self.n_coorb = n;
        self.dynamics = vec![[0.0; N_STATE]; n_ds];
        self.omega_cache = vec![f64::NAN; n_ds];
        self.series = vec![0.0; N_STATE * n];
        self.cos_mphi = vec![[0.0; LMAX]; n];
        self.sin_mphi = vec![[0.0; LMAX]; n];
        let m = n.max(n_ds);
        self.spline_y = vec![0.0; m];
        self.spline_c = vec![0.0; m];
        self.spline_z = vec![0.0; m];
        self.piece = vec![0.0; n];
        for buffer in [
            &mut self.h_coorb_re,
            &mut self.h_coorb_im,
            &mut self.h_re,
            &mut self.h_im,
        ] {
            *buffer = vec![0.0; N_MODES * n];
        }
        for buffer in [
            &mut self.hp_model,
            &mut self.hc_model,
            &mut self.c_plus,
            &mut self.c_cross,
            &mut self.coef_b,
            &mut self.coef_d,
        ] {
            *buffer = vec![0.0; n];
        }
    }
}

/// The NRSur7dq4 surrogate for one binary.
///
/// # Conventions
///
/// [`NRSur7dq4::polarizations`] reproduces `SimInspiralChooseTDWaveform(m1,
/// m2, *spin_1, *spin_2, distance, inclination, phi_ref, 0, 0, 0, delta_t,
/// f_min, f_ref, None, NRSur7dq4)` and [`NRSur7dq4::modes`]
/// `SimInspiralChooseTDModes(0, delta_t, m1, m2, *spin_1, *spin_2, f_min,
/// f_ref, distance, None, 4, NRSur7dq4)`, with masses in solar masses (LAL
/// receives `m * MSUN_SI`), the distance in Mpc (`d * 1e6 * PC_SI`) and the
/// dimensionless spins in the LAL source frame at the reference frequency.
///
/// * The binary is set up where the orbital frequency of the coprecessing
///   frame waveform is `pi f_ref` (`f_ref = 0` uses `f_min`, `f_min = 0`
///   too selects the start of the surrogate). There the spins are given in
///   the coorbital frame: `z` along the principal axis of the radiation,
///   `x` from the lighter to the heavier black hole.
/// * If `mass_1 < mass_2` the labels are swapped (with the in-plane spins
///   rotated by pi) and the odd-`m` modes negated, as in LAL.
/// * The waveform starts where the orbital frequency (of the reference
///   spins) is `pi f_min`, or at the start of the surrogate for `f_min = 0`;
///   a non-zero `f_min` below the surrogate's starting frequency
///   ([`NRSur7dq4::start_frequency`]) is an error. The epoch is the time of
///   the first sample relative to the peak.
/// * The polarizations are evaluated at `(inclination, pi/2 - phi_ref)`.
#[derive(Clone, Debug)]
pub struct NRSur7dq4<'a> {
    data: &'a NRSur7dq4Data,
    binary: PrecessingBinary,
    unlimited_extrapolation: bool,
    /// See [`NRSur7dq4::with_continuous`].
    continuous: bool,
    /// The active coorbital modes, by [`mode_index`].
    active: [bool; N_MODES],
}

/// The masses and spins after `PrecessingNRSur_switch_labels_if_needed`.
#[derive(Clone, Copy, Debug)]
struct Labelled {
    switched: bool,
    m1_si: f64,
    m2_si: f64,
    chi_a: [f64; 3],
    chi_b: [f64; 3],
    q: f64,
    m_tot_sec: f64,
}

/// The result of the core evaluation: the inertial modes are in the
/// workspace.
struct Core {
    labels: Labelled,
    /// Start of the output in units of the total mass.
    t0: f64,
}

/// `XLALGPSAdd(&zero, t)` followed by `XLALGPSGetREAL8`: `t` rounded to
/// nanoseconds.
fn gps_epoch(t: f64) -> f64 {
    let seconds = t.floor() as i32;
    let nanoseconds = ((t - seconds as f64) * 1e9).round_ties_even() as i32;
    let total = 1_000_000_000i64 * seconds as i64 + nanoseconds as i64;
    let (s, ns) = (total / 1_000_000_000, total % 1_000_000_000);
    s as f64 + ns as f64 / 1e9
}

/// `XLALSpinWeightedSphericalHarmonic(theta, phi, -2, ell, m)` for `ell <=
/// 4`.
pub(crate) fn spin_weighted_spherical_harmonic(
    theta: f64,
    phi: f64,
    ell: usize,
    m: i32,
) -> Complex<f64> {
    let half = theta / 2.0;
    let sin = f64::sin;
    let cos = f64::cos;
    let fac = match (ell, m) {
        (2, -2) => {
            let c = cos(theta);
            (5.0 / (64.0 * PI)).sqrt() * (1.0 - c) * (1.0 - c)
        }
        (2, -1) => {
            let (s, c) = f64::sin_cos(theta);
            (5.0 / (16.0 * PI)).sqrt() * s * (1.0 - c)
        }
        (2, 0) => {
            let s = sin(theta);
            (15.0 / (32.0 * PI)).sqrt() * s * s
        }
        (2, 1) => {
            let (s, c) = f64::sin_cos(theta);
            (5.0 / (16.0 * PI)).sqrt() * s * (1.0 + c)
        }
        (2, 2) => {
            let c = cos(theta);
            (5.0 / (64.0 * PI)).sqrt() * (1.0 + c) * (1.0 + c)
        }
        (3, -3) => {
            let (sh, ch) = f64::sin_cos(half);
            (21.0 / (2.0 * PI)).sqrt() * ch * sh.powf(5.0)
        }
        (3, -2) => (7.0 / (4.0 * PI)).sqrt() * fma(3.0, cos(theta), 2.0) * sin(half).powf(4.0),
        (3, -1) => {
            let (s, s2, s3) = (sin(theta), sin(2.0 * theta), sin(3.0 * theta));
            (35.0 / (2.0 * PI)).sqrt() * fma(-3.0, s3, fma(4.0, s2, s)) / 32.0
        }
        (3, 0) => {
            let (s, c) = f64::sin_cos(theta);
            ((105.0 / (2.0 * PI)).sqrt() * c * (s * s)) / 4.0
        }
        (3, 1) => {
            let (s, s2, s3) = (sin(theta), sin(2.0 * theta), sin(3.0 * theta));
            -(35.0 / (2.0 * PI)).sqrt() * fma(-3.0, s3, fma(-4.0, s2, s)) / 32.0
        }
        (3, 2) => (7.0 / PI).sqrt() * cos(half).powf(4.0) * fma(3.0, cos(theta), -2.0) / 2.0,
        (3, 3) => {
            let (sh, ch) = f64::sin_cos(half);
            -(21.0 / (2.0 * PI)).sqrt() * ch.powf(5.0) * sh
        }
        (4, -4) => {
            let (sh, ch) = f64::sin_cos(half);
            3.0 * (7.0 / PI).sqrt() * (ch * ch) * sh.powf(6.0)
        }
        (4, -3) => {
            let (sh, ch) = f64::sin_cos(half);
            3.0 * (7.0 / (2.0 * PI)).sqrt() * ch * fma(2.0, cos(theta), 1.0) * sh.powf(5.0)
        }
        (4, -2) => {
            let (c, c2, sh) = (cos(theta), cos(2.0 * theta), sin(half));
            (3.0 * fma(7.0, c2, fma(14.0, c, 9.0)) * sh.powf(4.0)) / (4.0 * PI.sqrt())
        }
        (4, -1) => {
            let (s, s2) = (sin(theta), sin(2.0 * theta));
            let (s3, s4) = (sin(3.0 * theta), sin(4.0 * theta));
            (3.0 * fma(-7.0, s4, fma(7.0, s3, fma(3.0, s, 2.0 * s2)))) / (32.0 * (2.0 * PI).sqrt())
        }
        (4, 0) => {
            let (c2, s) = (cos(2.0 * theta), sin(theta));
            (3.0 * (5.0 / (2.0 * PI)).sqrt() * fma(7.0, c2, 5.0) * (s * s)) / 16.0
        }
        (4, 1) => {
            let (s, s2) = (sin(theta), sin(2.0 * theta));
            let (s3, s4) = (sin(3.0 * theta), sin(4.0 * theta));
            (3.0 * fma(7.0, s4, fma(7.0, s3, fms(3.0, s, 2.0 * s2)))) / (32.0 * (2.0 * PI).sqrt())
        }
        (4, 2) => {
            let (ch, c, c2) = (cos(half), cos(theta), cos(2.0 * theta));
            (3.0 * ch.powf(4.0) * fma(7.0, c2, fma(-14.0, c, 9.0))) / (4.0 * PI.sqrt())
        }
        (4, 3) => {
            let (sh, ch) = f64::sin_cos(half);
            -3.0 * (7.0 / (2.0 * PI)).sqrt() * ch.powf(5.0) * fma(2.0, cos(theta), -1.0) * sh
        }
        (4, 4) => {
            let (sh, ch) = f64::sin_cos(half);
            3.0 * (7.0 / PI).sqrt() * ch.powf(6.0) * (sh * sh)
        }
        _ => panic!("spin-weighted spherical harmonic ({ell}, {m}) not implemented"),
    };
    if m != 0 {
        let z = cexp_i(m as f64 * phi);
        Complex::new(z.re * fac, z.im * fac)
    } else {
        Complex::new(fac, 0.0)
    }
}

impl<'a> NRSur7dq4<'a> {
    /// The surrogate for `binary` (masses in solar masses, distance in Mpc,
    /// spins at `binary.f_ref`).
    pub fn new(data: &'a NRSur7dq4Data, binary: &PrecessingBinary) -> Result<Self, NRSurError> {
        for m in [binary.mass_1, binary.mass_2] {
            if !(m > 0.0 && m.is_finite()) {
                return Err(NRSurError::InvalidMass);
            }
        }
        Ok(Self {
            data,
            binary: *binary,
            unlimited_extrapolation: false,
            continuous: false,
            active: [true; N_MODES],
        })
    }

    /// Make the waveform change continuously with the parameters (off by
    /// default, reproducing LAL). LAL's surrogate has two discontinuities:
    ///
    /// * The dynamics are integrated from the node nearest to the
    ///   reference frequency, and the integration error depends on that
    ///   node, so the waveform jumps (by ~1e-4) when the reference time
    ///   crosses the middle of an interval between nodes. With this option
    ///   the dynamics are integrated from both nodes of the interval and
    ///   blended smoothly between them: they are LAL's when the reference
    ///   time is at a node.
    /// * The waveform ends abruptly at the end of the surrogate, ~100 M
    ///   after the peak, so a sample jumps to zero when the end crosses it.
    ///   With this option the last [`END_TAPER`] M are tapered to zero.
    /// * The reference time is interpolated linearly between nodes, so it
    ///   has a kink as it crosses one. With this option the corners are
    ///   rounded off.
    /// * The epoch of [`Self::polarizations_with`] and
    ///   [`Self::modes_with`] is rounded to the nanosecond (as LAL's
    ///   `LIGOTimeGPS`), so it moves in steps as the masses change. With
    ///   this option it is exact.
    pub fn with_continuous(mut self, continuous: bool) -> Self {
        self.continuous = continuous;
        self
    }

    /// Allow extrapolation beyond `q = 6.01` and `|chi| = 1`
    /// (`unlimited_extrapolation` in the LALDict). Use at your own risk.
    pub fn with_unlimited_extrapolation(mut self, unlimited: bool) -> Self {
        self.unlimited_extrapolation = unlimited;
        self
    }

    /// Evaluate only the given coorbital-frame modes `(ell, m)`, `2 <= ell
    /// <= 4` (LAL's `ModeArray`). As in LAL, `(ell, m)` and `(ell, -m)` are
    /// evaluated together if either is listed, and the inertial modes of an
    /// `ell` none of whose coorbital modes is listed are left out.
    pub fn with_mode_array(mut self, modes: &[(u32, i32)]) -> Result<Self, NRSurError> {
        let mut active = [false; N_MODES];
        for &(ell, m) in modes {
            if !(2..=LMAX as u32).contains(&ell) || m.unsigned_abs() > ell {
                return Err(NRSurError::InvalidModeArray);
            }
            active[mode_index(ell as usize, m)] = true;
        }
        self.active = active;
        Ok(self)
    }

    pub fn binary(&self) -> &PrecessingBinary {
        &self.binary
    }

    pub fn data(&self) -> &'a NRSur7dq4Data {
        self.data
    }

    /// The epoch of the output for its first sample at `t` seconds from
    /// the peak: rounded to the nanosecond as LAL's `LIGOTimeGPS`, or
    /// exact with [`Self::with_continuous`] (the rounding makes a time
    /// shift by the epoch step as the masses change).
    fn epoch(&self, t: f64) -> f64 {
        if self.continuous {
            t
        } else {
            gps_epoch(t)
        }
    }

    fn intrinsic_key(&self, f_min: f64) -> IntrinsicKey {
        let b = &self.binary;
        let values = [
            b.mass_1,
            b.mass_2,
            b.spin_1[0],
            b.spin_1[1],
            b.spin_1[2],
            b.spin_2[0],
            b.spin_2[1],
            b.spin_2[2],
            b.f_ref,
            f_min,
        ];
        IntrinsicKey {
            data: self.data as *const NRSur7dq4Data as usize,
            bits: values.map(f64::to_bits),
            unlimited_extrapolation: self.unlimited_extrapolation,
            continuous: self.continuous,
            active: self.active,
        }
    }

    fn labelled(&self) -> Labelled {
        let b = &self.binary;
        let (mut m1, mut m2) = (b.mass_1 * MSUN_SI, b.mass_2 * MSUN_SI);
        let (mut s1, mut s2) = (b.spin_1, b.spin_2);
        let switched = m1 < m2;
        if switched {
            std::mem::swap(&mut m1, &mut m2);
            let (a, c) = (s1, s2);
            s1 = [-c[0], -c[1], c[2]];
            s2 = [-a[0], -a[1], a[2]];
        }
        let mass_a = m1 / MSUN_SI;
        let mass_b = m2 / MSUN_SI;
        let m_tot = mass_a + mass_b;
        Labelled {
            switched,
            m1_si: m1,
            m2_si: m2,
            chi_a: s1,
            chi_b: s2,
            q: mass_a / mass_b,
            m_tot_sec: m_tot * MTSUN_SI,
        }
    }

    /// `PrecessingNRSur_StartFrequency`: the lowest `f_min` (in Hz) for this
    /// binary, where the orbital frequency of the first dynamics node (with
    /// the reference spins) is reached. `f_min` must be zero or at least this.
    pub fn start_frequency(&self) -> f64 {
        let l = self.labelled();
        self.start_frequency_of(&l)
    }

    fn start_frequency_of(&self, l: &Labelled) -> f64 {
        let m_tot = (l.m1_si + l.m2_si) / MSUN_SI;
        let m_tot_sec = m_tot * MTSUN_SI;
        let q = l.m1_si / l.m2_si;
        let y0 = reference_state(&l.chi_a, &l.chi_b, &[1.0, 0.0, 0.0, 0.0], 0.0);
        let omega = Dynamics::new(self.data, q).omega(0, &y0);
        omega / (PI * m_tot_sec)
    }

    /// `PrecessingNRSur_core` up to the inertial modes (in `ws.h_re`,
    /// `ws.h_im`) and the start of the output (`t0`), with LAL's checks in
    /// LAL's order.
    fn core(&self, ws: &mut NRSurWorkspace, f_min: f64, delta_t: f64) -> Result<Core, NRSurError> {
        if !(delta_t > 0.0 && delta_t.is_finite()) {
            return Err(NRSurError::InvalidTimeStep(delta_t));
        }
        ws.prepare(self.data);
        let key = self.intrinsic_key(f_min);
        if let Some((cached, labels, t0)) = ws.cached {
            if cached == key {
                return Ok(Core { labels, t0 });
            }
        }
        ws.cached = None;
        let data = self.data;
        let l = self.labelled();
        // get_dimless_omega
        let omega_min = (l.m_tot_sec * f_min) * PI;
        let f_ref = self.binary.f_ref;
        let omega_ref = if f_ref == 0.0 {
            omega_min
        } else {
            (l.m_tot_sec * f_ref) * PI
        };
        if omega_ref + 1e-13 < omega_min {
            return Err(NRSurError::ReferenceBelowMinimum { f_ref, f_min });
        }

        // The inertial frame is aligned with the coorbital frame at the
        // reference point.
        let init_quat = [1.0, 0.0, 0.0, 0.0];
        let dynamics = Dynamics::new(data, l.q);
        ws.omega_cache.fill(f64::NAN);
        dynamics.integrate(
            &l.chi_a,
            &l.chi_b,
            omega_ref,
            0.0,
            &init_quat,
            self.unlimited_extrapolation,
            self.continuous,
            &mut ws.omega_cache,
            &mut ws.dynamics,
        )?;

        // The start of the output: checked after the modes in LAL, which
        // cannot fail once the dynamics are integrated.
        let start_freq = self.start_frequency_of(&l);
        let t0 = if f_min >= start_freq {
            let y0 = reference_state(&l.chi_a, &l.chi_b, &init_quat, 0.0);
            dynamics.t_ref(omega_min, &y0, &mut ws.omega_cache, self.continuous)?
        } else if f_min > 0.0 {
            return Err(NRSurError::StartFrequency {
                f_min,
                f_start: start_freq,
            });
        } else {
            data.t_coorb.x[0]
        };

        self.coorbital_dynamics(ws, &l);
        self.coorbital_modes(ws, l.q);
        let n = ws.n_coorb;
        let quat = {
            let s = &ws.series;
            [&s[0..n], &s[n..2 * n], &s[2 * n..3 * n], &s[3 * n..4 * n]]
        };
        coorbital_to_inertial(
            &data.wigner,
            n,
            quat,
            &ws.cos_mphi,
            &ws.sin_mphi,
            &ws.h_coorb_re,
            &ws.h_coorb_im,
            &mut ws.h_re,
            &mut ws.h_im,
        );
        ws.cached = Some((key, l, t0));
        Ok(Core { labels: l, t0 })
    }

    /// Interpolate the dynamics onto the coorbital times, normalise, and
    /// rotate the spins into the coorbital frame.
    fn coorbital_dynamics(&self, ws: &mut NRSurWorkspace, l: &Labelled) {
        let data = self.data;
        let knots = &data.t_ds;
        let n_ds = knots.len();
        let n = ws.n_coorb;
        let points = &data.coorb_in_ds;
        for comp in 0..N_STATE {
            let y = &mut ws.spline_y[..n_ds];
            for (y, state) in y.iter_mut().zip(&ws.dynamics) {
                *y = state[comp];
            }
            let c = &mut ws.spline_c[..n_ds];
            knots.second_derivatives(y, c, &mut ws.spline_z);
            let out = &mut ws.series[comp * n..(comp + 1) * n];
            // coefficients of each interval, computed when first needed
            let mut current = usize::MAX;
            let (mut b, mut d) = (0.0, 0.0);
            for (k, out) in out.iter_mut().enumerate() {
                let i = points.index[k];
                if i != current {
                    (b, d) = knots.interval(i, y, c);
                    current = i;
                }
                *out = cubic(y[i], b, c[i], d, points.delx[k]);
            }
        }

        // PrecessingNRSur_normalize_results
        let norm = |v: &[f64; 3]| fma(v[2], v[2], fma(v[0], v[0], v[1] * v[1])).sqrt();
        let (norm_a, norm_b) = (norm(&l.chi_a), norm(&l.chi_b));
        let (quat, rest) = ws.series.split_at_mut(4 * n);
        let (phi, spins) = rest.split_at_mut(n);
        let (chi_a, chi_b) = spins.split_at_mut(3 * n);
        for (chi, norm) in [(&mut *chi_a, norm_a), (&mut *chi_b, norm_b)] {
            if norm > 0.0 {
                let (x, rest) = chi.split_at_mut(n);
                let (y, z) = rest.split_at_mut(n);
                for ((x, y), z) in x.iter_mut().zip(y.iter_mut()).zip(z.iter_mut()) {
                    let n_chi = fma(*z, *z, fma(*x, *x, *y * *y)).sqrt();
                    *x *= norm / n_chi;
                    *y *= norm / n_chi;
                    *z *= norm / n_chi;
                }
            }
        }
        {
            let (q0, rest) = quat.split_at_mut(n);
            let (qx, rest) = rest.split_at_mut(n);
            let (qy, qz) = rest.split_at_mut(n);
            for t in 0..n {
                let n_q = fma(
                    qz[t],
                    qz[t],
                    fma(qy[t], qy[t], fma(q0[t], q0[t], qx[t] * qx[t])),
                )
                .sqrt();
                q0[t] /= n_q;
                qx[t] /= n_q;
                qy[t] /= n_q;
                qz[t] /= n_q;
            }
        }

        // PrecessingNRSur_rotate_spins, and cos(m phi), sin(m phi) for the
        // transformation to the coprecessing frame.
        for t in 0..n {
            let mut cos_m = [0.0; LMAX];
            let mut sin_m = [0.0; LMAX];
            for (k, (c, s)) in cos_m.iter_mut().zip(sin_m.iter_mut()).enumerate() {
                let arg = (k + 1) as f64 * phi[t];
                *c = f64::cos(arg);
                *s = f64::sin(arg);
            }
            let (sp, cp) = f64::sin_cos(phi[t]);
            for chi in [&mut *chi_a, &mut *chi_b] {
                let tmp = chi[t];
                let y = chi[n + t];
                chi[t] = fma(tmp, cp, y * sp);
                chi[n + t] = fma(-tmp, sp, y * cp);
            }
            ws.cos_mphi[t] = cos_m;
            ws.sin_mphi[t] = sin_m;
        }
    }

    /// Evaluate the coorbital modes (the empirical interpolants) into
    /// `ws.h_coorb_re`, `ws.h_coorb_im`.
    fn coorbital_modes(&self, ws: &mut NRSurWorkspace, q: f64) {
        let n = ws.n_coorb;
        let params = MassRatioParameters::nrsur7dq4(q);
        // modes that are not evaluated stay zero
        if self.active.iter().any(|&a| !a) {
            ws.h_coorb_re.fill(0.0);
            ws.h_coorb_im.fill(0.0);
        }
        let NRSurWorkspace {
            series,
            nodes,
            piece,
            h_coorb_re,
            h_coorb_im,
            ..
        } = ws;
        let spins = &series[5 * n..11 * n];
        let eval = |data: &DataPiece, nodes: &mut Vec<f64>, piece: &mut [f64]| {
            eval_data_piece(data, &params, q, spins, n, nodes, piece);
        };
        for (ell_index, ell_data) in self.data.modes.iter().enumerate() {
            let ell = ell_index + 2;
            let i0 = mode_index(ell, 0);
            if self.active[i0] {
                eval(&ell_data.m0_real, nodes, piece);
                h_coorb_re[i0 * n..(i0 + 1) * n].copy_from_slice(&piece[..n]);
                eval(&ell_data.m0_imag, nodes, piece);
                h_coorb_im[i0 * n..(i0 + 1) * n].copy_from_slice(&piece[..n]);
            }
            for m in 1..=ell {
                let (ip, im) = (mode_index(ell, m as i32), mode_index(ell, -(m as i32)));
                if !self.active[ip] && !self.active[im] {
                    continue;
                }
                // h^{l,m} = X_+ + X_-, h^{l,-m} = (X_+ - X_-)^*
                let (re_p, re_m) = two_rows(h_coorb_re, ip, im, n);
                eval(&ell_data.x_real_plus[m - 1], nodes, piece);
                re_p.copy_from_slice(&piece[..n]);
                re_m.copy_from_slice(&piece[..n]);
                eval(&ell_data.x_real_minus[m - 1], nodes, piece);
                for t in 0..n {
                    re_p[t] += piece[t];
                    re_m[t] -= piece[t];
                }
                let (im_p, im_m) = two_rows(h_coorb_im, ip, im, n);
                eval(&ell_data.x_imag_plus[m - 1], nodes, piece);
                for t in 0..n {
                    im_p[t] = piece[t];
                    im_m[t] = -piece[t];
                }
                eval(&ell_data.x_imag_minus[m - 1], nodes, piece);
                for t in 0..n {
                    im_p[t] += piece[t];
                    im_m[t] += piece[t];
                }
            }
        }
    }

    /// Whether any coorbital mode of `ell` is active (otherwise LAL skips
    /// all inertial modes of that `ell`).
    fn ell_active(&self, ell: usize) -> bool {
        (-(ell as i32)..=ell as i32).any(|m| self.active[mode_index(ell, m)])
    }

    /// The output sampling interval in units of the total mass and the
    /// number of samples, `ceil((t_end - t0) / dt)`.
    fn output_grid(&self, core: &Core, delta_t: f64) -> (f64, usize) {
        let dt = delta_t / core.labels.m_tot_sec;
        let tf = self.data.t_coorb.x[self.data.t_coorb.len() - 1];
        let nt = ((tf - core.t0) / dt).ceil();
        (dt, if nt > 0.0 { nt as usize } else { 0 })
    }

    /// `h_+` and `h_x` (`XLALSimInspiralPrecessingNRSurPolarizations`).
    pub fn polarizations(
        &self,
        delta_t: f64,
        f_min: f64,
    ) -> Result<TimeSeriesPolarizations, NRSurError> {
        let mut out = TimeSeriesPolarizations::default();
        self.polarizations_with(&mut NRSurWorkspace::new(), delta_t, f_min, &mut out)?;
        Ok(out)
    }

    /// [`Self::polarizations`] reusing the buffers of `ws` and `out`.
    pub fn polarizations_with(
        &self,
        ws: &mut NRSurWorkspace,
        delta_t: f64,
        f_min: f64,
        out: &mut TimeSeriesPolarizations,
    ) -> Result<(), NRSurError> {
        let core = self.core(ws, f_min, delta_t)?;
        self.model_polarizations(ws, &core);
        let (dt, nt) = self.output_grid(&core, delta_t);
        out.epoch = self.epoch(core.labels.m_tot_sec * core.t0);
        out.delta_t = delta_t;
        self.sample_polarizations(ws, core.t0, dt, nt, out, f64::NAN);
        Ok(())
    }

    /// `h_+` and `h_x` sampled every `delta_t` seconds from each of `starts`
    /// (in seconds relative to the peak of the waveform), `length` samples
    /// each, written into the matching entry of `outs` (whose `epoch` is
    /// set to the start).
    ///
    /// This is [`Self::polarizations_with`] on grids of the caller's choosing:
    /// the surrogate is evaluated (or taken from `ws`) once and the modes
    /// summed once, then the polarizations are interpolated onto each grid
    /// with the same splines. With a start on the grid of
    /// [`Self::polarizations_with`] the samples agree with it. Samples
    /// outside the surrogate (before its start, set by `f_min` as there, or
    /// after the end of the ringdown) are zero.
    ///
    /// For example, a detector at a delay `dt_d` from the geocentre with the
    /// peak at geocentre time `t_c` sees, on its samples `t_0 + j delta_t`,
    /// the polarizations from start `t_0 - t_c - dt_d`.
    ///
    /// Returns the span `(start, end)` of the surrogate in seconds relative
    /// to the peak. The surrogate starts abruptly, so a caller whose grid
    /// moves relative to it (e.g. as the arrival time changes) should taper
    /// the start in continuous time, measured from `start`: otherwise the
    /// first non-zero sample, and anything computed from the samples, jumps
    /// whenever the start crosses a sample.
    ///
    /// # Panics
    ///
    /// If `starts` and `outs` have different lengths.
    pub fn polarizations_on_grids_with(
        &self,
        ws: &mut NRSurWorkspace,
        f_min: f64,
        delta_t: f64,
        starts: &[f64],
        length: usize,
        outs: &mut [TimeSeriesPolarizations],
    ) -> Result<(f64, f64), NRSurError> {
        assert_eq!(starts.len(), outs.len(), "one output per start time");
        let core = self.core(ws, f_min, delta_t)?;
        self.model_polarizations(ws, &core);
        let m_tot_sec = core.labels.m_tot_sec;
        let dt = delta_t / m_tot_sec;
        // The first time of the output of `polarizations_with`: earlier
        // samples are before the start of the surrogate for this `f_min`.
        let t_first = core.t0;
        for (&start, out) in starts.iter().zip(outs.iter_mut()) {
            out.epoch = start;
            out.delta_t = delta_t;
            self.sample_polarizations(ws, start / m_tot_sec, dt, length, out, 0.0);
            // Zero the samples before the start of the surrogate.
            // (with a tolerance, so that a grid through `t_first` keeps that
            // sample despite rounding)
            let before = ((t_first - start / m_tot_sec) / dt - 1e-9).ceil();
            if before > 0.0 {
                let before = (before as usize).min(length);
                out.hplus[..before].fill(0.0);
                out.hcross[..before].fill(0.0);
            }
        }
        let t_end = self.data.t_coorb.x[self.data.t_coorb.len() - 1];
        Ok((m_tot_sec * t_first, m_tot_sec * t_end))
    }

    /// Sum the modes into `h_+` and `h_x` on the coorbital times
    /// (`ws.hp_model`, `ws.hc_model`) and prepare their splines
    /// (`ws.c_plus`, `ws.c_cross`).
    fn model_polarizations(&self, ws: &mut NRSurWorkspace, core: &Core) {
        let n = ws.n_coorb;
        let l = &core.labels;
        ws.hp_model.fill(0.0);
        ws.hc_model.fill(0.0);
        let phi = FRAC_PI_2 - self.binary.phi_ref;
        for ell in 2..=LMAX {
            if !self.ell_active(ell) {
                continue;
            }
            for m in -(ell as i32)..=ell as i32 {
                let y = spin_weighted_spherical_harmonic(self.binary.inclination, phi, ell, m);
                let (c_re, c_im) = (y.re, y.im);
                let prefactor = if m % 2 != 0 && l.switched { -1.0 } else { 1.0 };
                let i = mode_index(ell, m);
                let (re, im) = (&ws.h_re[i * n..(i + 1) * n], &ws.h_im[i * n..(i + 1) * n]);
                for (((hp, hc), &hre), &him) in ws
                    .hp_model
                    .iter_mut()
                    .zip(ws.hc_model.iter_mut())
                    .zip(re)
                    .zip(im)
                {
                    *hp = fma(prefactor, fms(c_re, hre, c_im * him), *hp);
                    *hc = fnms(prefactor, fma(c_im, hre, c_re * him), *hc);
                }
            }
        }
        let distance = self.binary.luminosity_distance * 1e6 * PC_SI;
        let a0 = l.m_tot_sec * C_SI / distance;
        for h in ws.hp_model.iter_mut().chain(ws.hc_model.iter_mut()) {
            *h *= a0;
        }
        let knots = &self.data.t_coorb;
        if self.continuous {
            // Taper the end of the surrogate (where the ringdown is still
            // ~5e-4 of the peak) to zero, so that samples cross it smoothly.
            let t_end = knots.x[knots.len() - 1];
            for (k, &t) in knots.x.iter().enumerate().rev() {
                let x = (t_end - t) / END_TAPER;
                if x >= 1.0 {
                    break;
                }
                let w = (0.5 * PI * x).sin().powi(2);
                ws.hp_model[k] *= w;
                ws.hc_model[k] *= w;
            }
        }
        knots.second_derivatives(&ws.hp_model, &mut ws.c_plus, &mut ws.spline_z);
        knots.second_derivatives(&ws.hc_model, &mut ws.c_cross, &mut ws.spline_z);
    }

    /// Interpolate the prepared `h_+` and `h_x` onto `t0 + dt j`, `j < nt`
    /// (in units of the total mass), with `outside` beyond the coorbital
    /// times (LAL's NaN, or zero).
    fn sample_polarizations(
        &self,
        ws: &mut NRSurWorkspace,
        t0: f64,
        dt: f64,
        nt: usize,
        out: &mut TimeSeriesPolarizations,
        outside: f64,
    ) {
        let knots = &self.data.t_coorb;
        ws.grid.locate(knots, t0, dt, nt);
        for (y, c, h) in [
            (&ws.hp_model, &ws.c_plus, &mut out.hplus),
            (&ws.hc_model, &ws.c_cross, &mut out.hcross),
        ] {
            h.clear();
            h.resize(nt, 0.0);
            ws.grid
                .evaluate(knots, y, c, &mut ws.coef_b, &mut ws.coef_d, |j, v| {
                    h[j] = if v.is_nan() { outside } else { v }
                });
        }
    }

    /// The inertial-frame modes `h_{ell m}`, `ell <= 4`
    /// (`XLALSimInspiralPrecessingNRSurModes`).
    pub fn modes(&self, delta_t: f64, f_min: f64) -> Result<TimeSeriesModes, NRSurError> {
        let mut out = TimeSeriesModes::default();
        self.modes_with(&mut NRSurWorkspace::new(), delta_t, f_min, &mut out)?;
        Ok(out)
    }

    /// [`Self::modes`] reusing the buffers of `ws` and `out`.
    pub fn modes_with(
        &self,
        ws: &mut NRSurWorkspace,
        delta_t: f64,
        f_min: f64,
        out: &mut TimeSeriesModes,
    ) -> Result<(), NRSurError> {
        let core = self.core(ws, f_min, delta_t)?;
        let n = ws.n_coorb;
        let l = &core.labels;
        let distance = self.binary.luminosity_distance * 1e6 * PC_SI;
        let a0 = l.m_tot_sec * C_SI / distance;
        let (dt, nt) = self.output_grid(&core, delta_t);
        out.epoch = self.epoch(l.m_tot_sec * core.t0);
        out.delta_t = delta_t;
        let knots = &self.data.t_coorb;
        ws.grid.locate(knots, core.t0, dt, nt);
        let mut count = 0;
        for ell in 2..=LMAX {
            if !self.ell_active(ell) {
                continue;
            }
            for m in -(ell as i32)..=ell as i32 {
                let prefactor = if m % 2 != 0 && l.switched { -1.0 } else { 1.0 };
                let scale = prefactor * a0;
                let i = mode_index(ell, m);
                // scaled copies (gsl_vector_scale), keeping the modes cached
                let (re, im) = (&mut ws.hp_model, &mut ws.hc_model);
                for (out, h) in re.iter_mut().zip(&ws.h_re[i * n..(i + 1) * n]) {
                    *out = h * scale;
                }
                for (out, h) in im.iter_mut().zip(&ws.h_im[i * n..(i + 1) * n]) {
                    *out = h * scale;
                }
                knots.second_derivatives(re, &mut ws.c_plus, &mut ws.spline_z);
                knots.second_derivatives(im, &mut ws.c_cross, &mut ws.spline_z);
                if out.modes.len() <= count {
                    out.modes.push(Mode::default());
                }
                let mode = &mut out.modes[count];
                mode.ell = ell as u32;
                mode.m = m;
                let data = &mut mode.data;
                data.clear();
                data.resize(nt, Complex::new(0.0, 0.0));
                let (b, d) = (&mut ws.coef_b, &mut ws.coef_d);
                ws.grid
                    .evaluate(knots, re, &ws.c_plus, b, d, |j, v| data[j].re = v);
                ws.grid
                    .evaluate(knots, im, &ws.c_cross, b, d, |j, v| data[j].im = v);
                count += 1;
            }
        }
        out.modes.truncate(count);
        Ok(())
    }

    /// The dynamics of this binary on the dynamics nodes, as
    /// `XLALPrecessingNRSurDynamics` with LAL's waveform conventions
    /// (identity quaternion and zero orbital phase at the reference
    /// frequency `f_ref`, or `f_min` if `f_ref` is zero). The spins are
    /// those after the label switch (heavier black hole first).
    pub fn dynamics(&self, f_min: f64) -> Result<NRSurDynamics, NRSurError> {
        let l = self.labelled();
        let f_ref = if self.binary.f_ref == 0.0 {
            f_min
        } else {
            self.binary.f_ref
        };
        let omega_ref = (l.m_tot_sec * f_ref) * PI;
        self.data.dynamics(
            l.q,
            l.chi_a,
            l.chi_b,
            omega_ref,
            [1.0, 0.0, 0.0, 0.0],
            0.0,
            self.unlimited_extrapolation,
        )
    }
}

impl NRSur7dq4Data {
    /// `XLALPrecessingNRSurDynamics`: the coprecessing-frame quaternion,
    /// orbital phase and spins at the dynamics nodes for mass ratio `q >= 1`,
    /// spins `chi_a0`, `chi_b0` in the coorbital frame at the dimensionless
    /// orbital frequency `omega_ref` (rad / M), where the coprecessing
    /// quaternion is `init_quat` and the orbital phase `init_orbphase`.
    #[allow(clippy::too_many_arguments)]
    pub fn dynamics(
        &self,
        q: f64,
        chi_a0: [f64; 3],
        chi_b0: [f64; 3],
        omega_ref: f64,
        init_quat: [f64; 4],
        init_orbphase: f64,
        unlimited_extrapolation: bool,
    ) -> Result<NRSurDynamics, NRSurError> {
        // coorbital -> coprecessing frame
        let (s, c) = (init_orbphase.sin(), init_orbphase.cos());
        let rotate = |v: [f64; 3]| [fms(v[0], c, v[1] * s), fma(v[1], c, v[0] * s), v[2]];
        let (chi_a, chi_b) = (rotate(chi_a0), rotate(chi_b0));
        let n = self.t_ds.len();
        let mut states = vec![[0.0; N_STATE]; n];
        let mut omega_cache = vec![f64::NAN; n];
        Dynamics::new(self, q).integrate(
            &chi_a,
            &chi_b,
            omega_ref,
            init_orbphase,
            &init_quat,
            unlimited_extrapolation,
            false,
            &mut omega_cache,
            &mut states,
        )?;
        let column = |j: usize| states.iter().map(|s| s[j]).collect::<Vec<f64>>();
        Ok(NRSurDynamics {
            t: self.t_ds.x.clone(),
            quaternion: [column(0), column(1), column(2), column(3)],
            orbital_phase: column(4),
            chi_a: [column(5), column(6), column(7)],
            chi_b: [column(8), column(9), column(10)],
        })
    }
}

/// The ODE state at the reference point.
fn reference_state(chi_a: &[f64; 3], chi_b: &[f64; 3], quat: &[f64; 4], orbphase: f64) -> State {
    let mut y = [0.0; N_STATE];
    y[..4].copy_from_slice(quat);
    y[4] = orbphase;
    y[5..8].copy_from_slice(chi_a);
    y[8..11].copy_from_slice(chi_b);
    y
}

/// Mutable rows `a` and `b` (distinct) of a `[row][n]` buffer.
fn two_rows(buffer: &mut [f64], a: usize, b: usize, n: usize) -> (&mut [f64], &mut [f64]) {
    assert_ne!(a, b);
    if a < b {
        let (lo, hi) = buffer.split_at_mut(b * n);
        (&mut lo[a * n..(a + 1) * n], &mut hi[..n])
    } else {
        let (lo, hi) = buffer.split_at_mut(a * n);
        (&mut hi[..n], &mut lo[b * n..(b + 1) * n])
    }
}

/// `PrecessingNRSur_eval_data_piece`: the fits at the empirical nodes
/// (with the coorbital spins at the node times) and the empirical
/// interpolant `basis^T nodes` (`gsl_blas_dgemv`, accumulated column by
/// column, which reads the basis contiguously), written to `out`.
fn eval_data_piece(
    data: &DataPiece,
    params: &MassRatioParameters,
    q: f64,
    spins: &[f64],
    n: usize,
    nodes: &mut Vec<f64>,
    out: &mut [f64],
) {
    nodes.clear();
    for (fit, &t) in data.fits.iter().zip(&data.node_indices) {
        let x = [
            q,
            spins[t],
            spins[n + t],
            spins[2 * n + t],
            spins[3 * n + t],
            spins[4 * n + t],
            spins[5 * n + t],
        ];
        nodes.push(fit.eval(&params.table(&x)));
    }
    let out = &mut out[..n];
    let mut rows = nodes.iter().zip(data.basis.chunks_exact(n));
    // y = x_0 A_0, then y += x_j A_j column by column
    match rows.next() {
        Some((&a, row)) => {
            for (o, &b) in out.iter_mut().zip(row) {
                *o = a * b;
            }
        }
        None => out.fill(0.0),
    }
    for (&a, row) in rows {
        for (o, &b) in out.iter_mut().zip(row) {
            *o = fma(a, b, *o);
        }
    }
}

/// The output times `t0 + dt j`, `j < nt`, located in the coorbital
/// times as `gsl_interp_accel_find` locates them.
#[derive(Clone, Debug, Default)]
struct OutputGrid {
    /// Interval of each time ([`usize::MAX`] outside the knots, where GSL
    /// returns NaN with a domain error).
    index: Vec<usize>,
    /// Offset from the left knot of the interval.
    delx: Vec<f64>,
    /// The first and last intervals used.
    first: usize,
    last: usize,
}

impl OutputGrid {
    fn locate(&mut self, knots: &FixedKnots, t0: f64, dt: f64, nt: usize) {
        self.index.clear();
        self.delx.clear();
        self.first = usize::MAX;
        self.last = 0;
        let mut hint = 0;
        for j in 0..nt {
            let t = fma(dt, j as f64, t0);
            if knots.contains(t) {
                hint = knots.locate_from(t, hint);
                self.index.push(hint);
                self.delx.push(t - knots.x[hint]);
                self.first = self.first.min(hint);
                self.last = self.last.max(hint);
            } else {
                self.index.push(usize::MAX);
                self.delx.push(f64::NAN);
            }
        }
    }

    /// Evaluate the natural spline `(y, c)` at the grid (`gsl_spline_eval`),
    /// passing each value to `set`. `b` and `d` hold the coefficients of the
    /// intervals.
    fn evaluate(
        &self,
        knots: &FixedKnots,
        y: &[f64],
        c: &[f64],
        b: &mut [f64],
        d: &mut [f64],
        mut set: impl FnMut(usize, f64),
    ) {
        if self.first <= self.last {
            knots.intervals(self.first..self.last + 1, y, c, b, d);
        }
        for (j, (&i, &delx)) in self.index.iter().zip(&self.delx).enumerate() {
            set(
                j,
                if i == usize::MAX {
                    f64::NAN
                } else {
                    cubic(y[i], b[i], c[i], d[i], delx)
                },
            );
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gps_epoch_rounds_like_lal() {
        assert_eq!(gps_epoch(-0.9586624661), -0.958662466);
        assert_eq!(gps_epoch(1.25), 1.25);
        assert_eq!(gps_epoch(-2.0000000004), -2.0);
    }

    #[test]
    fn harmonics_satisfy_the_addition_theorem() {
        // sum_m |_{-2}Y_{lm}|^2 = (2l + 1) / (4 pi)
        for theta in [0.0, 0.3, 1.2, 2.9] {
            for ell in 2..=4usize {
                let sum: f64 = (-(ell as i32)..=ell as i32)
                    .map(|m| spin_weighted_spherical_harmonic(theta, 0.7, ell, m).norm_sqr())
                    .sum();
                let expected = (2 * ell + 1) as f64 / (4.0 * PI);
                assert!((sum - expected).abs() < 1e-14, "{ell} {theta}");
            }
        }
    }

    /// Against `XLALSpinWeightedSphericalHarmonic`, to a few ulp.
    #[test]
    fn harmonics_match_lal() {
        #[derive(serde::Deserialize)]
        struct Case {
            theta: Vec<f64>,
            phi: Vec<f64>,
            values: Vec<(usize, i32, f64, f64)>,
        }
        #[derive(serde::Deserialize)]
        struct File {
            cases: Vec<Case>,
        }
        let path = crate::reference::data_path("spin_weighted_harmonics");
        let file: File = serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap();
        let case = &file.cases[0];
        let mut values = case.values.iter();
        let mut mismatches = Vec::new();
        for (&theta, &phi) in case.theta.iter().zip(&case.phi) {
            for ell in 2..=4usize {
                for m in -(ell as i32)..=ell as i32 {
                    let &(l, mm, re, im) = values.next().unwrap();
                    assert_eq!((l, mm), (ell, m));
                    let y = spin_weighted_spherical_harmonic(theta, phi, ell, m);
                    // |Y| <= 1.3, so to a few ulp
                    let matches = (y.re - re).abs().max((y.im - im).abs()) <= 4.0 * f64::EPSILON;
                    if !matches {
                        mismatches.push(format!("({ell}, {m}) theta={theta}"));
                    }
                }
            }
        }
        assert!(mismatches.is_empty(), "{}", mismatches.join("\n"));
    }
}

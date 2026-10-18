//! Adaptive Runge-Kutta-Fehlberg integration with Hermite-like resampling.
//!
//! This is a port of `XLALAdaptiveRungeKutta4Hermite` (lal 7.7,
//! `LALAdaptiveRungeKuttaIntegrator.c`) together with the pieces of GSL 2.8's
//! (legacy) `gsl_odeiv` interface it uses:
//!
//! - `gsl_odeiv_step_rkf45` (`ode-initval/rkf45.c`),
//! - `gsl_odeiv_control_y_new` (`ode-initval/cstd.c`, standard control with
//!   `a_y = 1`, `a_dydt = 0`),
//! - `gsl_odeiv_evolve_apply` (`ode-initval/evolve.c`).
//!
//! The operation order (including the fused multiply-adds emitted by clang,
//! see [`super::fp`]) follows the C sources so that the sequence of accepted
//! steps, and therefore the output, agrees with LAL to rounding precision.
//!
//! One deliberate (bit-for-bit neutral) optimisation: `gsl_odeiv_evolve_apply`
//! evaluates the derivative at the start of every step, which for this stepper
//! always equals the derivative evaluated at the end of the previous step.  We
//! reuse the latter instead of recomputing it.

use super::fp::{fma, fnms};

/// GSL success status.
pub const GSL_SUCCESS: i32 = 0;

/// A system of first order ODEs, `dy/dt = f(t, y)`, with LAL's stopping test.
///
/// Status codes follow the GSL/LAL convention: `0` for success, anything
/// else signals a failure (for [`OdeSystem::derivatives`]) or the reason to
/// stop (for [`OdeSystem::stop`]).
pub trait OdeSystem<const N: usize> {
    /// Evaluate the derivatives; returns `0` on success.
    fn derivatives(&self, t: f64, y: &[f64; N], dydt: &mut [f64; N]) -> i32;

    /// Stopping test called after every successful step with the state and
    /// its derivative; returns `0` to continue.
    fn stop(&mut self, t: f64, y: &[f64; N], dydt: &[f64; N]) -> i32;
}

// Runge-Kutta-Fehlberg coefficients from GSL's rkf45.c.
const AH: [f64; 5] = [1.0 / 4.0, 3.0 / 8.0, 12.0 / 13.0, 1.0, 1.0 / 2.0];
const B3: [f64; 2] = [3.0 / 32.0, 9.0 / 32.0];
const B4: [f64; 3] = [1932.0 / 2197.0, -7200.0 / 2197.0, 7296.0 / 2197.0];
const B5: [f64; 4] = [
    8341.0 / 4104.0,
    -32832.0 / 4104.0,
    29440.0 / 4104.0,
    -845.0 / 4104.0,
];
const B6: [f64; 5] = [
    -6080.0 / 20520.0,
    41040.0 / 20520.0,
    -28352.0 / 20520.0,
    9295.0 / 20520.0,
    -5643.0 / 20520.0,
];
const C1: f64 = 902880.0 / 7618050.0;
const C3: f64 = 3953664.0 / 7618050.0;
const C4: f64 = 3855735.0 / 7618050.0;
const C5: f64 = -1371249.0 / 7618050.0;
const C6: f64 = 277020.0 / 7618050.0;
const EC: [f64; 7] = [
    0.0,
    1.0 / 360.0,
    0.0,
    -128.0 / 4275.0,
    -2197.0 / 75240.0,
    1.0 / 50.0,
    2.0 / 55.0,
];
/// `rkf45_order`
const RKF45_ORDER: u32 = 5;

/// State of the `rkf45` stepper (`rkf45_state_t`).
#[derive(Clone, Debug)]
pub struct Rkf45<const N: usize> {
    pub k1: [f64; N],
    pub k2: [f64; N],
    pub k3: [f64; N],
    pub k4: [f64; N],
    pub k5: [f64; N],
    pub k6: [f64; N],
    pub y0: [f64; N],
    ytmp: [f64; N],
}

impl<const N: usize> Default for Rkf45<N> {
    fn default() -> Self {
        Self {
            k1: [0.0; N],
            k2: [0.0; N],
            k3: [0.0; N],
            k4: [0.0; N],
            k5: [0.0; N],
            k6: [0.0; N],
            y0: [0.0; N],
            ytmp: [0.0; N],
        }
    }
}

impl<const N: usize> Rkf45<N> {
    /// Port of `rkf45_apply` (always called with `dydt_in` and `dydt_out`).
    ///
    /// On failure `y` is left unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn apply<S: OdeSystem<N>>(
        &mut self,
        t: f64,
        h: f64,
        y: &mut [f64; N],
        yerr: &mut [f64; N],
        dydt_in: &[f64; N],
        dydt_out: &mut [f64; N],
        sys: &S,
    ) -> i32 {
        self.y0 = *y;
        self.k1 = *dydt_in;
        let k1 = &self.k1;

        // k2 step
        let ah0h = AH[0] * h;
        for i in 0..N {
            self.ytmp[i] = fma(ah0h, k1[i], y[i]);
        }
        let s = sys.derivatives(fma(AH[0], h, t), &self.ytmp, &mut self.k2);
        if s != GSL_SUCCESS {
            return s;
        }
        let k2 = &self.k2;

        // k3 step
        for i in 0..N {
            self.ytmp[i] = fma(h, fma(B3[0], k1[i], B3[1] * k2[i]), y[i]);
        }
        let s = sys.derivatives(fma(AH[1], h, t), &self.ytmp, &mut self.k3);
        if s != GSL_SUCCESS {
            return s;
        }
        let k3 = &self.k3;

        // k4 step
        for i in 0..N {
            let d = fma(B4[2], k3[i], fma(B4[0], k1[i], B4[1] * k2[i]));
            self.ytmp[i] = fma(h, d, y[i]);
        }
        let s = sys.derivatives(fma(AH[2], h, t), &self.ytmp, &mut self.k4);
        if s != GSL_SUCCESS {
            return s;
        }
        let k4 = &self.k4;

        // k5 step
        for i in 0..N {
            let mut d = fma(B5[0], k1[i], B5[1] * k2[i]);
            d = fma(B5[2], k3[i], d);
            d = fma(B5[3], k4[i], d);
            self.ytmp[i] = fma(h, d, y[i]);
        }
        let s = sys.derivatives(fma(AH[3], h, t), &self.ytmp, &mut self.k5);
        if s != GSL_SUCCESS {
            return s;
        }
        let k5 = &self.k5;

        // k6 step
        for i in 0..N {
            let mut d = fma(B6[0], k1[i], B6[1] * k2[i]);
            d = fma(B6[2], k3[i], d);
            d = fma(B6[3], k4[i], d);
            d = fma(B6[4], k5[i], d);
            self.ytmp[i] = fma(h, d, y[i]);
        }
        let s = sys.derivatives(fma(AH[4], h, t), &self.ytmp, &mut self.k6);
        if s != GSL_SUCCESS {
            return s;
        }
        let k6 = &self.k6;

        // final sum
        for i in 0..N {
            let mut d = fma(C1, k1[i], C3 * k3[i]);
            d = fma(C4, k4[i], d);
            d = fma(C5, k5[i], d);
            d = fma(C6, k6[i], d);
            y[i] = fma(h, d, y[i]);
        }

        // derivatives at output
        let s = sys.derivatives(t + h, y, dydt_out);
        if s != GSL_SUCCESS {
            *y = self.y0;
            return s;
        }

        // difference between 4th and 5th order
        for i in 0..N {
            let mut d = fma(EC[1], k1[i], EC[3] * k3[i]);
            d = fma(EC[4], k4[i], d);
            d = fma(EC[5], k5[i], d);
            d = fma(EC[6], k6[i], d);
            yerr[i] = h * d;
        }
        GSL_SUCCESS
    }
}

/// Result of [`StandardControl::hadjust`] (`GSL_ODEIV_HADJ_*`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StepAdjustment {
    Increased,
    Unchanged,
    Decreased,
}

/// The standard GSL step size control, `gsl_odeiv_control_standard`.
#[derive(Clone, Copy, Debug)]
pub struct StandardControl {
    pub eps_abs: f64,
    pub eps_rel: f64,
    pub a_y: f64,
    pub a_dydt: f64,
}

impl StandardControl {
    /// `gsl_odeiv_control_y_new`
    pub fn y_new(eps_abs: f64, eps_rel: f64) -> Self {
        Self {
            eps_abs,
            eps_rel,
            a_y: 1.0,
            a_dydt: 0.0,
        }
    }

    /// Port of `std_control_hadjust`.
    #[allow(clippy::manual_clamp)] // mirrors the C branches
    pub fn hadjust<const N: usize>(
        &self,
        order: u32,
        y: &[f64; N],
        yerr: &[f64; N],
        yp: &[f64; N],
        h: &mut f64,
    ) -> StepAdjustment {
        const S: f64 = 0.9;
        let h_old = *h;
        let mut rmax = f64::MIN_POSITIVE; // DBL_MIN
        for i in 0..N {
            let d0 = fma(
                self.eps_rel,
                fma(self.a_y, y[i].abs(), self.a_dydt * (h_old * yp[i]).abs()),
                self.eps_abs,
            );
            let r = yerr[i].abs() / d0.abs();
            // GSL_MAX_DBL: NaN never replaces the running maximum
            if r > rmax {
                rmax = r;
            }
        }
        if rmax > 1.1 {
            // decrease step, no more than factor of 5, but a fraction S more
            // than scaling suggests (for better accuracy)
            let mut r = S / rmax.powf(1.0 / order as f64);
            if r < 0.2 {
                r = 0.2;
            }
            *h = r * h_old;
            StepAdjustment::Decreased
        } else if rmax < 0.5 {
            // increase step, no more than factor of 5
            let mut r = S / rmax.powf(1.0 / (order as f64 + 1.0));
            if r > 5.0 {
                r = 5.0;
            }
            if r < 1.0 {
                // don't allow any decrease caused by S<1
                r = 1.0;
            }
            *h = r * h_old;
            StepAdjustment::Increased
        } else {
            StepAdjustment::Unchanged
        }
    }
}

/// `gsl_odeiv_evolve` with an `rkf45` stepper and standard control.
#[derive(Clone, Debug)]
pub struct Evolve<const N: usize> {
    pub step: Rkf45<N>,
    pub control: StandardControl,
    y0: [f64; N],
    yerr: [f64; N],
    dydt_in: [f64; N],
    pub dydt_out: [f64; N],
    /// Whether `dydt_in` holds the derivative at the current state.
    dydt_in_valid: bool,
    pub count: u64,
    pub failed_steps: u64,
    pub last_step: f64,
}

impl<const N: usize> Evolve<N> {
    pub fn new(control: StandardControl) -> Self {
        Self {
            step: Rkf45::default(),
            control,
            y0: [0.0; N],
            yerr: [0.0; N],
            dydt_in: [0.0; N],
            dydt_out: [0.0; N],
            dydt_in_valid: false,
            count: 0,
            failed_steps: 0,
            last_step: 0.0,
        }
    }

    /// `gsl_odeiv_evolve_reset` (plus `gsl_odeiv_step_reset`).
    pub fn reset(&mut self) {
        self.step = Rkf45::default();
        self.dydt_in_valid = false;
        self.count = 0;
        self.failed_steps = 0;
        self.last_step = 0.0;
    }

    /// Port of `gsl_odeiv_evolve_apply`.
    pub fn apply<S: OdeSystem<N>>(
        &mut self,
        sys: &S,
        t: &mut f64,
        t1: f64,
        h: &mut f64,
        y: &mut [f64; N],
    ) -> i32 {
        let t0 = *t;
        let mut h0 = *h;
        let dt = t1 - t0;

        if (dt < 0.0 && h0 > 0.0) || (dt > 0.0 && h0 < 0.0) {
            // "step direction must match interval direction"
            return GSL_EINVAL;
        }

        self.y0 = *y;

        // Calculate initial dydt once (rkf45 can use dydt_in).  The
        // derivative at the end of the previous successful step is the
        // derivative at the current state, so reuse it.
        if !self.dydt_in_valid {
            let status = sys.derivatives(t0, y, &mut self.dydt_in);
            if status != GSL_SUCCESS {
                return status;
            }
            self.dydt_in_valid = true;
        }

        loop {
            let final_step = if (dt >= 0.0 && h0 > dt) || (dt < 0.0 && h0 < dt) {
                h0 = dt;
                true
            } else {
                false
            };

            let step_status = self.step.apply(
                t0,
                h0,
                y,
                &mut self.yerr,
                &self.dydt_in,
                &mut self.dydt_out,
                sys,
            );

            if step_status != GSL_SUCCESS {
                *h = h0; // notify user of step-size which caused the failure
                *t = t0; // restore original t value
                return step_status;
            }

            self.count += 1;
            self.last_step = h0;
            *t = if final_step { t1 } else { t0 + h0 };

            // Check error and attempt to adjust the step.
            let h_old = h0;
            let adjustment =
                self.control
                    .hadjust(RKF45_ORDER, y, &self.yerr, &self.dydt_out, &mut h0);
            if adjustment == StepAdjustment::Decreased {
                // Check that the reported status is correct (i.e. an actual
                // decrease in h0 occured) and the suggested h0 will change
                // the time by at least 1 ulp
                let t_curr = *t;
                let t_next = *t + h0;
                if h0.abs() < h_old.abs() && t_next != t_curr {
                    // Step was decreased. Undo step, and try again with new h0.
                    *y = self.y0;
                    self.failed_steps += 1;
                    continue;
                } else {
                    h0 = h_old; // keep current step size
                }
            }
            break;
        }

        // The state has moved on: the next initial derivative is dydt_out.
        self.dydt_in = self.dydt_out;
        *h = h0; // suggest step size for next time-step
        GSL_SUCCESS
    }
}

/// GSL's `GSL_EINVAL`.
pub const GSL_EINVAL: i32 = 4;

/// Errors from [`adaptive_rk4_hermite`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum IntegratorError {
    /// `(tend_in - tinit)` and `deltat` have different signs.
    InvalidInterval,
}

/// Evenly sampled output of [`adaptive_rk4_hermite`].
#[derive(Clone, Debug)]
pub struct HermiteOutput<const N: usize> {
    /// Sample times.
    pub t: Vec<f64>,
    /// State at each sample time.
    pub y: Vec<[f64; N]>,
    /// The integrator's return code: `0`, the failing status of the
    /// derivative function after exhausting the retries, or the value
    /// returned by the stopping test.
    pub return_code: i32,
}

/// Settings of the `LALAdaptiveRungeKuttaIntegrator` struct.
#[derive(Clone, Copy, Debug)]
pub struct IntegratorSettings {
    pub eps_abs: f64,
    pub eps_rel: f64,
    /// Number of retries with a ten times smaller step after a failure of
    /// the derivative function (LAL default: 6).
    pub retries: u32,
    /// Stop only on the stopping test, `tend` is only used to size buffers.
    pub stop_on_test_only: bool,
}

/// Where [`adaptive_rk4_hermite_into`] records the evenly spaced samples.
pub(crate) trait HermiteSink<const N: usize> {
    /// Called once before the first sample with LAL's estimate of the number
    /// of samples (which can be exceeded).
    fn reserve(&mut self, n: usize);
    /// The state `y` at time `t`.
    fn push(&mut self, t: f64, y: &[f64; N]);
}

impl<const N: usize> HermiteSink<N> for HermiteOutput<N> {
    fn reserve(&mut self, n: usize) {
        self.t.reserve(n);
        self.y.reserve(n);
    }

    #[inline]
    fn push(&mut self, t: f64, y: &[f64; N]) {
        self.t.push(t);
        self.y.push(*y);
    }
}

/// Port of `XLALAdaptiveRungeKutta4Hermite`.
///
/// Integrates from `tinit` with adaptive RKF45 steps and records the state
/// every `deltat` using the fourth-order interpolant of the RKF45 step.  The
/// first sample is `yinit`; on return `yinit` holds the last *interpolated*
/// sample, as in LAL.
pub fn adaptive_rk4_hermite<S: OdeSystem<N>, const N: usize>(
    sys: &mut S,
    settings: &IntegratorSettings,
    yinit: &mut [f64; N],
    tinit: f64,
    tend_in: f64,
    deltat: f64,
) -> Result<HermiteOutput<N>, IntegratorError> {
    let mut output = HermiteOutput {
        t: Vec::new(),
        y: Vec::new(),
        return_code: 0,
    };
    output.return_code =
        adaptive_rk4_hermite_into(sys, settings, yinit, tinit, tend_in, deltat, &mut output)?;
    Ok(output)
}

/// [`adaptive_rk4_hermite`] recording the samples in `sink`; returns the
/// integrator's return code.
pub(crate) fn adaptive_rk4_hermite_into<S: OdeSystem<N>, K: HermiteSink<N>, const N: usize>(
    sys: &mut S,
    settings: &IntegratorSettings,
    yinit: &mut [f64; N],
    tinit: f64,
    tend_in: f64,
    deltat: f64,
    sink: &mut K,
) -> Result<i32, IntegratorError> {
    let mut tend = tend_in;
    if settings.stop_on_test_only {
        tend = if tend < tinit {
            f64::NEG_INFINITY
        } else {
            f64::INFINITY
        };
    }

    // outputlen = ((int)(tend_in - tinit) / deltat): note the cast binds to
    // the numerator only.
    let length_estimate = ((tend_in - tinit) as i32 as f64 / deltat) as i64;
    if length_estimate < 0 {
        return Err(IntegratorError::InvalidInterval);
    }
    // LAL's estimate is often slightly short (its callers estimate the
    // duration at Newtonian order); a margin avoids growing the buffers.
    let capacity = (length_estimate + 2).clamp(2, 1 << 24) as usize;
    sink.reserve(capacity + capacity / 8);

    let mut evolve = Evolve::<N>::new(StandardControl::y_new(settings.eps_abs, settings.eps_rel));
    let mut return_code = 0;
    let mut retries = settings.retries;
    let mut t = tinit;
    let mut tintp = tinit;
    let mut h = deltat;

    // Copy over first step.
    sink.push(tinit, yinit);
    let mut last = *yinit;

    evolve.reset();

    // Enter evolution loop.  NOTE: we *always* take at least one step.
    loop {
        let told = t;
        let status = evolve.apply(sys, &mut t, tend, &mut h, yinit);

        // Check for failure, retry if haven't retried too many times already.
        if status != GSL_SUCCESS {
            if retries > 0 {
                retries -= 1;
                // Retries to spare; reduce h, try again.
                h /= 10.0;
                continue;
            } else {
                // Out of retries, bail with status code.
                return_code = status;
                break;
            }
        } else {
            // Successful step, reset retry counter.
            retries = settings.retries;
        }

        // Now interpolate until we would go past the current integrator
        // time, t.  Note we square to get an absolute value, because we may
        // be integrating t in the positive or negative direction.
        let rkf = &evolve.step;
        while (tintp + deltat) * (tintp + deltat) < t * t {
            tintp += deltat;
            // tintp = told + (t-told)*theta, 0 <= theta <= 1.  We have to
            // compute h = (t-told) because the integrator returns a
            // suggested next h, not the actual stepsize taken.
            let h_used = t - told;
            let theta = (tintp - told) / h_used;

            // These are the interpolating coefficients for y(t + h*theta) =
            // ynew + i1*h*k1 + i5*h*k5 + i6*h*k6 + O(h^4).
            let i0 = fma(theta * theta, fnms(4.0, theta, 3.0), 1.0);
            let i1 = -theta * (theta - 1.0);
            let i6 = -4.0 * theta * theta * (theta - 1.0);
            let iend = theta * theta * fma(4.0, theta, -3.0);

            let h_i1 = h_used * i1;
            let h_i6 = h_used * i6;
            let mut ytemp = [0.0; N];
            for i in 0..N {
                let mut x = fma(i0, rkf.y0[i], iend * yinit[i]);
                x = fma(h_i1, rkf.k1[i], x);
                ytemp[i] = fma(h_i6, rkf.k6[i], x);
            }
            sink.push(tintp, &ytemp);
            last = ytemp;
        }

        // Now that we have recorded the last interpolated step that we
        // could, check for termination criteria.
        if !settings.stop_on_test_only && t >= tend {
            break;
        }

        // Call the stopping function with the last value of y and dydt from
        // the integrator.
        let status = sys.stop(t, yinit, &evolve.dydt_out);
        if status != GSL_SUCCESS {
            return_code = status;
            break;
        }
    }

    // Store the final *interpolated* sample in yinit.
    *yinit = last;

    Ok(return_code)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Simple harmonic oscillator, y'' = -y.
    struct Oscillator {
        stop_time: f64,
    }

    impl OdeSystem<2> for Oscillator {
        fn derivatives(&self, _t: f64, y: &[f64; 2], dydt: &mut [f64; 2]) -> i32 {
            dydt[0] = y[1];
            dydt[1] = -y[0];
            GSL_SUCCESS
        }

        fn stop(&mut self, t: f64, _y: &[f64; 2], _dydt: &[f64; 2]) -> i32 {
            if t.abs() > self.stop_time {
                1
            } else {
                0
            }
        }
    }

    #[test]
    fn rkf45_step_is_accurate() {
        let sys = Oscillator { stop_time: 0.0 };
        let mut step = Rkf45::<2>::default();
        let mut y = [1.0, 0.0];
        let mut yerr = [0.0; 2];
        let mut dydt_in = [0.0; 2];
        let mut dydt_out = [0.0; 2];
        sys.derivatives(0.0, &y, &mut dydt_in);
        let h = 0.01;
        assert_eq!(
            step.apply(0.0, h, &mut y, &mut yerr, &dydt_in, &mut dydt_out, &sys),
            GSL_SUCCESS
        );
        // fifth order accurate step: local error O(h^6)
        assert!((y[0] - h.cos()).abs() < 1e-13);
        assert!((y[1] + h.sin()).abs() < 1e-13);
        assert!(yerr.iter().all(|e| e.abs() < 1e-10));
        assert_eq!(dydt_out, [y[1], -y[0]]);
    }

    #[test]
    fn hadjust_matches_gsl_rules() {
        let control = StandardControl::y_new(1e-6, 0.0);
        let y = [1.0];
        let yp = [0.0];
        // large error: decrease, at most by a factor 5
        let mut h = 1.0;
        let adj = control.hadjust(5, &y, &[1.0], &yp, &mut h);
        assert_eq!(adj, StepAdjustment::Decreased);
        assert_eq!(h, 0.2);
        // tiny error: increase, at most by a factor 5
        let mut h = 1.0;
        let adj = control.hadjust(5, &y, &[1e-20], &yp, &mut h);
        assert_eq!(adj, StepAdjustment::Increased);
        assert_eq!(h, 5.0);
        // error within [0.5, 1.1] of the tolerance: unchanged
        let mut h = 1.0;
        let adj = control.hadjust(5, &y, &[1e-6], &yp, &mut h);
        assert_eq!(adj, StepAdjustment::Unchanged);
        assert_eq!(h, 1.0);
        // NaN errors are ignored, as GSL_MAX_DBL never selects them
        let mut h = 1.0;
        let adj = control.hadjust(5, &[1.0, 1.0], &[f64::NAN, 1e-20], &[0.0; 2], &mut h);
        assert_eq!(adj, StepAdjustment::Increased);
    }

    #[test]
    fn hermite_output_is_uniform_and_accurate() {
        for direction in [1.0, -1.0] {
            let mut sys = Oscillator { stop_time: 10.0 };
            let settings = IntegratorSettings {
                eps_abs: 1e-12,
                eps_rel: 1e-12,
                retries: 6,
                stop_on_test_only: true,
            };
            let mut y = [1.0, 0.0];
            let dt = direction * 0.05;
            let out = adaptive_rk4_hermite(&mut sys, &settings, &mut y, 0.0, 10.0 * direction, dt)
                .unwrap();
            assert_eq!(out.return_code, 1);
            assert_eq!(out.t.len(), out.y.len());
            let mut expected_t = 0.0;
            let mut max_error: f64 = 0.0;
            for (t, y) in out.t.iter().zip(&out.y) {
                assert_eq!(*t, expected_t);
                expected_t += dt;
                max_error = max_error
                    .max((y[0] - t.cos()).abs())
                    .max((y[1] + t.sin()).abs());
            }
            // the interpolant is fourth order accurate in the step size
            println!("max interpolation error {max_error:e}");
            assert!(max_error < 1e-8, "{max_error:e}");
            assert!(out.t.last().unwrap().abs() > 9.9);
            assert_eq!(y, *out.y.last().unwrap());
        }
    }

    #[test]
    fn mismatched_direction_is_an_error() {
        let mut sys = Oscillator { stop_time: 1.0 };
        let settings = IntegratorSettings {
            eps_abs: 1e-12,
            eps_rel: 1e-12,
            retries: 6,
            stop_on_test_only: true,
        };
        let mut y = [1.0, 0.0];
        assert!(adaptive_rk4_hermite(&mut sys, &settings, &mut y, 0.0, 10.0, -0.1).is_err());
    }
}

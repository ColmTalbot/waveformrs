//! The dynamics surrogate: integration of the ODE for the coprecessing-frame
//! quaternion, the orbital phase and the spins
//! (`PrecessingNRSur_IntegrateDynamics` and the functions it calls).
//!
//! The state is `y = [q0, qx, qy, qz, orbphase, chiAx, chiAy, chiAz, chiBx,
//! chiBy, chiBz]` with the spins in the coprecessing frame. The time
//! derivatives are fits evaluated at the dynamics nodes; the ODE is started
//! with three RK4 steps and integrated with a 4th-order Adams-Bashforth
//! scheme on the (non-uniform) nodes.

use crate::nrsur::data::{DynamicsNode, NRSur7dq4Data};
use crate::nrsur::fit::{MassRatioParameters, PowerTable};
use crate::nrsur::spline::polynomial4;
use crate::nrsur::trig::sincos;
use crate::nrsur::NRSurError;
use crate::spintaylor::fp::{fma, fms, fnms};

/// Length of the ODE state.
pub(crate) const N_STATE: usize = 11;
pub(crate) type State = [f64; N_STATE];

/// `omega_ref` above which `PrecessingNRSur_get_t_ref` refuses (the merger
/// frequency is larger for all masses and spins).
pub(crate) const OMEGA_REF_MAX: f64 = 0.201;

/// `PrecessingNRSur_ds_fit_x`: the fit parameters with the spins rotated
/// into the coorbital frame.
#[inline]
fn fit_x(q: f64, y: &State, sp: f64, cp: f64) -> [f64; 7] {
    [
        q,
        fma(y[5], cp, y[6] * sp),
        fma(-y[5], sp, y[6] * cp),
        y[7],
        fma(y[8], cp, y[9] * sp),
        fma(-y[8], sp, y[9] * cp),
        y[10],
    ]
}

/// `PrecessingNRSur_normalize_y`: keep the quaternion a unit quaternion and
/// the spin magnitudes constant.
#[inline]
pub(crate) fn normalize_y(norm_a: f64, norm_b: f64, y: &mut State) {
    let mut sum = 0.0;
    for v in &y[0..4] {
        sum = fma(*v, *v, sum);
    }
    let n_q = sum.sqrt();
    let mut sum = 0.0;
    for v in &y[5..8] {
        sum = fma(*v, *v, sum);
    }
    let n_a = sum.sqrt();
    let mut sum = 0.0;
    for v in &y[8..11] {
        sum = fma(*v, *v, sum);
    }
    let n_b = sum.sqrt();
    for v in &mut y[0..4] {
        *v /= n_q;
    }
    for v in &mut y[5..8] {
        *v = *v * norm_a / n_a;
    }
    for v in &mut y[8..11] {
        *v = *v * norm_b / n_b;
    }
}

/// `NRSur_ab4_dy`: the 4th-order Adams-Bashforth increment from the
/// derivatives `k1..k4` at four nodes separated by `t12, t23, t34` to the
/// next node at `t45`.
#[allow(clippy::too_many_arguments)]
#[inline]
fn ab4_dy(
    k1: &State,
    k2: &State,
    k3: &State,
    k4: &State,
    t12: f64,
    t23: f64,
    t34: f64,
    t45: f64,
) -> State {
    let t13 = t12 + t23;
    let t14 = t13 + t34;
    let t24 = t23 + t34;
    let denom1 = t12 * t13 * t14;
    let denom2 = t12 * t23 * t24;
    let denom3 = t13 * t23 * t34;
    let dx0 = t45;
    let dx1 = t45 * t45 / 2.0;
    let dx2 = t45 * t45 * t45 / 3.0;
    let dx3 = dx1 * dx1;
    let a1 = -1.0 * t24 * t34 / denom1;
    let a2 = t14 * t34 / denom2;
    let a3 = -1.0 * t14 * t24 / denom3;
    let a4 = -(a1 + a2 + a3);
    let b1 = -1.0 * (t24 + t34) / denom1;
    let b2 = (t14 + t34) / denom2;
    let b3 = -1.0 * (t14 + t24) / denom3;
    let b4 = -(b1 + b2 + b3);
    let c1 = -1.0 / denom1;
    let c2 = 1.0 / denom2;
    let c3 = -1.0 / denom3;
    let c4 = -(c1 + c2 + c3);
    let mut dy = [0.0; N_STATE];
    for i in 0..N_STATE {
        let mut tmp = k4[i] * dx0;
        let a = fma(a4, k4[i], fma(a3, k3[i], fma(a1, k1[i], a2 * k2[i])));
        tmp = fma(a, dx1, tmp);
        let b = fma(b4, k4[i], fma(b3, k3[i], fma(b1, k1[i], b2 * k2[i])));
        tmp = fma(b, dx2, tmp);
        let c = fma(c4, k4[i], fma(c3, k3[i], fma(c1, k1[i], c2 * k2[i])));
        tmp = fma(c, dx3, tmp);
        dy[i] = tmp;
    }
    dy
}

/// The RK4 combination `y + h/6 (k1 + 2 k2 + 2 k3 + k4)` (`h` may be
/// negative only through `sign`, as LAL writes the backward steps with a
/// minus sign).
#[inline]
fn rk4_combine(
    y: &State,
    h: f64,
    backwards: bool,
    k1: &State,
    k2: &State,
    k3: &State,
    k4: &State,
) -> State {
    let mut out = [0.0; N_STATE];
    let h6 = h / 6.0;
    for j in 0..N_STATE {
        let sum = fma(2.0, k3[j], fma(2.0, k2[j], k1[j])) + k4[j];
        out[j] = if backwards {
            fnms(h6, sum, y[j])
        } else {
            fma(h6, sum, y[j])
        };
    }
    out
}

/// `y + a k` (or `y - a k`) as in the RK4 substeps.
#[inline]
fn axpy(y: &State, a: f64, backwards: bool, k: &State) -> State {
    let mut out = [0.0; N_STATE];
    for j in 0..N_STATE {
        out[j] = if backwards {
            fnms(a, k[j], y[j])
        } else {
            fma(a, k[j], y[j])
        };
    }
    out
}

/// Evaluation of the dynamics surrogate for fixed mass ratio.
pub(crate) struct Dynamics<'d> {
    data: &'d NRSur7dq4Data,
    q: f64,
    params: MassRatioParameters,
}

impl<'d> Dynamics<'d> {
    pub(crate) fn new(data: &'d NRSur7dq4Data, q: f64) -> Self {
        Self {
            data,
            q,
            params: MassRatioParameters::nrsur7dq4(q),
        }
    }

    fn n_nodes(&self) -> usize {
        self.data.ds_nodes.len()
    }

    fn t(&self, i: usize) -> f64 {
        self.data.t_ds.x[i]
    }

    /// The power table of the fits at the state `y`, and `sin`, `cos` of
    /// the orbital phase.
    #[inline]
    fn table(&self, y: &State) -> (PowerTable, f64, f64) {
        let (sp, cp) = sincos(y[4]);
        (self.params.table(&fit_x(self.q, y, sp, cp)), sp, cp)
    }

    /// `PrecessingNRSur_get_omega`: the orbital frequency fit at a node.
    pub(crate) fn omega(&self, node: usize, y: &State) -> f64 {
        self.data.ds_nodes[node].omega().eval(&self.table(y).0)
    }

    /// `PrecessingNRSur_get_time_deriv_from_index` (with
    /// `PrecessingNRSur_assemble_dydt`) for the fits of `node`.
    fn deriv_at(&self, node: &DynamicsNode, y: &State) -> State {
        let (table, sp, cp) = self.table(y);
        let v: [f64; 9] = std::array::from_fn(|i| node.fits[i].eval(&table));
        let omega = v[0];
        let (oc, ca, cb) = ([v[1], v[2]], [v[3], v[4], v[5]], [v[6], v[7], v[8]]);
        let oqx = fms(oc[0], cp, oc[1] * sp);
        let oqy = fma(oc[0], sp, oc[1] * cp);
        [
            fms(-0.5 * y[1], oqx, 0.5 * y[2] * oqy),
            fma(-0.5 * y[3], oqy, 0.5 * y[0] * oqx),
            fma(0.5 * y[3], oqx, 0.5 * y[0] * oqy),
            fms(0.5 * y[1], oqy, 0.5 * y[2] * oqx),
            omega,
            fms(ca[0], cp, ca[1] * sp),
            fma(ca[0], sp, ca[1] * cp),
            ca[2],
            fms(cb[0], cp, cb[1] * sp),
            fma(cb[0], sp, cb[1] * cp),
            cb[2],
        ]
    }

    /// The derivative at dynamics node `i`.
    fn deriv_node(&self, i: usize, y: &State) -> State {
        self.deriv_at(&self.data.ds_nodes[i], y)
    }

    /// The derivative at the half node `t_{i+1/2}`, `i < 3`.
    fn deriv_half_node(&self, i: usize, y: &State) -> State {
        self.deriv_at(&self.data.ds_half_nodes[i], y)
    }

    /// `PrecessingNRSur_get_time_deriv`: the derivative at any time, by
    /// polynomial interpolation of the derivatives at four nearby nodes.
    fn deriv_time(&self, t: f64, y: &State) -> Result<State, NRSurError> {
        let imax = self.n_nodes() - 1;
        if t < self.t(0) || t > self.t(imax) {
            return Err(NRSurError::Failure(format!(
                "tried to get dydt at t={t:.12}, not in the valid range [{:.12}, {:.12}]",
                self.t(0),
                self.t(imax)
            )));
        }
        let mut i1 = 0;
        while self.t(i1 + 1) < t {
            i1 += 1;
        }
        let i0 = i1.saturating_sub(1).min(imax - 3);
        let times = [self.t(i0), self.t(i0 + 1), self.t(i0 + 2), self.t(i0 + 3)];
        let d = [
            self.deriv_node(i0, y),
            self.deriv_node(i0 + 1, y),
            self.deriv_node(i0 + 2, y),
            self.deriv_node(i0 + 3, y),
        ];
        let mut dydt = [0.0; N_STATE];
        for (j, out) in dydt.iter_mut().enumerate() {
            *out = polynomial4(&times, &[d[0][j], d[1][j], d[2][j], d[3][j]], t);
        }
        Ok(dydt)
    }

    /// `PrecessingNRSur_get_t_ref`: the time at which the orbital frequency
    /// fit (evaluated with the reference state `y0` at every node) equals
    /// `omega_ref`, by linear interpolation between nodes. `omega` caches
    /// the fit at the nodes for this `y0` (NaN when not yet evaluated).
    pub(crate) fn t_ref(
        &self,
        omega_ref: f64,
        y0: &State,
        omega: &mut [f64],
    ) -> Result<f64, NRSurError> {
        if omega_ref.abs() < 1e-10 {
            // LAL warns that omega_ref = 0 selects the start of the surrogate.
            return Ok(self.data.version.start_time());
        }
        if omega_ref > OMEGA_REF_MAX {
            return Err(NRSurError::ReferenceFrequencyTooHigh { omega: omega_ref });
        }
        let mut omega_at = |i: usize| {
            if omega[i].is_nan() {
                omega[i] = self.omega(i, y0);
            }
            omega[i]
        };
        let mut omega_min = omega_at(0);
        if omega_ref < omega_min {
            return Err(NRSurError::FrequencyBelowSurrogate {
                omega: omega_ref,
                omega_min,
            });
        }
        let n = self.n_nodes();
        let mut imax = 1;
        let mut omega_max = omega_at(imax);
        while omega_max <= omega_ref && imax < n {
            imax += 1;
            omega_min = omega_max;
            if imax == n {
                // LAL reads one node past the end here.
                return Err(NRSurError::Failure(format!(
                    "tried all nodes and still have omega={omega_min:.4} <= \
                     omega_ref={omega_ref:.4}"
                )));
            }
            omega_max = omega_at(imax);
        }
        if omega_max <= omega_ref {
            return Err(NRSurError::Failure(format!(
                "tried all nodes and still have omega={omega_min:.4} <= omega_ref={omega_ref:.4}"
            )));
        }
        let t_min = self.t(imax - 1);
        let t_max = self.t(imax);
        Ok(fma(
            t_min,
            omega_max - omega_ref,
            t_max * (omega_ref - omega_min),
        ) / (omega_max - omega_min))
    }

    /// `PrecessingNRSur_initialize_at_dynamics_node`: step from `t_ref` to
    /// the nearest node, which is returned with the state there.
    fn initialize_at_node(
        &self,
        t_ref: f64,
        y_ref: &State,
        norm_a: f64,
        norm_b: f64,
    ) -> Result<(usize, State), NRSurError> {
        let mut imin = 0;
        let mut imax = self.n_nodes() - 1;
        let mut tmin = self.t(imin);
        let mut tmax = self.t(imax);
        while imax - imin > 1 {
            let i0 = (imax + imin) / 2;
            let t0 = self.t(i0);
            if t0 > t_ref {
                imax = i0;
                tmax = t0;
            } else {
                imin = i0;
                tmin = t0;
            }
        }
        let (i0, dt) = if (t_ref - tmin).abs() < (tmax - t_ref).abs() {
            (imin, tmin - t_ref)
        } else {
            (imax, tmax - t_ref)
        };
        let dydt = self.deriv_time(t_ref, y_ref)?;
        let mut y = *y_ref;
        for j in 0..N_STATE {
            y[j] = fma(dt, dydt[j], y[j]);
        }
        normalize_y(norm_a, norm_b, &mut y);
        Ok((i0, y))
    }

    /// `PrecessingNRSur_IntegrateDynamics`: the state at every dynamics node,
    /// written to `out` (length `n_nodes`), from the state at the reference
    /// frequency `omega_ref`. `omega_cache` is as for [`Self::t_ref`].
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn integrate(
        &self,
        chi_a0: &[f64; 3],
        chi_b0: &[f64; 3],
        omega_ref: f64,
        init_orbphase: f64,
        init_quat: &[f64; 4],
        unlimited_extrapolation: bool,
        omega_cache: &mut [f64],
        out: &mut [State],
    ) -> Result<(), NRSurError> {
        let norm_a = fma(
            chi_a0[2],
            chi_a0[2],
            fma(chi_a0[0], chi_a0[0], chi_a0[1] * chi_a0[1]),
        )
        .sqrt();
        let norm_b = fma(
            chi_b0[2],
            chi_b0[2],
            fma(chi_b0[0], chi_b0[0], chi_b0[1] * chi_b0[1]),
        )
        .sqrt();
        let version = self.data.version;
        let q = self.q;
        if q < 0.999 {
            return Err(NRSurError::InvalidMassRatio(q));
        }
        if q > version.q_max() && !unlimited_extrapolation {
            return Err(NRSurError::MassRatioTooLarge {
                q,
                q_max: version.q_max(),
            });
        }
        for (label, norm) in [("A", norm_a), ("B", norm_b)] {
            if norm > version.chi_max() && !unlimited_extrapolation {
                return Err(NRSurError::SpinTooLarge {
                    label,
                    magnitude: norm,
                    chi_max: version.chi_max(),
                });
            }
        }

        let mut y0 = [0.0; N_STATE];
        y0[..4].copy_from_slice(init_quat);
        y0[4] = init_orbphase;
        y0[5..8].copy_from_slice(chi_a0);
        y0[8..11].copy_from_slice(chi_b0);
        let t_ref = self.t_ref(omega_ref, &y0, omega_cache)?;
        let (i0, y_i0) = self.initialize_at_node(t_ref, &y0, norm_a, norm_b)?;
        out[i0] = y_i0;

        let (i_start, time_steps, dydt) = if i0 == 0 {
            self.initialize_rk4_half_nodes(out, norm_a, norm_b)
        } else {
            self.initialize_rk4(out, norm_a, norm_b, i0)?
        };
        self.integrate_ab4(out, &time_steps, &dydt, norm_a, norm_b, i_start);
        Ok(())
    }

    /// `PrecessingNRSur_initialize_RK4_with_half_nodes`: three RK4 steps from
    /// node 0 using the fits at the half nodes for the substeps.
    fn initialize_rk4_half_nodes(
        &self,
        out: &mut [State],
        norm_a: f64,
        norm_b: f64,
    ) -> (usize, [f64; 3], [State; 4]) {
        let mut time_steps = [0.0; 3];
        for (i, ts) in time_steps.iter_mut().enumerate() {
            *ts = self.t(i + 1) - self.t(i);
        }
        let mut dydt = [[0.0; N_STATE]; 4];
        for i in 0..3 {
            let node = out[i];
            let h = time_steps[i];
            dydt[i] = self.deriv_node(i, &node);
            let k2 = self.deriv_half_node(i, &axpy(&node, 0.5 * h, false, &dydt[i]));
            let k3 = self.deriv_half_node(i, &axpy(&node, 0.5 * h, false, &k2));
            let k4 = self.deriv_node(i + 1, &axpy(&node, h, false, &k3));
            let mut next = rk4_combine(&node, h, false, &dydt[i], &k2, &k3, &k4);
            normalize_y(norm_a, norm_b, &mut next);
            out[i + 1] = next;
        }
        dydt[3] = self.deriv_node(3, &out[3]);
        (0, time_steps, dydt)
    }

    /// `PrecessingNRSur_initialize_RK4`: three RK4 steps (backwards if
    /// possible) from node `i0`, with the substeps interpolated in time.
    fn initialize_rk4(
        &self,
        out: &mut [State],
        norm_a: f64,
        norm_b: f64,
        i0: usize,
    ) -> Result<(usize, [f64; 3], [State; 4]), NRSurError> {
        let mut dydt = [[0.0; N_STATE]; 4];
        let mut time_steps = [0.0; 3];
        if i0 > 2 {
            let i_start = i0 - 3;
            for (i, ts) in time_steps.iter_mut().enumerate() {
                *ts = self.t(i_start + i + 1) - self.t(i_start + i);
            }
            for i in 0..3 {
                let node = out[i0 - i];
                let h = time_steps[2 - i];
                dydt[3 - i] = self.deriv_node(i0 - i, &node);
                let t_half = 0.5 * (self.t(i0 - i) + self.t(i0 - i - 1));
                let k2 = self.deriv_time(t_half, &axpy(&node, 0.5 * h, true, &dydt[3 - i]))?;
                let k3 = self.deriv_time(t_half, &axpy(&node, 0.5 * h, true, &k2))?;
                let k4 = self.deriv_node(i0 - i - 1, &axpy(&node, h, true, &k3));
                let mut next = rk4_combine(&node, h, true, &dydt[3 - i], &k2, &k3, &k4);
                normalize_y(norm_a, norm_b, &mut next);
                out[i0 - i - 1] = next;
            }
            dydt[0] = self.deriv_node(i0 - 3, &out[i0 - 3]);
            Ok((i_start, time_steps, dydt))
        } else {
            for (i, ts) in time_steps.iter_mut().enumerate() {
                *ts = self.t(i0 + i + 1) - self.t(i0 + i);
            }
            for i in 0..3 {
                let node = out[i0 + i];
                let h = time_steps[i];
                dydt[i] = self.deriv_node(i0 + i, &node);
                let t_half = 0.5 * (self.t(i + i0) + self.t(i + i0 + 1));
                let k2 = self.deriv_time(t_half, &axpy(&node, 0.5 * h, false, &dydt[i]))?;
                let k3 = self.deriv_time(t_half, &axpy(&node, 0.5 * h, false, &k2))?;
                let k4 = self.deriv_node(i0 + i + 1, &axpy(&node, h, false, &k3));
                let mut next = rk4_combine(&node, h, false, &dydt[i], &k2, &k3, &k4);
                normalize_y(norm_a, norm_b, &mut next);
                out[i0 + i + 1] = next;
            }
            dydt[3] = self.deriv_node(i0 + 3, &out[i0 + 3]);
            Ok((i0, time_steps, dydt))
        }
    }

    /// `PrecessingNRSur_integrate_AB4`: Adams-Bashforth steps forwards from
    /// node `i_start + 3` and backwards from `i_start`.
    fn integrate_ab4(
        &self,
        out: &mut [State],
        time_steps: &[f64; 3],
        dydt: &[State; 4],
        norm_a: f64,
        norm_b: f64,
        i_start: usize,
    ) {
        let n = self.n_nodes();
        let (mut dt1, mut dt2, mut dt3) = (time_steps[0], time_steps[1], time_steps[2]);
        let [mut k1, mut k2, mut k3, mut k4] = *dydt;
        for i in i_start + 4..n {
            let dt4 = self.t(i) - self.t(i - 1);
            let mut next = ab4_dy(&k1, &k2, &k3, &k4, dt1, dt2, dt3, dt4);
            for (y, prev) in next.iter_mut().zip(&out[i - 1]) {
                *y += prev;
            }
            normalize_y(norm_a, norm_b, &mut next);
            out[i] = next;
            if i < n - 1 {
                (dt1, dt2, dt3) = (dt2, dt3, dt4);
                (k1, k2, k3) = (k2, k3, k4);
                k4 = self.deriv_node(i, &out[i]);
            }
        }

        let (mut dt1, mut dt2, mut dt3) = (-time_steps[2], -time_steps[1], -time_steps[0]);
        let [mut k4, mut k3, mut k2, mut k1] = *dydt;
        for i in (0..i_start).rev() {
            let dt4 = self.t(i) - self.t(i + 1);
            let mut next = ab4_dy(&k1, &k2, &k3, &k4, dt1, dt2, dt3, dt4);
            for (y, prev) in next.iter_mut().zip(&out[i + 1]) {
                *y += prev;
            }
            normalize_y(norm_a, norm_b, &mut next);
            out[i] = next;
            if i > 0 {
                (dt1, dt2, dt3) = (dt2, dt3, dt4);
                (k1, k2, k3) = (k2, k3, k4);
                k4 = self.deriv_node(i, &out[i]);
            }
        }
    }
}

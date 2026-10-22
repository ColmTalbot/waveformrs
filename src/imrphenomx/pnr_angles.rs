// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2024-2026 Colm Talbot. Parts of waveformrs are translated
// from LALSuite, GSL and fdlibm: see NOTICE.md for their authors.
//! The model-independent parts of the PNR precession angles of
//! IMRPhenomXPNR (Hamilton et al. 2021, arXiv:2107.08876), following
//! `LALSimIMRPhenomX_PNR_alpha.c`, `LALSimIMRPhenomX_PNR_beta.c` and
//! `LALSimIMRPhenomX_PNR_internals.c`:
//!
//! - the merger-ringdown (MR) expressions of `alpha` and `beta` and their
//!   connection to the inspiral (C1 interpolation, rescaling of the inspiral
//!   `beta`, connection frequencies from the turning points of the MR
//!   `beta`);
//! - the "waveform" `beta` of the inspiral (Eq. 41 of arXiv:2107.08876);
//! - the windows that blend PNR and PN angles outside the calibration
//!   region and switch off the co-precessing deviations;
//! - the linear frequency map of the angles of the higher modes.
//!
//! Products and sums are fused where the hardware has fused multiply-adds
//! (see [`crate::spintaylor::fp`]). LAL computes the inflection points of the
//! MR `beta` in complex arithmetic and passes complex arguments to `acos` and
//! `cos`, which C converts to their real parts; the same is done here, with
//! a scaled complex division and the principal square root.

use std::f64::consts::PI;

use crate::imrphenomx::precession::{atan2_tol, MAX_TOL_ATAN};
use crate::spintaylor::fp::{fma, fnms};

/// Calibration region of the PNR angles: `q <= 8.5` and `chi <= 0.85`.
pub(crate) const PNR_Q_WINDOW_LOWER: f64 = 8.5;
pub(crate) const PNR_CHI_WINDOW_LOWER: f64 = 0.85;
/// End of the transition to the PN angles: `q <= 12` and `chi <= 1.2`.
pub(crate) const PNR_Q_WINDOW_UPPER: f64 = 12.0;
pub(crate) const PNR_CHI_WINDOW_UPPER: f64 = 1.2;

/// Frequency step of the finite differences of the connections (`dMf`).
pub(crate) const PNR_DMF: f64 = 0.0005;

/// `IMRPhenomX_PNR_alpha_parameters`: the coefficients of the MR `alpha`
/// (Eq. 48 of arXiv:2107.08876) and of its connection to the inspiral.
#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct AlphaParams {
    pub(crate) a1: f64,
    pub(crate) a2: f64,
    pub(crate) a3: f64,
    pub(crate) a4: f64,
    pub(crate) mf_alpha_lower: f64,
    pub(crate) mf_alpha_upper: f64,
    pub(crate) alpha_lower: f64,
    pub(crate) alpha_upper: f64,
    pub(crate) derivative_alpha_lower: f64,
    pub(crate) derivative_alpha_upper: f64,
    pub(crate) alpha_interp: [f64; 4],
    pub(crate) alpha_mr_offset: f64,
}

/// `IMRPhenomX_PNR_beta_parameters`: the coefficients of the MR `beta`
/// (Eq. 49 of arXiv:2107.08876), the connection frequencies, the rescaling
/// of the inspiral `beta` and (version 330) the cubic connecting the
/// SpinTaylor `beta` to the MSA `beta`.
#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct BetaParams {
    pub(crate) b0: f64,
    pub(crate) b1: f64,
    pub(crate) b2: f64,
    pub(crate) b3: f64,
    pub(crate) b4: f64,
    pub(crate) b5: f64,
    pub(crate) mf_beta_lower: f64,
    pub(crate) mf_beta_upper: f64,
    pub(crate) beta_lower: f64,
    pub(crate) beta_upper: f64,
    pub(crate) derivative_beta_lower: f64,
    pub(crate) derivative_beta_upper: f64,
    pub(crate) beta_rescale_1: f64,
    pub(crate) beta_rescale_2: f64,
    pub(crate) beta_interp: [f64; 4],
}

// ---------------------------------------------------------------------------
// alpha
// ---------------------------------------------------------------------------

/// `IMRPhenomX_PNR_alpha_interpolation_0`, Eq. 52 of arXiv:2107.08876.
pub(crate) fn alpha_interpolation_0(
    mf1: f64,
    mf2: f64,
    alpha1: f64,
    alpha2: f64,
    dalpha1: f64,
    dalpha2: f64,
) -> f64 {
    let d = (mf2 - mf1) * (mf2 - mf1) * (mf2 - mf1);
    let n = fma(
        2.0,
        fma(mf1, alpha1, -(mf2 * alpha2)),
        -((mf1 - mf2) * (fma(mf1, dalpha1, mf2 * dalpha2) + (alpha1 + alpha2))),
    );
    n / d
}

/// `IMRPhenomX_PNR_alpha_interpolation_1`.
pub(crate) fn alpha_interpolation_1(
    mf1: f64,
    mf2: f64,
    alpha1: f64,
    alpha2: f64,
    dalpha1: f64,
    dalpha2: f64,
) -> f64 {
    let d = (mf2 - mf1) * (mf2 - mf1) * (mf2 - mf1);
    let n = fma(
        3.0 * (mf1 + mf2),
        fma(mf2, alpha2, -(mf1 * alpha1)),
        (mf1 - mf2)
            * fma(
                fma(2.0, mf2, mf1),
                fma(mf1, dalpha1, alpha1),
                fma(2.0, mf1, mf2) * fma(mf2, dalpha2, alpha2),
            ),
    );
    n / d
}

/// `IMRPhenomX_PNR_alpha_interpolation_2`.
pub(crate) fn alpha_interpolation_2(
    mf1: f64,
    mf2: f64,
    alpha1: f64,
    alpha2: f64,
    dalpha1: f64,
    dalpha2: f64,
) -> f64 {
    let d = (mf2 - mf1) * (mf2 - mf1) * (mf2 - mf1);
    let n = fma(
        6.0 * mf1 * mf2,
        fma(mf1, alpha1, -(mf2 * alpha2)),
        -((mf1 - mf2)
            * fma(
                mf2 * fma(2.0, mf1, mf2),
                fma(mf1, dalpha1, alpha1),
                mf1 * fma(2.0, mf2, mf1) * fma(mf2, dalpha2, alpha2),
            )),
    );
    n / d
}

/// `IMRPhenomX_PNR_alpha_interpolation_3`.
pub(crate) fn alpha_interpolation_3(
    mf1: f64,
    mf2: f64,
    alpha1: f64,
    alpha2: f64,
    dalpha1: f64,
    dalpha2: f64,
) -> f64 {
    let d = (mf2 - mf1) * (mf2 - mf1) * (mf2 - mf1);
    let n = fma(
        mf1 * mf2 * (mf1 - mf2),
        fma(
            mf2,
            fma(mf1, dalpha1, alpha1),
            mf1 * fma(mf2, dalpha2, alpha2),
        ),
        fma(
            mf1 * mf2 * mf2 * fma(-3.0, mf1, mf2),
            alpha1,
            -(mf1 * mf1 * mf2 * fma(-3.0, mf2, mf1) * alpha2),
        ),
    );
    n / d
}

/// `IMRPhenomX_PNR_intermediate_alpha_expression`, Eq. 51 of
/// arXiv:2107.08876.
#[inline]
pub(crate) fn intermediate_alpha(mf: f64, p: &AlphaParams) -> f64 {
    let [a0, a1, a2, a3] = p.alpha_interp;
    fma(a0 * mf, mf, a1 * mf) + a2 + a3 / mf
}

/// `IMRPhenomX_PNR_MR_alpha_expression`, Eq. 48 of arXiv:2107.08876.
#[inline]
pub(crate) fn mr_alpha(mf: f64, p: &AlphaParams) -> f64 {
    -(p.a1 / mf + (p.a2 * p.a3.sqrt()) / fma(mf - p.a4, mf - p.a4, p.a3))
}

// ---------------------------------------------------------------------------
// beta
// ---------------------------------------------------------------------------

/// `IMRPhenomX_PNR_MR_beta_expression`, Eq. 49 of arXiv:2107.08876.
#[inline]
pub(crate) fn mr_beta(mf: f64, p: &BetaParams) -> f64 {
    p.b0 + fma(p.b3 * mf, mf, fma(p.b2, mf, p.b1)) / fma(p.b4 * (mf + p.b5), mf + p.b5, 1.0)
}

/// `IMRPhenomX_PNR_MR_dbeta_expression`.
pub(crate) fn mr_dbeta(mf: f64, p: &BetaParams) -> f64 {
    let den = fma(p.b4 * (mf + p.b5), mf + p.b5, 1.0);
    (fma(
        fma(2.0 * p.b3 * p.b4, p.b5, -(p.b2 * p.b4)) * mf,
        mf,
        fma(
            2.0 * p.b3 * p.b4 * p.b5,
            p.b5,
            fma(2.0, p.b3, -(2.0 * p.b1 * p.b4)),
        ) * mf,
    ) + fma(
        p.b2 * p.b4 * p.b5,
        p.b5,
        fma(-(2.0 * p.b1 * p.b4), p.b5, p.b2),
    )) / (den * den)
}

/// `IMRPhenomX_PNR_MR_ddbeta_expression`.
pub(crate) fn mr_ddbeta(mf: f64, p: &BetaParams) -> f64 {
    let (b1, b2, b3, b4, b5) = (p.b1, p.b2, p.b3, p.b4, p.b5);
    let a = fma(b2 * b4, b4, -(2.0 * b3 * b4 * b4 * b5));
    let b = fma(
        -(3.0 * b3 * b4 * b4 * b5),
        b5,
        fma(-3.0 * b3, b4, 3.0 * b1 * b4 * b4),
    );
    let c = fma(
        -(3.0 * b2 * b4 * b4 * b5),
        b5,
        fma(-3.0 * b2, b4, 6.0 * b1 * b4 * b4 * b5),
    );
    let d = fma(
        b3 * b4 * b4 * b5 * b5 * b5,
        b5,
        fma(
            -(2.0 * b2 * b4 * b4 * b5 * b5),
            b5,
            fma(
                3.0 * b1 * b4 * b4 * b5,
                b5,
                fma(
                    2.0 * b3 * b4 * b5,
                    b5,
                    fma(-(2.0 * b2 * b4), b5, fma(-b1, b4, b3)),
                ),
            ),
        ),
    );
    2.0 * (fma(c, mf, fma(a * mf * mf, mf, b * mf * mf)) + d)
        / fma(b4 * (b5 + mf), b5 + mf, 1.0).powf(3.0)
}

/// `IMRPhenomX_PNR_MR_dddbeta_expression`.
pub(crate) fn mr_dddbeta(mf: f64, p: &BetaParams) -> f64 {
    let (b1, b2, b3, b4, b5) = (p.b1, p.b2, p.b3, p.b4, p.b5);
    let a = fma(-b2 * b4, b4, 2.0 * b3 * b4 * b4 * b5);
    let b = fma(
        4.0 * b3 * b4 * b4 * b5,
        b5,
        fma(4.0 * b3, b4, -(4.0 * b1 * b4 * b4)),
    );
    let c = fma(
        6.0 * b2 * b4 * b4 * b5,
        b5,
        fma(6.0 * b2, b4, -(12.0 * b1 * b4 * b4 * b5)),
    );
    let d = fma(
        -(4.0 * b3 * b4 * b4 * b5 * b5 * b5),
        b5,
        fma(
            8.0 * b2 * b4 * b4 * b5 * b5,
            b5,
            fma(
                -(12.0 * b1 * b4 * b4 * b5),
                b5,
                fma(
                    -(8.0 * b3 * b4 * b5),
                    b5,
                    fma(8.0 * b2 * b4, b5, fma(-4.0, b3, 4.0 * b1 * b4)),
                ),
            ),
        ),
    );
    let e = fma(
        -(2.0 * b3 * b4 * b4 * b5 * b5 * b5 * b5),
        b5,
        fma(
            3.0 * b2 * b4 * b4 * b5 * b5 * b5,
            b5,
            fma(
                -(4.0 * b1 * b4 * b4 * b5 * b5),
                b5,
                fma(
                    -(4.0 * b3 * b4 * b5 * b5),
                    b5,
                    fma(
                        2.0 * b2 * b4 * b5,
                        b5,
                        fma(4.0 * b1 * b4, b5, fma(-(2.0 * b3), b5, -b2)),
                    ),
                ),
            ),
        ),
    );
    6.0 * b4
        * (fma(
            d,
            mf,
            fma(c * mf, mf, fma(a * mf * mf * mf, mf, b * mf * mf * mf)),
        ) + e)
        / fma(b4 * (b5 + mf), b5 + mf, 1.0).powf(4.0)
}

/// `IMRPhenomX_PNR_beta_rescaling_1`, Eq. 54 of arXiv:2107.08876.
pub(crate) fn beta_rescaling_1(mf: f64, beta1: f64, beta2: f64, dbeta1: f64, dbeta2: f64) -> f64 {
    let d = beta1 * beta1 * mf;
    let n = fma(
        -2.0 * beta1,
        beta2 - beta1,
        fma(beta1, dbeta2, -(beta2 * dbeta1)) * mf,
    );
    -n / d
}

/// `IMRPhenomX_PNR_beta_rescaling_2`, Eq. 55 of arXiv:2107.08876.
pub(crate) fn beta_rescaling_2(mf: f64, beta1: f64, beta2: f64, dbeta1: f64, dbeta2: f64) -> f64 {
    let d = beta1 * mf * (beta1 * mf);
    let n = fma(
        beta1,
        beta2 - beta1,
        -(fma(beta1, dbeta2, -(beta2 * dbeta1)) * mf),
    );
    -n / d
}

/// `IMRPhenomX_PNR_rescale_beta_expression`, Eq. 53 of arXiv:2107.08876.
#[inline]
pub(crate) fn rescale_beta(mf: f64, p: &BetaParams) -> f64 {
    fma(p.beta_rescale_2 * mf, mf, fma(p.beta_rescale_1, mf, 1.0))
}

/// `IMRPhenomX_ST_PNR_beta_coeffs_0`: cubic connecting the SpinTaylor and
/// PNR `beta` (version 330).
pub(crate) fn st_beta_coeffs_0(
    mf1: f64,
    mf2: f64,
    beta1: f64,
    beta2: f64,
    dbeta1: f64,
    dbeta2: f64,
) -> f64 {
    let mf1_2 = mf1 * mf1;
    let mf1_3 = mf1_2 * mf1;
    let mf2_2 = mf2 * mf2;
    let mf2_3 = mf2_2 * mf2;
    let n = fma(
        -(dbeta1 * mf1),
        mf2_3,
        fma(
            beta1,
            mf2_3,
            fma(
                -(dbeta2 * mf1_2),
                mf2_2,
                fma(
                    dbeta1 * mf1_2,
                    mf2_2,
                    fma(
                        -(3.0 * beta1 * mf1),
                        mf2_2,
                        fma(
                            dbeta2 * mf1_3,
                            mf2,
                            fma(-beta2, mf1_3, 3.0 * beta2 * mf1_2 * mf2),
                        ),
                    ),
                ),
            ),
        ),
    );
    let d = (mf1 - mf2) * (mf1 - mf2) * (mf1 - mf2);
    -n / d
}

/// `IMRPhenomX_ST_PNR_beta_coeffs_1`.
pub(crate) fn st_beta_coeffs_1(
    mf1: f64,
    mf2: f64,
    beta1: f64,
    beta2: f64,
    dbeta1: f64,
    dbeta2: f64,
) -> f64 {
    let mf1_2 = mf1 * mf1;
    let mf1_3 = mf1_2 * mf1;
    let mf2_2 = mf2 * mf2;
    let mf2_3 = mf2_2 * mf2;
    let n = fma(
        dbeta1,
        mf2_3,
        fma(
            2.0 * dbeta2 * mf1,
            mf2_2,
            fma(
                dbeta1 * mf1,
                mf2_2,
                fma(
                    -(dbeta2 * mf1_2),
                    mf2,
                    fma(
                        -(2.0 * dbeta1 * mf1_2),
                        mf2,
                        fma(
                            -(6.0 * beta2 * mf1),
                            mf2,
                            fma(-dbeta2, mf1_3, 6.0 * beta1 * mf1 * mf2),
                        ),
                    ),
                ),
            ),
        ),
    );
    let d = (mf1 - mf2) * (mf1 - mf2) * (mf1 - mf2);
    -n / d
}

/// `IMRPhenomX_ST_PNR_beta_coeffs_2`.
pub(crate) fn st_beta_coeffs_2(
    mf1: f64,
    mf2: f64,
    beta1: f64,
    beta2: f64,
    dbeta1: f64,
    dbeta2: f64,
) -> f64 {
    let mf1_2 = mf1 * mf1;
    let mf2_2 = mf2 * mf2;
    let n = fma(
        -fma(2.0, dbeta1, dbeta2),
        mf2_2,
        fma(
            (dbeta1 - dbeta2) * mf1,
            mf2,
            fma(
                -3.0 * (beta1 - beta2),
                mf1 + mf2,
                fma(2.0, dbeta2, dbeta1) * mf1_2,
            ),
        ),
    );
    let d = (mf1 - mf2) * (mf1 - mf2) * (mf1 - mf2);
    -n / d
}

/// `IMRPhenomX_ST_PNR_beta_coeffs_3`.
pub(crate) fn st_beta_coeffs_3(
    mf1: f64,
    mf2: f64,
    beta1: f64,
    beta2: f64,
    dbeta1: f64,
    dbeta2: f64,
) -> f64 {
    let n = fma(2.0, beta1 - beta2, -((dbeta1 + dbeta2) * (mf1 - mf2)));
    let d = (mf1 - mf2) * (mf1 - mf2) * (mf1 - mf2);
    -n / d
}

/// `beta_connection` of `LALSimIMRPhenomX_precession.c`: the cubic between
/// the SpinTaylor and the PNR `beta`.
#[inline]
pub(crate) fn beta_connection(mf: f64, p: &BetaParams) -> f64 {
    let [b0, b1, b2, b3] = p.beta_interp;
    fma(b3 * mf * mf, mf, fma(b2 * mf, mf, fma(b1, mf, b0)))
}

/// `IMRPhenomX_PNR_arctan_window`, Eq. 62 of arXiv:2107.08876: keep `beta`
/// away from 0 and pi.
pub(crate) fn arctan_window(beta: f64) -> f64 {
    let window_border = 0.01;
    if beta <= window_border || beta >= PI - window_border {
        let p = 0.002;
        // (LAL's literal is two ulps above pi/2)
        let pi_by_2 = std::f64::consts::FRAC_PI_2;
        let pi_by_2_1mp = 1.569378278348018;
        let pi_by_2_1oq = 7.308338225719002e97;
        let sign = if beta < pi_by_2 { -1.0 } else { 1.0 };
        return fma(
            sign * pi_by_2_1mp,
            (beta - pi_by_2).powf(1.0 / p).atan2(pi_by_2_1oq).powf(p),
            pi_by_2,
        );
    }
    beta
}

/// `IMRPhenomX_PNR_AttachMRBeta`: whether the MR `beta` is attached.
pub(crate) fn attach_mr_beta(p: &BetaParams) -> bool {
    !(p.mf_beta_lower < 0.009
        || p.mf_beta_lower == 100.0
        || p.beta_upper <= 0.0
        || p.beta_upper > 5.0 * (p.b0 + 0.1))
}

/// A complex number with the arithmetic of C's `double complex` (scaled
/// divisions, the algorithm of compiler-rt's `__divdc3`).
#[derive(Clone, Copy, Debug)]
struct Cplx {
    re: f64,
    im: f64,
}

impl Cplx {
    fn real(re: f64) -> Self {
        Self { re, im: 0.0 }
    }

    /// `real * complex`, componentwise.
    fn scale(self, x: f64) -> Self {
        Self {
            re: x * self.re,
            im: x * self.im,
        }
    }

    fn mul(self, o: Self) -> Self {
        Self {
            re: self.re * o.re - self.im * o.im,
            im: self.re * o.im + self.im * o.re,
        }
    }

    fn sub(self, o: Self) -> Self {
        Self {
            re: self.re - o.re,
            im: self.im - o.im,
        }
    }

    fn add(self, o: Self) -> Self {
        Self {
            re: self.re + o.re,
            im: self.im + o.im,
        }
    }

    /// `complex / real`, componentwise.
    fn div_real(self, x: f64) -> Self {
        Self {
            re: self.re / x,
            im: self.im / x,
        }
    }

    /// `complex / complex`: compiler-rt's `__divdc3` (scaling by `logb` of
    /// the denominator; the NaN recovery branches are not needed here).
    fn div(self, o: Self) -> Self {
        let (a, b) = (self.re, self.im);
        let (mut c, mut d) = (o.re, o.im);
        let logbw = c.abs().max(d.abs()).log2().floor();
        let mut ilogbw = 0;
        if logbw.is_finite() {
            ilogbw = logbw as i32;
            c = scalbn(c, -ilogbw);
            d = scalbn(d, -ilogbw);
        }
        let denom = fma(c, c, d * d);
        Self {
            re: scalbn(fma(a, c, b * d) / denom, -ilogbw),
            im: scalbn(fma(b, c, -(a * d)) / denom, -ilogbw),
        }
    }

    /// C99 `csqrt` for an argument with vanishing imaginary part.
    fn sqrt_of_real_axis(self) -> Self {
        let x = self.re;
        if x.is_nan() {
            Self {
                re: f64::NAN,
                im: f64::NAN,
            }
        } else if x >= 0.0 {
            Self {
                re: x.sqrt(),
                im: self.im,
            }
        } else {
            Self {
                re: 0.0,
                im: (-x).sqrt().copysign(self.im),
            }
        }
    }
}

/// `scalbn` (exact scaling by a power of two for normal results).
fn scalbn(x: f64, n: i32) -> f64 {
    // Split the exponent so that the power of two is representable.
    let mut x = x;
    let mut n = n;
    while n > 1000 {
        x *= 2f64.powi(1000);
        n -= 1000;
    }
    while n < -1000 {
        x *= 2f64.powi(-1000);
        n += 1000;
    }
    x * 2f64.powi(n)
}

/// The cubic `a x^3 + b x^2 + c x + d` whose roots are the inflection
/// points of the MR `beta` (Eq. 68 of arXiv:2107.08876).
fn inflection_cubic(p: &BetaParams) -> [f64; 4] {
    let (b1, b2, b3, b4, b5) = (p.b1, p.b2, p.b3, p.b4, p.b5);
    let a = 2.0 * fma(b2 * b4, b4, -(2.0 * b3 * b4 * b4 * b5));
    let b = 6.0 * fma(-(b3 * b4 * b4 * b5), b5, fma(-b3, b4, b1 * b4 * b4));
    let c = 6.0
        * fma(
            -(b2 * b4 * b4 * b5),
            b5,
            fma(-b2, b4, 2.0 * b1 * b4 * b4 * b5),
        );
    let d = 2.0
        * fma(
            b3 * b4 * b4 * b5 * b5 * b5,
            b5,
            fma(
                -(2.0 * b2 * b4 * b4 * b5 * b5),
                b5,
                fma(
                    3.0 * b1 * b4 * b4 * b5,
                    b5,
                    fma(
                        2.0 * b3 * b4 * b5,
                        b5,
                        fma(-(2.0 * b2 * b4), b5, fma(-b1, b4, b3)),
                    ),
                ),
            ),
        );
    [a, b, c, d]
}

/// `IMRPhenomX_PNR_three_inflection_points`: the roots of the depressed
/// cubic with the trigonometric formula, in complex arithmetic.
fn three_inflection_points(p: &BetaParams) -> [Cplx; 3] {
    let [a, b, c, d] = inflection_cubic(p);
    let (ca, cb, cc, cd) = (Cplx::real(a), Cplx::real(b), Cplx::real(c), Cplx::real(d));
    let three_a = ca.scale(3.0);
    // r = (3ac - b^2) / (3a^2), s = (2b^3 - 9abc + 27a^2 d) / (27a^3)
    let r = three_a.mul(cc).sub(cb.mul(cb)).div(three_a.mul(ca));
    let s_num = cb
        .scale(2.0)
        .mul(cb)
        .mul(cb)
        .sub(ca.scale(9.0).mul(cb).mul(cc))
        .add(ca.scale(27.0).mul(ca).mul(cd));
    let s = s_num.div(ca.scale(27.0).mul(ca).mul(ca));
    // phi = acos(Re[(3 s / (2 r)) sqrt(-3 / r)])
    let q1 = s.scale(3.0).div(r.scale(2.0));
    let w1 = Cplx::real(-3.0).div(r).sqrt_of_real_axis();
    let phi = q1.mul(w1).re.acos();
    // 2 sqrt(-r / 3) cos((phi - 2 n pi) / 3) - b / (3a), with -r/3 = r/(-3)
    let two_w2 = r.div_real(-3.0).sqrt_of_real_axis().scale(2.0);
    let q = cb.div(three_a);
    let shifts = [0.0, -2.0 * PI, -4.0 * PI];
    shifts.map(|shift| {
        let cosv = if shift == 0.0 {
            (phi / 3.0).cos()
        } else {
            ((phi + shift) / 3.0).cos()
        };
        Cplx {
            re: cosv * two_w2.re,
            im: cosv * two_w2.im,
        }
        .sub(q)
    })
}

/// `IMRPhenomX_PNR_single_inflection_point`: the inflection frequency used
/// for the connection of `beta` (Eq. 69 of arXiv:2107.08876).
pub(crate) fn single_inflection_point(p: &BetaParams) -> f64 {
    let [a, b, c, d] = inflection_cubic(p);
    let mut f_inflection = 0.0;
    if a.abs() < 1e-10 {
        if b.abs() < 2e-10 {
            f_inflection = -d / c;
        } else {
            // IMRPhenomX_PNR_two_inflection_points
            let w = Cplx::real(fma(c, c, -(4.0 * b * d))).sqrt_of_real_axis();
            let k = 1.0 / (2.0 * b);
            let roots = [
                Cplx {
                    re: -c + w.re,
                    im: w.im,
                }
                .scale(k),
                Cplx {
                    re: -c - w.re,
                    im: -w.im,
                }
                .scale(k),
            ];
            // choose the inflection point with negative gradient
            for root in roots {
                if mr_dbeta(root.re, p) < 0.0 {
                    f_inflection = root.re;
                }
            }
        }
    } else {
        let roots = three_inflection_points(p);
        let mut w = 0;
        let mut f_temp = 0.0;
        // a substantial imaginary part (sic: without fabs) marks a complex root
        for root in roots {
            if root.im < 1e-10 {
                w += 1;
                f_temp = root.re;
            }
        }
        if w == 1 {
            f_inflection = f_temp;
        } else if a < 0.0 {
            // take the central root
            f_inflection = roots[1].re;
        } else if b / (3.0 * a) > p.b5 / 2.0 - (2141.0 / 90988.0) {
            f_inflection = roots[0].re;
        } else {
            f_inflection = roots[2].re;
        }
    }
    f_inflection
}

/// `IMRPhenomX_PNR_BetaConnectionFrequencies`: the lower and upper
/// connection frequencies of `beta` (Sec. 8B of arXiv:2107.08876).
pub(crate) fn beta_connection_frequencies(p: &mut BetaParams) {
    let (b1, b2, b3, b4, b5) = (p.b1, p.b2, p.b3, p.b4, p.b5);
    let mf_inflection = single_inflection_point(p);
    let dbeta_inflection = mr_dbeta(mf_inflection, p);

    // the two extrema of beta
    let k = fma(b3 * b4 * b5, b5, fma(-b1, b4, b3));
    let root_term = fma(
        b4 * fma(-(2.0 * b3), b5, b2),
        fma(b2 * b4 * b5, b5, fma(-(2.0 * b1 * b4), b5, b2)),
        k * k,
    );
    let den = b4 * fma(-(2.0 * b3), b5, b2);
    let mf_plus = (k + root_term.sqrt()) / den;
    let mf_minus = (k - root_term.sqrt()) / den;

    let (mf_at_minimum, mf_at_maximum) = if mr_ddbeta(mf_plus, p) > 0.0 {
        (mf_plus, mf_minus)
    } else {
        (mf_minus, mf_plus)
    };
    let ddbeta = mr_ddbeta(mf_at_maximum, p);
    let dddbeta = mr_dddbeta(mf_at_maximum, p);
    let sign = if dbeta_inflection > 0.0 { 1.0 } else { -1.0 };
    let chosen_dbeta = sign * (dbeta_inflection / 100.0 * (dbeta_inflection / 100.0)) * 25.0;
    let root = fma(ddbeta, ddbeta, 2.0 * dddbeta * chosen_dbeta).sqrt();
    let d1 = 1.0 / dddbeta * (-ddbeta + root);
    let d2 = 1.0 / dddbeta * (-ddbeta - root);
    let delta = if d1 > 0.0 {
        if d2 > 0.0 {
            if d1 < d2 {
                d1
            } else {
                d2
            }
        } else {
            d1
        }
    } else {
        d2
    };

    let mf_low = if mf_inflection >= 0.06 {
        mf_inflection - 0.03
    } else {
        3.0 * mf_inflection / 5.0
    };
    let mf_im = if mf_at_minimum > mf_at_maximum || mf_inflection > mf_at_maximum {
        if mf_at_maximum >= mf_low {
            mf_at_maximum + delta
        } else {
            mf_low
        }
    } else if mf_at_minimum > 0.06 {
        mf_at_minimum - 0.03
    } else {
        3.0 * mf_at_minimum / 5.0
    };
    let mf_mr = if mf_at_minimum > mf_inflection {
        mf_at_minimum
    } else {
        100.
    };
    if mf_im < 0.0 || mf_im.is_nan() {
        p.mf_beta_lower = 100.0;
        p.mf_beta_upper = 100.0;
    } else {
        p.mf_beta_lower = mf_im;
        p.mf_beta_upper = mf_mr;
    }
}

// ---------------------------------------------------------------------------
// "waveform" beta of the inspiral
// ---------------------------------------------------------------------------

/// `IMRPhenomX_PNR_chi_calc`: the effective spin magnitude from the
/// frequency-dependent in-plane spin (Eqs. 47 and 18 of arXiv:2107.08876).
pub(crate) fn chi_calc(m1: f64, l: f64, j0: f64, l0: f64, chi_parr: f64, beta: f64) -> f64 {
    let s_perp = (j0 - (l0 - l)) * beta.sin();
    let chi_p = s_perp / (m1 * m1);
    fma(chi_parr, chi_parr, chi_p * chi_p).sqrt()
}

/// `IMRPhenomX_PNR_PNWaveformBeta`, Eq. 41 of arXiv:2107.08876: the opening
/// angle of the optimal emission direction from the dynamical `iota`.
pub(crate) fn pn_waveform_beta(
    mf: f64,
    iota: f64,
    m1: f64,
    m2: f64,
    chi: f64,
    costheta: f64,
) -> f64 {
    let m = m1 + m2;
    let nu = (m1 * m2) / (m * m);
    let delta = (m1 - m2) / m;
    let theta = costheta.acos();
    let w_orb = PI * mf;
    let v = w_orb.cbrt();
    let v2 = v * v;
    let v3 = v * v2;
    let cos_iota = iota.cos();
    let sin_iota = iota.sin();
    let cos_half_iota = (iota / 2.0).cos();
    let sin_half_iota = (iota / 2.0).sin();
    let cos_theta = theta.cos();
    let sin_theta = theta.sin();

    // numerator expansion of h21/h22
    let n0 = 84.0 * sin_iota;
    let n2 = 2.0 * fma(55.0, nu, -107.0) * sin_iota;
    let n3 = fma(
        -7.0 * (fma(5.0, nu, 6.0 * delta) + 6.0) * chi * fma(2.0, cos_iota, -1.0),
        sin_theta,
        56.0 * fma(-((1.0 + delta - nu) * chi), cos_theta, 3.0 * PI) * sin_iota,
    );
    let n = fma(n3, v3, fma(n2, v2, n0)) / cos_half_iota;

    // denominator expansion of h21/h22
    let d0 = 84.0 * cos_half_iota;
    let d2 = 2.0 * fma(55.0, nu, -107.0) * cos_half_iota;
    let d3 = fma(
        56.0 * fma((nu - 1.0 - delta) * chi, cos_theta, 3.0 * PI),
        cos_half_iota,
        14.0 * fma(5.0, nu, fma(6.0, delta, 6.0)) * chi * sin_theta * sin_half_iota,
    );
    let d = fma(d3, v3, fma(d2, v2, d0));
    2.0 * atan2_tol(n, 2.0 * d, MAX_TOL_ATAN)
}

// ---------------------------------------------------------------------------
// windows and frequency maps
// ---------------------------------------------------------------------------

/// `0.5 cos(pi x) + 0.5` for `0 < x <= 1`, 0 above and 1 below.
fn cosine_window(argument: f64) -> f64 {
    if argument > 0.0 && argument <= 1.0 {
        fma(0.5, (argument * PI).cos(), 0.5)
    } else if argument > 1.0 {
        0.0
    } else {
        1.0
    }
}

/// `IMRPhenomX_PNR_AnglesWindow`: blending of the PNR and PN angles in the
/// transition region (`8.5 < q <= 12`, `0.85 < chi <= 1.2`).
pub(crate) fn angles_window(q: f64, chi_single_spin: f64) -> f64 {
    let window_q_argument = (q - 8.5) / 3.5;
    let window_chi_argument = (chi_single_spin - 0.85) / 0.35;
    cosine_window(window_q_argument) * cosine_window(window_chi_argument)
}

/// `IMRPhenomX_PNR_CoprecWindow`: the co-precessing deviations are switched
/// off between `q = 10` and `q = 20`.
pub(crate) fn coprec_window(q: f64) -> f64 {
    cosine_window((q - 10.0) / 10.0)
}

/// `IMRPhenomX_PNR_LinearFrequencySlope`.
pub(crate) fn linear_frequency_slope(
    emm: f64,
    mf_lower: f64,
    mf_upper: f64,
    mf_rd_22: f64,
    mf_rd_lm: f64,
) -> f64 {
    (mf_upper - (mf_rd_lm - mf_rd_22) - 2.0 * mf_lower / emm) / (mf_upper - mf_lower)
}

/// `IMRPhenomX_PNR_LinearFrequencyMap`: the frequency of the (2,2) angles
/// used for mode `(l, m)` (inspiral scaling `2 Mf / m`, ringdown shift by
/// the difference of the ringdown frequencies, linear in between). Returns
/// `None` where LAL fails (`Mf_lower > Mf_upper` in the transition).
#[allow(clippy::too_many_arguments)]
pub(crate) fn linear_frequency_map(
    mf: f64,
    ell: u32,
    emm: u32,
    mf_lower: f64,
    mf_upper: f64,
    mf_rd_22: f64,
    mf_rd_lm: f64,
    inspiral: bool,
) -> Option<f64> {
    if ell == 2 && emm == 2 {
        return Some(mf);
    }
    let emm = emm as f64;
    if mf <= mf_lower || inspiral {
        return Some(2.0 * mf / emm);
    }
    if mf_lower > mf_upper {
        return None;
    }
    if mf > mf_upper {
        return Some(mf - (mf_rd_lm - mf_rd_22));
    }
    let a = linear_frequency_slope(emm, mf_lower, mf_upper, mf_rd_22, mf_rd_lm);
    let b = 2.0 * mf_lower / emm;
    Some(fma(a, mf - mf_lower, b))
}

/// `IMRPhenomX_PNR_LinearFrequencyMapTransitionFrequencies`: the transition
/// frequencies `(Mf_low, Mf_high)` of the map for `m'`.
pub(crate) fn linear_frequency_map_transition_frequencies(
    emmprime: u32,
    mf_rd_22: f64,
    mf_rd_lm: f64,
    pnr_hm_mflow: f64,
    pnr_hm_mfhigh: f64,
) -> (f64, f64) {
    let mf_low = 0.65 * pnr_hm_mflow * emmprime as f64 / 2.0;
    let mut mf_high = 1.1 * (pnr_hm_mfhigh + (mf_rd_lm - mf_rd_22));
    // cases with the (2,1) mode where Mf_high might misbehave
    if mf_high < 0.0 || (mf_rd_lm - mf_rd_22 < 0.0 && mf_high < pnr_hm_mfhigh / 2.0) {
        mf_high = pnr_hm_mfhigh;
    }
    (mf_low, mf_high)
}

/// `IMRPhenomX_PNR_alphadot_cosbeta` integrated over `[t1, t2]` with Boole's
/// rule as in `IMRPhenomX_PNR_GeneratePNRGamma_FromInterpolants`, for the
/// integrand `f(t) = -alpha'(t) cos(beta(t))`.
#[inline]
pub(crate) fn boole_step(gamma_prev: f64, t1: f64, t2: f64, f: impl Fn(f64) -> f64) -> f64 {
    let sum = fma(
        7.0,
        f(t2),
        fma(
            32.0,
            f(fma(3.0, t1, t2) / 4.0),
            fma(
                12.0,
                f((t1 + t2) / 2.0),
                fma(7.0, f(t1), 32.0 * f(fma(3.0, t2, t1) / 4.0)),
            ),
        ),
    );
    fma(1.0 / 90.0 * (t2 - t1), sum, gamma_prev)
}

/// `-alpha'(f) cos(beta(f))`.
#[inline]
pub(crate) fn alphadot_cosbeta(alphadot: f64, beta: f64) -> f64 {
    -(alphadot * beta.cos())
}

/// `sqrt(1 - x^2)` as contracted by clang.
#[inline]
pub(crate) fn sqrt_one_minus_square(x: f64) -> f64 {
    fnms(x, x, 1.0).sqrt()
}

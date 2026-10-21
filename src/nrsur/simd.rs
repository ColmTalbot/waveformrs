// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2024-2026 Colm Talbot. Parts of waveformrs are translated
// from LALSuite, GSL and fdlibm: see NOTICE.md for their authors.
//! Four `f64` lanes with element-wise, separately rounded arithmetic (NEON
//! on aarch64, plain arrays elsewhere), used to evaluate independent time
//! samples together. Every operation rounds exactly as the scalar operation
//! does; nothing is fused.

use std::ops::{Add, Div, Mul, Neg, Sub};

#[cfg(target_arch = "aarch64")]
mod pair {
    use std::arch::aarch64::{
        float64x2_t, vaddq_f64, vdivq_f64, vdupq_n_f64, vld1q_f64, vmulq_f64, vnegq_f64, vst1q_f64,
        vsubq_f64,
    };

    #[derive(Clone, Copy)]
    pub(crate) struct F64x2(float64x2_t);

    // SAFETY (all): NEON is part of the aarch64 baseline, and the loads and
    // stores are from arrays of two elements.
    #[allow(unused_unsafe)]
    impl F64x2 {
        #[inline(always)]
        pub(crate) fn splat(x: f64) -> Self {
            Self(unsafe { vdupq_n_f64(x) })
        }
        #[inline(always)]
        pub(crate) fn from_array(a: [f64; 2]) -> Self {
            Self(unsafe { vld1q_f64(a.as_ptr()) })
        }
        #[inline(always)]
        pub(crate) fn to_array(self) -> [f64; 2] {
            let mut out = [0.0; 2];
            unsafe { vst1q_f64(out.as_mut_ptr(), self.0) };
            out
        }
        #[inline(always)]
        pub(crate) fn add(self, b: Self) -> Self {
            Self(unsafe { vaddq_f64(self.0, b.0) })
        }
        #[inline(always)]
        pub(crate) fn sub(self, b: Self) -> Self {
            Self(unsafe { vsubq_f64(self.0, b.0) })
        }
        #[inline(always)]
        pub(crate) fn mul(self, b: Self) -> Self {
            Self(unsafe { vmulq_f64(self.0, b.0) })
        }
        #[inline(always)]
        pub(crate) fn div(self, b: Self) -> Self {
            Self(unsafe { vdivq_f64(self.0, b.0) })
        }
        #[inline(always)]
        pub(crate) fn neg(self) -> Self {
            Self(unsafe { vnegq_f64(self.0) })
        }
    }
}

#[cfg(not(target_arch = "aarch64"))]
mod pair {
    #[derive(Clone, Copy)]
    pub(crate) struct F64x2([f64; 2]);

    impl F64x2 {
        #[inline(always)]
        pub(crate) fn splat(x: f64) -> Self {
            Self([x; 2])
        }
        #[inline(always)]
        pub(crate) fn from_array(a: [f64; 2]) -> Self {
            Self(a)
        }
        #[inline(always)]
        pub(crate) fn to_array(self) -> [f64; 2] {
            self.0
        }
        #[inline(always)]
        pub(crate) fn add(self, b: Self) -> Self {
            Self([self.0[0] + b.0[0], self.0[1] + b.0[1]])
        }
        #[inline(always)]
        pub(crate) fn sub(self, b: Self) -> Self {
            Self([self.0[0] - b.0[0], self.0[1] - b.0[1]])
        }
        #[inline(always)]
        pub(crate) fn mul(self, b: Self) -> Self {
            Self([self.0[0] * b.0[0], self.0[1] * b.0[1]])
        }
        #[inline(always)]
        pub(crate) fn div(self, b: Self) -> Self {
            Self([self.0[0] / b.0[0], self.0[1] / b.0[1]])
        }
        #[inline(always)]
        pub(crate) fn neg(self) -> Self {
            Self([-self.0[0], -self.0[1]])
        }
    }
}

use pair::F64x2;

/// Four `f64` lanes.
#[derive(Clone, Copy)]
pub(crate) struct F64x4([F64x2; 2]);

impl F64x4 {
    pub(crate) const LANES: usize = 4;

    #[inline(always)]
    pub(crate) fn splat(x: f64) -> Self {
        Self([F64x2::splat(x); 2])
    }

    #[inline(always)]
    pub(crate) fn from_array(a: [f64; 4]) -> Self {
        Self([
            F64x2::from_array([a[0], a[1]]),
            F64x2::from_array([a[2], a[3]]),
        ])
    }

    /// The first four elements of `s`.
    #[inline(always)]
    pub(crate) fn from_slice(s: &[f64]) -> Self {
        Self::from_array([s[0], s[1], s[2], s[3]])
    }

    #[inline(always)]
    pub(crate) fn to_array(self) -> [f64; 4] {
        let [a, b] = [self.0[0].to_array(), self.0[1].to_array()];
        [a[0], a[1], b[0], b[1]]
    }
}

macro_rules! binary {
    ($trait:ident, $method:ident) => {
        impl $trait for F64x4 {
            type Output = Self;
            #[inline(always)]
            fn $method(self, b: Self) -> Self {
                Self([self.0[0].$method(b.0[0]), self.0[1].$method(b.0[1])])
            }
        }
    };
}

binary!(Add, add);
binary!(Sub, sub);
binary!(Mul, mul);
binary!(Div, div);

impl Neg for F64x4 {
    type Output = Self;
    #[inline(always)]
    fn neg(self) -> Self {
        Self([self.0[0].neg(), self.0[1].neg()])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lane_wise() {
        let a = F64x4::from_array([1.0, -2.5, 3.25, 1e-300]);
        let b = F64x4::from_array([0.3, 7.0, -0.1, 1e10]);
        let x = [1.0, -2.5, 3.25, 1e-300];
        let y = [0.3, 7.0, -0.1, 1e10];
        let r = ((a * b + a) / b - a).to_array();
        for l in 0..4 {
            assert_eq!(r[l], (x[l] * y[l] + x[l]) / y[l] - x[l]);
        }
        assert_eq!((-a).to_array()[1], 2.5);
    }
}

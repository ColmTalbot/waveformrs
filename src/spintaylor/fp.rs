//! Floating point helpers that reproduce the arithmetic of the C code.
//!
//! LALSuite and GSL are plain C, and C compilers are allowed to contract an
//! expression `a * b + c` into a fused multiply-add (a single rounding).
//! Clang does so by default (`-ffp-contract=on`) on targets with hardware FMA,
//! which includes the arm64 (Apple silicon) builds of `lalsimulation` and
//! `gsl`.  The adaptive step size control of the integrator amplifies
//! last-bit differences in the derivatives into differences in the accepted
//! step sizes, so to agree with LAL to rounding precision we mirror clang's
//! contraction rule explicitly:
//!
//! For each `+`/`-` (including `+=`/`-=`) whose left operand is a product
//! computed in the same expression, the product is fused with the other
//! operand; otherwise, if the right operand is such a product, that one is
//! fused.  Products of compile-time constants are folded first and never
//! fused.
//!
//! With this the output of [`super::SpinTaylorT4::evolve_orbit`] is bit for
//! bit identical to the conda-forge osx-arm64 build of lalsimulation 6.2.0
//! (verified with the disassembly of `liblalsimulation`/`libgsl` and on
//! hundreds of evolutions).  Without it, differences of order 1e-7 (relative
//! to the largest value of a series) appear in the spin components.
//!
//! The generic x86-64 builds of LAL (GCC without `-mfma`) cannot contain
//! fused multiply-adds, so [`fma`] falls back to separately rounded
//! operations there (this also avoids the slow software `fma` of targets
//! without hardware support).  Other targets/compilers may contract
//! differently; expect agreement at the ~1e-7 level there.

/// Whether `a * b + c` is evaluated with a single rounding, mirroring how the
/// LAL/GSL binaries for this target are compiled.
pub const FUSED_MULTIPLY_ADD: bool = cfg!(target_arch = "aarch64");

/// `a * b + c`, fused if and only if the target LAL build fuses it.
#[inline(always)]
pub(crate) fn fma(a: f64, b: f64, c: f64) -> f64 {
    if FUSED_MULTIPLY_ADD {
        a.mul_add(b, c)
    } else {
        a * b + c
    }
}

/// `c - a * b` as compiled from the C expression `c - a * b`.
#[inline(always)]
pub(crate) fn fnms(a: f64, b: f64, c: f64) -> f64 {
    if FUSED_MULTIPLY_ADD {
        (-a).mul_add(b, c)
    } else {
        c - a * b
    }
}

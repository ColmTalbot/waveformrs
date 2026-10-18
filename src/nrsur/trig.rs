//! `sin` and `cos` as LALSimulation evaluates them.
//!
//! On Darwin, clang combines `sin(x)` and `cos(x)` of the same argument in
//! one basic block into a call to `__sincos_stret`, whose `sin` differs from
//! the C library's `sin` in the last bit for about 2 in 10^5 arguments.
//! The disassembly of the conda-forge osx-arm64 `liblalsimulation` shows
//! `__sincos_stret` in `PrecessingNRSur_ds_fit_x`/`assemble_dydt` (inlined
//! into `get_time_deriv_from_index` and `get_t_ref`) and in
//! `PrecessingNRSur_rotate_spins`, and separate `cos` and `sin` calls in
//! `TransformModesCoorbitalToInertial`. LLVM may do the same combination of
//! Rust's `f64::sin` and `f64::cos`, so both variants are called explicitly
//! here.

#[cfg(target_vendor = "apple")]
mod imp {
    #[repr(C)]
    struct SinCos {
        sin: f64,
        cos: f64,
    }

    extern "C" {
        fn __sincos_stret(x: f64) -> SinCos;
        fn sin(x: f64) -> f64;
        fn cos(x: f64) -> f64;
    }

    #[inline]
    pub(crate) fn sincos(x: f64) -> (f64, f64) {
        // SAFETY: pure functions of their argument.
        let sc = unsafe { __sincos_stret(x) };
        (sc.sin, sc.cos)
    }

    #[inline]
    pub(crate) fn separate_sin(x: f64) -> f64 {
        // SAFETY: as above. Declared without `readnone`, LLVM does not
        // combine these calls.
        unsafe { sin(x) }
    }

    #[inline]
    pub(crate) fn separate_cos(x: f64) -> f64 {
        // SAFETY: as above.
        unsafe { cos(x) }
    }
}

#[cfg(not(target_vendor = "apple"))]
mod imp {
    #[inline]
    pub(crate) fn sincos(x: f64) -> (f64, f64) {
        (x.sin(), x.cos())
    }

    #[inline]
    pub(crate) fn separate_sin(x: f64) -> f64 {
        x.sin()
    }

    #[inline]
    pub(crate) fn separate_cos(x: f64) -> f64 {
        x.cos()
    }
}

/// `(sin(x), cos(x))` where clang combines the two calls.
pub(crate) use imp::sincos;
/// `sin(x)` and `cos(x)` where they are separate calls in LAL.
pub(crate) use imp::{separate_cos, separate_sin};

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn close_to_std() {
        for x in [0.0, 0.3, -1.7, 123.456, -9876.5] {
            let (s, c) = sincos(x);
            assert!((s - x.sin()).abs() <= f64::EPSILON);
            assert!((c - x.cos()).abs() <= f64::EPSILON);
            assert!((separate_sin(x) - x.sin()).abs() <= f64::EPSILON);
            assert!((separate_cos(x) - x.cos()).abs() <= f64::EPSILON);
        }
    }
}

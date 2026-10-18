//! The C library's `cbrt`, which `f64::cbrt` is not on every platform.
//!
//! On linux-gnu, Rust links its `compiler-builtins` ahead of the system libm,
//! and that crate defines (weak) software versions of the exactly rounded
//! functions (`fma`, `floor`, `sqrt`, ...) and of `cbrt`. Its `cbrt` is
//! correctly rounded; glibc's is not, and differs from it in the last bit
//! for over half of all arguments. lalsimulation calls glibc's `cbrt` for
//! every power of the frequency, the reference velocity and the SpinTaylor
//! velocity, and the ill-conditioned parts of the models (the IMRPhenomXHM
//! (3,2) mode, the SpinTaylor step size control) turn these differences into
//! relative errors of up to 1e-3. [`LibmCbrt::libm_cbrt`] therefore calls
//! glibc's `cbrt` through its versioned symbol on linux-gnu (aarch64 and
//! x86-64), which the linker resolves in libm whatever the link order;
//! elsewhere it is `f64::cbrt` (the C library's on macOS).
//!
//! With this the IMRPhenomXHM modes and the SpinTaylor angles agree with the
//! conda-forge linux-64 build of lalsimulation bit for bit (docs/linux.md).

/// `cbrt` as the C library computes it.
pub(crate) trait LibmCbrt {
    fn libm_cbrt(self) -> f64;
}

#[cfg(all(
    target_os = "linux",
    target_env = "gnu",
    any(target_arch = "aarch64", target_arch = "x86_64")
))]
mod imp {
    extern "C" {
        // the first (and only) version of cbrt in glibc for the architecture
        #[cfg_attr(target_arch = "aarch64", link_name = "cbrt@GLIBC_2.17")]
        #[cfg_attr(target_arch = "x86_64", link_name = "cbrt@GLIBC_2.2.5")]
        fn cbrt(x: f64) -> f64;
    }

    #[inline]
    pub(super) fn libm_cbrt(x: f64) -> f64 {
        // SAFETY: a pure function of its argument.
        unsafe { cbrt(x) }
    }
}

#[cfg(not(all(
    target_os = "linux",
    target_env = "gnu",
    any(target_arch = "aarch64", target_arch = "x86_64")
)))]
mod imp {
    #[inline]
    pub(super) fn libm_cbrt(x: f64) -> f64 {
        x.cbrt()
    }
}

impl LibmCbrt for f64 {
    #[inline]
    fn libm_cbrt(self) -> f64 {
        imp::libm_cbrt(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn close_to_the_correctly_rounded_cube_root() {
        for x in [
            1e-300,
            2.5e-7,
            0.3,
            1.0,
            8.0,
            27.000_000_000_000_004,
            3.7e12,
            1e300,
        ] {
            let y = x.libm_cbrt();
            assert!((y - x.cbrt()).abs() <= f64::EPSILON * y, "{x}: {y}");
            assert_eq!((-x).libm_cbrt(), -y);
        }
        assert_eq!(0.0f64.libm_cbrt(), 0.0);
        assert!(f64::NAN.libm_cbrt().is_nan());
        assert_eq!(f64::INFINITY.libm_cbrt(), f64::INFINITY);
    }
}

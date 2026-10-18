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
//! relative errors of up to 1e-3.
//!
//! [`LibmCbrt::libm_cbrt`] is therefore glibc's algorithm on linux-gnu
//! (aarch64 and x86-64), [`glibc_cbrt`]: a port of
//! `sysdeps/ieee754/dbl-64/s_cbrt.c` with the fused multiply-adds of glibc's
//! aarch64 build, which the tests check bit for bit against the C library's
//! own `cbrt` (called through its versioned symbol, which the linker resolves
//! in libm whatever the link order). Inlined, it is also about twice as fast
//! as the call. Elsewhere it is `f64::cbrt`, which on macOS is the C
//! library's.
//!
//! With this the IMRPhenomXHM modes and the SpinTaylor angles agree with the
//! conda-forge linux-64 build of lalsimulation bit for bit (docs/linux.md).

/// `cbrt` as the C library computes it.
pub(crate) trait LibmCbrt {
    fn libm_cbrt(self) -> f64;
}

/// Whether the target is one whose glibc `cbrt` [`glibc_cbrt`] reproduces.
const GLIBC: bool = cfg!(all(
    target_os = "linux",
    target_env = "gnu",
    any(target_arch = "aarch64", target_arch = "x86_64")
));

impl LibmCbrt for f64 {
    #[inline(always)]
    fn libm_cbrt(self) -> f64 {
        if GLIBC {
            glibc_cbrt::<{ cfg!(target_arch = "aarch64") }>(self)
        } else {
            self.cbrt()
        }
    }
}

/// glibc's `cbrt` (`sysdeps/ieee754/dbl-64/s_cbrt.c`): `|x| = xm 2^xe` with
/// `xm` in [1/2, 1), a degree-6 polynomial for the cube root of `xm`, one
/// Halley step, the factor `2^((xe mod 3) / 3)` and `ldexp` by `xe / 3`. The
/// polynomial is evaluated with fused multiply-adds if `FUSED` (GCC contracts
/// it on aarch64, and cannot on baseline x86-64).
#[allow(clippy::excessive_precision)]
#[inline(always)]
fn glibc_cbrt<const FUSED: bool>(x: f64) -> f64 {
    const CBRT2: f64 = 1.259_921_049_894_873_164_8;
    const SQR_CBRT2: f64 = 1.587_401_051_968_199_474_8;
    const FACTOR: [f64; 5] = [1.0 / SQR_CBRT2, 1.0 / CBRT2, 1.0, CBRT2, SQR_CBRT2];
    let mul_add = |a: f64, b: f64, c: f64| if FUSED { a.mul_add(b, c) } else { a * b + c };
    let ax = x.abs();
    if !ax.is_finite() || ax == 0.0 {
        return x + x;
    }
    // frexp, with subnormals scaled by 2^54 first
    let (bits, adjust) = if ax < f64::MIN_POSITIVE {
        ((ax * f64::from_bits(0x4350_0000_0000_0000)).to_bits(), -54)
    } else {
        (ax.to_bits(), 0)
    };
    let xe = ((bits >> 52) & 0x7ff) as i32 - 1022 + adjust;
    let xm = f64::from_bits((bits & !(0x7ff << 52)) | (1022 << 52));
    let u = mul_add(-0.145_263_899_385_486_377, xm, 0.784_932_344_976_639_262);
    let u = mul_add(u, xm, -1.834_692_774_836_130_86);
    let u = mul_add(u, xm, 2.446_931_225_635_344_30);
    let u = mul_add(u, xm, -2.114_994_941_673_712_87);
    let u = mul_add(u, xm, 1.508_191_937_815_848_96);
    let u = mul_add(u, xm, 0.354_895_765_043_919_860);
    let t2 = u * u * u;
    // (the products by 2 are exact: fusing them changes nothing)
    let ym = u * (t2 + 2.0 * xm) / (2.0 * t2 + xm) * FACTOR[(2 + xe % 3) as usize];
    // ldexp: one (correctly rounded) multiplication by 2^(xe / 3), which is
    // normal since |xe / 3| <= 358
    let scale = f64::from_bits(((xe / 3 + 1023) as u64) << 52);
    (if x > 0.0 { ym } else { -ym }) * scale
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
            4.9e-324,
            2.2e-308,
        ] {
            for y in [x.libm_cbrt(), glibc_cbrt::<true>(x), glibc_cbrt::<false>(x)] {
                assert!((y - x.cbrt()).abs() <= 2.0 * f64::EPSILON * y, "{x}: {y}");
            }
            assert_eq!((-x).libm_cbrt(), -x.libm_cbrt());
        }
        for x in [0.0, -0.0, f64::INFINITY, f64::NEG_INFINITY] {
            assert_eq!(x.libm_cbrt().to_bits(), x.to_bits());
            assert_eq!(glibc_cbrt::<false>(x).to_bits(), x.to_bits());
        }
        assert!(f64::NAN.libm_cbrt().is_nan());
        assert!(glibc_cbrt::<true>(f64::NAN).is_nan());
    }

    /// The port against the C library's `cbrt` on linux-gnu: random bit
    /// patterns (every exponent, subnormals included) and a dense range.
    #[cfg(all(
        target_os = "linux",
        target_env = "gnu",
        any(target_arch = "aarch64", target_arch = "x86_64")
    ))]
    #[test]
    fn matches_glibc_bit_for_bit() {
        extern "C" {
            // the first (and only) version of cbrt in glibc for the architecture
            #[cfg_attr(target_arch = "aarch64", link_name = "cbrt@GLIBC_2.17")]
            #[cfg_attr(target_arch = "x86_64", link_name = "cbrt@GLIBC_2.2.5")]
            fn cbrt(x: f64) -> f64;
        }
        let mut state = 0x9E37_79B9_7F4A_7C15u64;
        let mut differ = Vec::new();
        for i in 0..2_000_000u32 {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            let x = if i % 2 == 0 {
                f64::from_bits(state)
            } else {
                1e-4 + 1e-6 * f64::from(i)
            };
            // SAFETY: a pure function of its argument.
            let reference = unsafe { cbrt(x) };
            let port = x.libm_cbrt();
            if port.to_bits() != reference.to_bits() && !(x.is_nan() && port.is_nan()) {
                differ.push(x);
            }
        }
        assert!(
            differ.is_empty(),
            "{} arguments differ, e.g. {:?}",
            differ.len(),
            &differ[..differ.len().min(5)]
        );
    }
}

// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2024-2026 Colm Talbot. Parts of waveformrs are translated
// from LALSuite, GSL, glibc and fdlibm: see NOTICE.md for their authors.
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
//! # Choice per platform
//!
//! [`FUSED_MULTIPLY_ADD`] depends on the architecture only, and was chosen
//! by comparing with the conda-forge builds of lalsimulation 6.2.0 on each
//! platform (docs/linux.md, report/data/linux_validation_*.json):
//!
//! * aarch64, macOS: fused. Clang (`-ffp-contract=on`) builds LAL and GSL,
//!   OpenBLAS's kernels fuse, and with the same libm the port is bit for
//!   bit identical.
//! * aarch64, Linux: fused. GCC 14 builds LAL and GSL with its default
//!   `-ffp-contract=fast`, which fuses every product that clang fuses within
//!   an expression (with the same choice of operand in `a * b + c * d`) and,
//!   in addition, products and sums in different statements (about 11% more
//!   fused instructions in `liblalsimulation.so` than in the osx-arm64
//!   dylib). Mirroring clang's subset is not bit for bit, but it is closer
//!   to the Linux build than not fusing at all (over random draws, median
//!   IMRPhenomXHM errors 2e-7 against 5e-7, NRSur7dq4 1e-14 against 1e-13),
//!   and OpenBLAS's neoversen1 kernels fuse exactly as mirrored here (the
//!   GSL LU solves agree bit for bit). glibc's libm, not Apple's, makes up
//!   the rest of the difference.
//! * x86-64 (and other architectures): separately rounded. The conda-forge
//!   linux-64 builds use `-march=nocona`, which has no FMA instruction:
//!   `liblalsimulation.so`, `liblal.so` and `libgsl.so` contain none, and the
//!   remnant and PN fits agree with them bit for bit. This also avoids the
//!   software `fma` call that `f64::mul_add` becomes without the `fma` target
//!   feature. OpenBLAS selects its kernels at run time; on CPUs with FMA
//!   (Haswell and later) they fuse, so the LU solves of
//!   `imrphenomx::linalg` then agree to rounding only.
//!
//! `--cfg waveformrs_fma="on"` or `"off"` (in `RUSTFLAGS`) overrides the
//! choice, e.g. for a LAL built with clang on Linux or with `-march=haswell`.

/// Whether `a * b + c` is evaluated with a single rounding, mirroring how the
/// LAL/GSL binaries for this target are compiled (see the module docs).
pub const FUSED_MULTIPLY_ADD: bool = cfg!(any(
    waveformrs_fma = "on",
    all(not(waveformrs_fma = "off"), target_arch = "aarch64")
));

/// `a * b + c`, fused if and only if the target LAL build fuses it.
#[inline(always)]
pub(crate) fn fma(a: f64, b: f64, c: f64) -> f64 {
    if FUSED_MULTIPLY_ADD {
        a.mul_add(b, c)
    } else {
        a * b + c
    }
}

/// `a * b - c` as compiled from the C expression `a * b - c`.
#[inline(always)]
pub(crate) fn fms(a: f64, b: f64, c: f64) -> f64 {
    if FUSED_MULTIPLY_ADD {
        a.mul_add(b, -c)
    } else {
        a * b - c
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

/// Dot product as written in LAL, `ax*bx + ay*by + az*bz`.
#[inline(always)]
pub(crate) fn dot(a: &[f64; 3], b: &[f64; 3]) -> f64 {
    fma(a[2], b[2], fma(a[0], b[0], a[1] * b[1]))
}

/// Cross product as written in `XLALSimInspiralVectorCrossProduct`.
#[inline(always)]
pub(crate) fn cross(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    [
        fms(a[1], b[2], a[2] * b[1]),
        fms(a[2], b[0], a[0] * b[2]),
        fms(a[0], b[1], a[1] * b[0]),
    ]
}

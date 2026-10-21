// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2024-2026 Colm Talbot. Parts of waveformrs are translated
// from LALSuite, GSL, glibc and fdlibm: see NOTICE.md for their authors.
//! The IMRPhenomX family of frequency-domain phenomenological models.
//!
//! This is a port of the IMRPhenomX implementation in LALSimulation
//! (v6.2.0) for the default model versions. The module layout mirrors the
//! LAL sources so that higher-mode and precessing extensions can follow the
//! C code:
//!
//! | module         | LALSimulation source                                   |
//! |----------------|--------------------------------------------------------|
//! | `xas`          | `LALSimIMRPhenomX.c` (IMRPhenomXAS driver)             |
//! | `internals`    | `LALSimIMRPhenomX_internals.c` (waveform struct, powers) |
//! | `amplitude`    | `LALSimIMRPhenomX_internals.c` (22 amplitude)          |
//! | `phase`        | `LALSimIMRPhenomX_internals.c` (22 phase)              |
//! | `inspiral`     | `LALSimIMRPhenomX_inspiral.c`                          |
//! | `intermediate` | `LALSimIMRPhenomX_intermediate.c`                      |
//! | `ringdown`     | `LALSimIMRPhenomX_ringdown.c`                          |
//! | `qnm`          | `LALSimIMRPhenomX_qnm.c`                               |
//! | `utilities`    | `LALSimIMRPhenomXUtilities.c`                          |
//! | `linalg`       | the GSL LU solves used for the collocation systems     |
//! | `xp`           | `LALSimIMRPhenomX.c` (IMRPhenomXP driver)              |
//! | `precessing_grid` | IMRPhenomXP/XPHM on a prepared frequency grid (no LAL counterpart) |
//! | `precession`   | `LALSimIMRPhenomX_precession.c` (precession struct, J frame, NNLO angles) |
//! | `spintaylor_angles` | `LALSimIMRPhenomX_precession.c` (SpinTaylor Euler angles) |
//! | `msa`          | `LALSimIMRPhenomX_precession.c` (MSA Euler angles)     |
//! | `wigner`       | Wigner-d matrices and twisting up, spherical harmonics |
//! | `xhm`          | `LALSimIMRPhenomXHM.c` (IMRPhenomXHM driver, mode sum) |
//! | `xphm`         | `LALSimIMRPhenomXPHM.c` (IMRPhenomXPHM driver, twisting up) |
//! | `xpnr`         | the IMRPhenomXPNR flags of `LALSimInspiralGeneratorLegacy.c` on top of `LALSimIMRPhenomXPHM.c` |
//! | `pnr`          | `LALSimIMRPhenomX_PNR*.c`, `LALSimIMRPhenomX_AntisymmetricWaveform.c` (PNR set-up and angles) |
//! | `pnr_angles`   | `LALSimIMRPhenomX_PNR_alpha.c`, `_beta.c`, `_internals.c` (model-independent PNR helpers) |
//! | `pnr_fits`     | `LALSimIMRPhenomX_PNR_coefficients.c`, `_deviations.c`  |
//! | `hm_internals` | `LALSimIMRPhenomXHM_internals.c` (mode struct, coefficients) |
//! | `hm_inspiral`  | `LALSimIMRPhenomXHM_inspiral.c`                        |
//! | `hm_intermediate` | `LALSimIMRPhenomXHM_intermediate.c`                 |
//! | `hm_ringdown`  | `LALSimIMRPhenomXHM_ringdown.c`                        |
//! | `hm_qnm`       | `LALSimIMRPhenomXHM_qnm.c`                             |
//! | `hm_amp_122019` | amplitude of the original IMRPhenomXHM release (122019) |
//! | `hm_fits_122018` | parameter-space fits of the 122019 amplitude        |
//! | `hm_poly_122019` | intermediate amplitude polynomial of the 122019 release |
//! | `poly_roots`   | GSL's `gsl_poly_complex_solve`                         |
//! | `libm`         | the C library's `cexp`                                 |

mod amplitude;
mod hm_amp_122019;
mod hm_fits_122018;
mod hm_inspiral;
mod hm_intermediate;
mod hm_internals;
mod hm_poly_122019;
mod hm_qnm;
mod hm_ringdown;
mod inspiral;
mod intermediate;
mod internals;
pub(crate) mod libm;
mod linalg;
mod msa;
mod nrtidal;
mod phase;
mod pnr;
mod pnr_angles;
mod pnr_fits;
mod poly_roots;
mod precessing_grid;
mod precession;
mod qnm;
mod ringdown;
mod spintaylor_angles;
mod utilities;
mod version_fallback;
mod wigner;
mod xas;
mod xhm;
mod xp;
mod xphm;
mod xpnr;

pub use internals::IMRPhenomXError;
pub use precession::{PrecessingBinary, PrecessionVersion};
pub use spintaylor_angles::EulerAngles;
pub use xas::IMRPhenomXAS;
pub use xhm::{IMRPhenomXHM, IMRPHENOMXHM_MODES};
pub use xp::{IMRPhenomXP, Polarizations};
pub use xphm::IMRPhenomXPHM;

pub(crate) use xhm::spin_weighted_spherical_harmonic;
pub use xpnr::IMRPhenomXPNR;

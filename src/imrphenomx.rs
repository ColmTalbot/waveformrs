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
//! | `precession`   | `LALSimIMRPhenomX_precession.c` (precession struct, J frame, NNLO angles) |
//! | `spintaylor_angles` | `LALSimIMRPhenomX_precession.c` (SpinTaylor Euler angles) |
//! | `msa`          | `LALSimIMRPhenomX_precession.c` (MSA Euler angles)     |
//! | `wigner`       | Wigner-d matrices and twisting up, spherical harmonics |
//! | `xhm`          | `LALSimIMRPhenomXHM.c` (IMRPhenomXHM driver, mode sum) |
//! | `xphm`         | `LALSimIMRPhenomXPHM.c` (IMRPhenomXPHM driver, twisting up) |
//! | `hm_internals` | `LALSimIMRPhenomXHM_internals.c` (mode struct, coefficients) |
//! | `hm_inspiral`  | `LALSimIMRPhenomXHM_inspiral.c`                        |
//! | `hm_intermediate` | `LALSimIMRPhenomXHM_intermediate.c`                 |
//! | `hm_ringdown`  | `LALSimIMRPhenomXHM_ringdown.c`                        |
//! | `hm_qnm`       | `LALSimIMRPhenomXHM_qnm.c`                             |
//! | `libm`         | the C library's `cexp`                                 |

mod amplitude;
mod hm_inspiral;
mod hm_intermediate;
mod hm_internals;
mod hm_qnm;
mod hm_ringdown;
mod inspiral;
mod intermediate;
mod internals;
mod libm;
mod linalg;
mod msa;
mod phase;
mod precession;
mod qnm;
mod ringdown;
mod spintaylor_angles;
mod utilities;
mod wigner;
mod xas;
mod xhm;
mod xp;
mod xphm;

pub use internals::IMRPhenomXError;
pub use precession::{PrecessingBinary, PrecessionVersion};
pub use spintaylor_angles::EulerAngles;
pub use xas::IMRPhenomXAS;
pub use xp::{IMRPhenomXP, Polarizations};
pub use xhm::{IMRPhenomXHM, IMRPHENOMXHM_MODES};
pub use xphm::IMRPhenomXPHM;

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
//! | `wigner`       | Wigner-d matrices and twisting up, spherical harmonics |

mod amplitude;
mod inspiral;
mod intermediate;
mod internals;
mod linalg;
mod phase;
mod precession;
mod qnm;
mod ringdown;
mod spintaylor_angles;
mod utilities;
mod wigner;
mod xas;
mod xp;

pub use internals::IMRPhenomXError;
pub use precession::{PrecessingBinary, PrecessionVersion};
pub use spintaylor_angles::EulerAngles;
pub use xas::IMRPhenomXAS;
pub use xp::{IMRPhenomXP, Polarizations};

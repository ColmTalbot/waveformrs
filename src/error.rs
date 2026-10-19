// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2024-2026 Colm Talbot. Parts of waveformrs are translated
// from LALSuite, GSL, glibc and fdlibm: see NOTICE.md for their authors.
//! The error type of the common model interface ([`crate::model`]) and the
//! batch evaluation ([`crate::batch`]).
//!
//! [`Error`] wraps the errors of the individual model families
//! ([`IMRPhenomXError`], [`NRSurError`]) and adds the failures of the
//! interface itself, such as in-plane spins passed to an aligned-spin
//! approximant (which `SimInspiralChooseFDWaveform` rejects too).

use std::fmt;

use crate::imrphenomx::IMRPhenomXError;
use crate::model::Approximant;
use crate::nrsur::NRSurError;

/// `Result` with [`Error`].
pub type Result<T, E = Error> = std::result::Result<T, E>;

/// Errors of the common model interface.
#[derive(Debug)]
#[non_exhaustive]
pub enum Error {
    /// An error of the IMRPhenomX family (parameters LALSimulation rejects,
    /// failed Euler angles, ...).
    IMRPhenomX(IMRPhenomXError),
    /// An error of NRSur7dq4 (data file, extrapolation limits, ...).
    NRSur(NRSurError),
    /// A parameter LALSimulation rejects for this approximant.
    InvalidParameter {
        approximant: Approximant,
        message: &'static str,
    },
    /// Non-zero in-plane spin components for an aligned-spin approximant
    /// (`checkTransverseSpinsZero` in `SimInspiralChooseFDWaveform`).
    InPlaneSpins(Approximant),
    /// The approximant does not support the request (a mode array for a
    /// model without one, `ChooseFDModes` for a model LAL does not support,
    /// a frequency-domain request for a time-domain model, ...).
    Unsupported {
        approximant: Approximant,
        what: &'static str,
    },
    /// The name is not one of [`Approximant::ALL`].
    UnknownApproximant(String),
    /// A request the interface cannot serve, such as building a boxed model
    /// with `from_binary`.
    InvalidRequest(&'static str),
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::IMRPhenomX(err) => write!(f, "{err}"),
            Self::NRSur(err) => write!(f, "{err}"),
            Self::InvalidParameter {
                approximant,
                message,
            } => write!(f, "{approximant}: {message}"),
            Self::InPlaneSpins(approximant) => write!(
                f,
                "{approximant}: non-zero transverse spins were given, but this is a \
                 non-precessing approximant"
            ),
            Self::Unsupported { approximant, what } => {
                write!(f, "{approximant} does not support {what}")
            }
            Self::UnknownApproximant(name) => {
                write!(f, "unknown approximant {name:?}; known: ")?;
                for (i, approximant) in Approximant::ALL.iter().enumerate() {
                    let separator = if i == 0 { "" } else { ", " };
                    write!(f, "{separator}{approximant}")?;
                }
                Ok(())
            }
            Self::InvalidRequest(why) => write!(f, "{why}"),
        }
    }
}

impl std::error::Error for Error {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::IMRPhenomX(err) => Some(err),
            Self::NRSur(err) => Some(err),
            _ => None,
        }
    }
}

impl From<IMRPhenomXError> for Error {
    fn from(err: IMRPhenomXError) -> Self {
        Self::IMRPhenomX(err)
    }
}

impl From<NRSurError> for Error {
    fn from(err: NRSurError) -> Self {
        Self::NRSur(err)
    }
}

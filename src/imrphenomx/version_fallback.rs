// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2024-2026 Colm Talbot. Parts of waveformrs are translated
// from LALSuite, GSL and fdlibm: see NOTICE.md for their authors.
//! The automatic changes of the Euler-angle prescription of IMRPhenomXP and
//! IMRPhenomXPHM, and the errors that replace them when they are switched
//! off (`with_version_fallback(false)`).
//!
//! Like LALSimulation, by default the models
//!
//! * use the NNLO angles (102) instead of the SpinTaylor angles (320) when
//!   the in-plane spins are below `1e-7`;
//! * use the MSA angles (223) when the SpinTaylor angles cannot be computed
//!   (the PN integration does not reach the first frequency);
//! * use the NNLO angles (102) when the MSA system cannot be initialised.
//!
//! Without fallbacks each of these is an
//! [`IMRPhenomXError::PrecessionVersionUnavailable`] instead.

use crate::imrphenomx::internals::IMRPhenomXError;
use crate::imrphenomx::precession::{PrecessionStruct, PrecessionVersion};

/// An error unless the precession struct uses the `requested` prescription
/// (it falls back to the NNLO angles for negligible in-plane spins with the
/// SpinTaylor angles and for a failed MSA system).
pub(crate) fn check_version(
    prec: &PrecessionStruct,
    requested: PrecessionVersion,
) -> Result<(), IMRPhenomXError> {
    let flag = requested.lal_flag();
    if prec.version == flag {
        return Ok(());
    }
    let reason = match requested {
        PrecessionVersion::SpinTaylor => "the in-plane spins are below 1e-7",
        PrecessionVersion::Msa => "the MSA system cannot be initialised",
        _ => "the prescription is not available for this binary",
    };
    Err(IMRPhenomXError::PrecessionVersionUnavailable {
        requested: flag,
        fallback: prec.version,
        reason,
    })
}

/// The failure of the SpinTaylor angles (which LAL answers with the MSA
/// angles) as an error.
pub(crate) fn forbid_spintaylor_fallback<T>(
    result: Result<T, IMRPhenomXError>,
) -> Result<T, IMRPhenomXError> {
    match result {
        Err(IMRPhenomXError::SpinTaylorFallbackToMsa(reason)) => {
            Err(IMRPhenomXError::PrecessionVersionUnavailable {
                requested: PrecessionVersion::SpinTaylor.lal_flag(),
                fallback: PrecessionVersion::Msa.lal_flag(),
                reason,
            })
        }
        result => result,
    }
}

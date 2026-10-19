// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2024-2026 Colm Talbot. Parts of waveformrs are translated
// from LALSuite, GSL, glibc and fdlibm: see NOTICE.md for their authors.
// Constants are quoted with the digits used in LALConstants.h.
#![allow(clippy::excessive_precision)]

pub const PC_SI: f64 = 3.085677581491367278913937957796471611e16;
pub const MPC_SI: f64 = 1.0e6 * PC_SI;
pub const SPEED_OF_LIGHT_IN_VACUUM: f64 = 299_792_458.0;
pub const GMSUN_SI: f64 = 1.3271244e20;
pub const SOLAR_RADIUS_IN_M: f64 = GMSUN_SI / SPEED_OF_LIGHT_IN_VACUUM / SPEED_OF_LIGHT_IN_VACUUM;
pub const SOLAR_RADIUS_IN_S: f64 = SOLAR_RADIUS_IN_M / SPEED_OF_LIGHT_IN_VACUUM;
pub const EGAMMA: f64 = 0.577215664901532860606512090082402431;

// Values as defined in LALConstants.h (lal 7.7), quoted to full double precision.
/// Solar mass in kg.
pub const MSUN_SI: f64 = 1.988409870698050731911960804878414216e30;
/// Geometrized solar mass in metres, G M_sun / c^2.
pub const MRSUN_SI: f64 = 1.476625038050124729627979840144936351e3;
/// Geometrized solar mass in seconds, G M_sun / c^3.
pub const MTSUN_SI: f64 = 4.925490947641266978197229498498379006e-6;
/// Newton's gravitational constant, N m^2 kg^-2.
pub const G_SI: f64 = 6.67430e-11;
/// Speed of light in vacuum, m s^-1.
pub const C_SI: f64 = SPEED_OF_LIGHT_IN_VACUUM;
/// Single precision machine epsilon, used for frequency comparisons in LAL.
pub const REAL4_EPS: f64 = 1.1920928955078125e-07;

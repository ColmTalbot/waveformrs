// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2024-2026 Colm Talbot. Parts of waveformrs are translated
// from LALSuite, GSL, glibc and fdlibm: see NOTICE.md for their authors.
use std::f64::consts::PI;

use once_cell::sync::Lazy;

use crate::imrphenomd::utils::Powers;

pub const PHI_FJOIN_INS: f64 = 0.018;
pub const AMP_FJOIN_INS: f64 = 0.014;
pub static POWERS_OF_PI: Lazy<Powers> = Lazy::new(|| Powers::new(PI));

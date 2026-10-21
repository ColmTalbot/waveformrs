// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2024-2026 Colm Talbot. Parts of waveformrs are translated
// from LALSuite, GSL and fdlibm: see NOTICE.md for their authors.
//! How the unit tests compare with lalsimulation's values.
//!
//! The values hard-coded in the tests and the data in `tests/data` come from
//! the conda-forge osx-arm64 build of lalsimulation 6.2.0. The crate does not
//! reproduce any build of LAL operation by operation: it uses Rust's
//! elementary functions and fuses multiply-adds where the hardware can, so
//! the tests compare to a few units in the last place.
//! `WAVEFORMRS_REFERENCE_DATA` points the tests that read reference data to a
//! directory generated with this machine's lalsimulation (see docs/linux.md).

use std::path::PathBuf;

/// The reference data `<name>.json` (`tests/data` or
/// `WAVEFORMRS_REFERENCE_DATA`).
pub(crate) fn data_path(name: &str) -> PathBuf {
    std::env::var_os("WAVEFORMRS_REFERENCE_DATA")
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data"))
        .join(format!("{name}.json"))
}

/// `value` against lalsimulation's `expected`, within `ulps` units of
/// `max(|expected|, 1)` in the last place.
#[track_caller]
pub(crate) fn assert_matches_lal(name: &str, value: f64, expected: f64, ulps: f64) {
    let tolerance = ulps * f64::EPSILON * expected.abs().max(1.0);
    assert!(
        (value - expected).abs() <= tolerance,
        "{name}: {value} differs from lalsimulation's {expected} by more than {ulps} ulp"
    );
}

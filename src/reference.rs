//! How the unit tests compare with lalsimulation's values.
//!
//! The values hard-coded in the tests and the data in `tests/data` come from
//! the conda-forge osx-arm64 build of lalsimulation 6.2.0, which this crate
//! reproduces bit for bit on macOS arm64: it mirrors clang's fused
//! multiply-adds (`spintaylor::fp`) and calls the same libm. On other
//! platforms the C library differs (glibc's `sin`, `pow`, `cbrt`, ... are not
//! Apple's), so the same tests compare to a few units in the last place.
//! `WAVEFORMRS_REFERENCE_DATA` points the tests that read reference data to a
//! directory generated with this machine's lalsimulation (see docs/linux.md).

use std::path::PathBuf;

/// Whether the hard-coded lalsimulation values can be reproduced to the last
/// bit (same fused multiply-adds, same libm).
pub(crate) const BIT_EXACT: bool = cfg!(all(target_arch = "aarch64", target_os = "macos"));

/// The reference data `<name>.json` (`tests/data` or
/// `WAVEFORMRS_REFERENCE_DATA`).
pub(crate) fn data_path(name: &str) -> PathBuf {
    std::env::var_os("WAVEFORMRS_REFERENCE_DATA")
        .map(PathBuf::from)
        .unwrap_or_else(|| PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/data"))
        .join(format!("{name}.json"))
}

/// Whether the reference data come from the lalsimulation of this platform.
pub(crate) fn same_platform() -> bool {
    BIT_EXACT || std::env::var_os("WAVEFORMRS_REFERENCE_DATA").is_some()
}

/// `value` against lalsimulation's `expected`: identical where [`BIT_EXACT`],
/// within `ulps` units of `max(|expected|, 1)` in the last place elsewhere.
#[track_caller]
pub(crate) fn assert_matches_lal(name: &str, value: f64, expected: f64, ulps: f64) {
    if BIT_EXACT {
        assert_eq!(value, expected, "{name} differs from lalsimulation");
    } else {
        let tolerance = ulps * f64::EPSILON * expected.abs().max(1.0);
        assert!(
            (value - expected).abs() <= tolerance,
            "{name}: {value} differs from lalsimulation's {expected} by more than {ulps} ulp"
        );
    }
}

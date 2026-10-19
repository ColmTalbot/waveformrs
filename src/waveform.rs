// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2024-2026 Colm Talbot. Parts of waveformrs are translated
// from LALSuite, GSL, glibc and fdlibm: see NOTICE.md for their authors.
use num_complex::Complex;

use crate::grid::FrequencyGrid;

/// The aligned-spin (2,2)-mode models TaylorF2, IMRPhenomD and
/// IMRPhenomXAS, evaluated as `h = A exp(-i Psi)`.
///
/// # Reference phase
///
/// In every method `phi_c` is the orbital reference phase, LALSimulation's
/// `phi_ref` (the argument `phiRef` of `SimInspiralChooseFDWaveform`), for
/// every model: the (2,2) mode is shifted by `2 phi_c`, so that
/// [`Waveform::waveform_modes`] with `phi_c = phi_ref` reproduces LAL's
/// polarizations. (IMRPhenomD and IMRPhenomXAS used to take the
/// gravitational-wave phase `2 phi_ref` here; see `CHANGELOG.md`.)
///
/// For an interface with all parameters, including `phi_ref` and the
/// inclination, fixed at construction, see
/// [`crate::model::FrequencyDomainModel`].
pub trait Waveform {
    /// The frequency variable `v` of [`Waveform::phase`] and
    /// [`Waveform::amplitude`] at `frequency` (Hz): the orbital speed
    /// `(pi M f)^(1/3)` for TaylorF2, the dimensionless frequency `Mf` for
    /// IMRPhenomD and IMRPhenomXAS.
    fn orbital_speed(&self, frequency: f64) -> f64;
    /// The phase `Psi(v)` for the orbital reference phase `phi_c`.
    fn phase(&self, v: f64, phi_c: f64) -> f64;
    /// The amplitude `A(v)`.
    fn amplitude(&self, v: f64) -> f64;

    /// The face-on `h_+` at `frequency` (Hz) for the orbital reference phase
    /// `phi_c`.
    fn waveform_single_frequency(&self, frequency: f64, phi_c: f64) -> Complex<f64> {
        let orbital_speed = self.orbital_speed(frequency);
        self.amplitude(orbital_speed)
            * Complex::<f64> {
                re: 0.0,
                im: -self.phase(orbital_speed, phi_c),
            }
            .exp()
    }
    /// The face-on `h_+` at `frequencies` (Hz) for the orbital reference
    /// phase `phi_c`.
    fn waveform(&self, frequencies: &[f64], phi_c: f64) -> Vec<Complex<f64>> {
        frequencies
            .iter()
            .map(|&f| self.waveform_single_frequency(f, phi_c))
            .collect()
    }

    /// The polarizations `(h_+, h_x)` at `frequencies` (Hz) for the orbital
    /// reference phase `phi_c` and the angle `theta_jn` between the orbital
    /// angular momentum and the line of sight: `h_+ = (1 + cos^2) / 2 h`,
    /// `h_x = -i cos h` with `h` from [`Waveform::waveform`], as
    /// `SimInspiralChooseFDWaveform` builds them.
    fn waveform_modes(
        &self,
        frequencies: &[f64],
        phi_c: f64,
        theta_jn: f64,
    ) -> (Vec<Complex<f64>>, Vec<Complex<f64>>) {
        let mut plus = self.waveform(frequencies, phi_c);
        let cos_theta_jn = theta_jn.cos();
        let plus_factor = (1.0 + cos_theta_jn * cos_theta_jn) / 2.0;
        // multiplying by -i cos(theta_jn) without a full complex product
        let cross = plus
            .iter()
            .map(|w| Complex::new(w.im * cos_theta_jn, -w.re * cos_theta_jn))
            .collect();
        plus.iter_mut().for_each(|w| *w *= plus_factor);
        (plus, cross)
    }

    /// [`Waveform::waveform`] on a prepared [`FrequencyGrid`], written into
    /// `out` (which must have the length of the grid).
    ///
    /// For evaluating many waveforms on the same frequencies: the models of
    /// this crate reuse the powers of the frequencies stored in the grid and
    /// allocate nothing. The result agrees with [`Waveform::waveform`] to
    /// rounding.
    fn waveform_on_grid(&self, grid: &FrequencyGrid, phi_c: f64, out: &mut [Complex<f64>]) {
        out.copy_from_slice(&self.waveform(grid.frequencies(), phi_c));
    }

    /// [`Waveform::waveform_modes`] on a prepared [`FrequencyGrid`], written
    /// into `plus` and `cross` (each with the length of the grid).
    fn waveform_modes_on_grid(
        &self,
        grid: &FrequencyGrid,
        phi_c: f64,
        theta_jn: f64,
        plus: &mut [Complex<f64>],
        cross: &mut [Complex<f64>],
    ) {
        self.waveform_on_grid(grid, phi_c, plus);
        let cos_theta_jn = theta_jn.cos();
        let plus_factor = (1.0 + cos_theta_jn * cos_theta_jn) / 2.0;
        for (p, c) in plus.iter_mut().zip(cross.iter_mut()) {
            *c = Complex::new(p.im * cos_theta_jn, -p.re * cos_theta_jn);
            *p *= plus_factor;
        }
    }
}

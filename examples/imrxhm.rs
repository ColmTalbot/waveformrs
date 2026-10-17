//! Timing of IMRPhenomXHM on a 0.25 Hz grid from 20 to 1500 Hz for
//! M = 90 Msun, q = 0.5 (as `examples/imrxas.rs`).
//!
//! Compare with lalsimulation (single thread):
//!
//! ```python
//! import lal, lalsimulation as ls, numpy as np, timeit
//! p = lal.CreateDict(); ls.SimInspiralWaveformParamsInsertPhenomXHMThresholdMband(p, 0)
//! f = lal.CreateREAL8Vector(5921); f.data = 20 + 0.25 * np.arange(5921)
//! args = (0.1, 60 * lal.MSUN_SI, 30 * lal.MSUN_SI, 0, 0, 0.0, 0, 0, 0.0, 20.0,
//!         1e8 * lal.PC_SI, 0.7, p, ls.IMRPhenomXHM, f)
//! timeit.timeit(lambda: ls.SimInspiralChooseFDWaveformSequence(*args), number=100) / 100
//! ```
//!
//! `SimInspiralChooseFDWaveformSequence` evaluates the higher modes on the
//! whole grid, as `IMRPhenomXHM::polarizations` does; on a uniform grid
//! `SimInspiralChooseFDWaveform` stops at Mf = 0.3 (and may multiband), so
//! compare it with the polarizations below that frequency. The single modes
//! correspond to `SimIMRPhenomXHMFrequencySequenceOneMode`.
use std::time::Instant;

use waveformrs::imrphenomx::IMRPhenomXHM;

fn time<F: FnMut()>(label: &str, iterations: u32, mut f: F) {
    f();
    let start = Instant::now();
    for _ in 0..iterations {
        f();
    }
    println!("{label:<45} {:>10.1?}", start.elapsed() / iterations);
}

fn main() {
    let total_mass = 90.0;
    let mass_ratio = 0.5;
    let (chi_1, chi_2) = (0.0, 0.0);
    let luminosity_distance = 100.0;
    let (phi_ref, inclination) = (0.1, 0.7);

    let delta_f = 0.25;
    let (f_lower, f_max) = (20.0, 1500.0);
    let n = ((f_max - f_lower) / delta_f) as usize + 1;
    let frequencies: Vec<f64> = (0..n).map(|i| f_lower + i as f64 * delta_f).collect();
    println!("{n} frequencies, M = {total_mass} Msun, q = {mass_ratio}");

    // `new` only stores the parameters: each call sets up the coefficients.
    let model = IMRPhenomXHM::new(total_mass, mass_ratio, chi_1, chi_2, luminosity_distance)
        .with_reference_frequency(f_lower);
    time("polarizations", 1000, || {
        std::hint::black_box(model.polarizations(&frequencies, phi_ref, inclination));
    });
    // LAL's uniform-grid path (ChooseFDWaveform) stops at Mf = 0.3.
    let f_cut = 0.3 / (total_mass * 4.925490947641267e-6);
    let below_cut: Vec<f64> = frequencies
        .iter()
        .copied()
        .filter(|&f| f <= f_cut)
        .collect();
    time(&format!("polarizations, f <= {f_cut:.0} Hz"), 1000, || {
        std::hint::black_box(model.polarizations(&below_cut, phi_ref, inclination));
    });
    time("polarizations, 2 frequencies (setup cost)", 10000, || {
        std::hint::black_box(model.polarizations(&frequencies[..2], phi_ref, inclination));
    });
    for (l, m) in [(2, -2), (2, -1), (3, -3), (3, -2), (4, -4)] {
        time(&format!("mode ({l},{m})"), 1000, || {
            std::hint::black_box(model.mode(l, m, &frequencies, phi_ref));
        });
    }
}

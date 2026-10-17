//! Timing of IMRPhenomXP on a 0.25 Hz grid from 20 to 1500 Hz, for the
//! SpinTaylor (PrecVersion 320) and NNLO (102) Euler angles.
//!
//! Compare with lalsimulation (single thread, multibanding off):
//!
//! ```python
//! import lal, lalsimulation as ls, numpy as np, timeit
//! p = lal.CreateDict(); ls.SimInspiralWaveformParamsInsertPhenomXPrecVersion(p, 320)
//! f = lal.CreateREAL8Vector(5921); f.data = 20 + 0.25 * np.arange(5921)
//! args = (0.1, 40 * lal.MSUN_SI, 20 * lal.MSUN_SI, 0.4, -0.3, 0.2, 0.1, 0.5, -0.3, 20.0,
//!         1e8 * lal.PC_SI, 0.7, p, ls.IMRPhenomXP, f)
//! timeit.timeit(lambda: ls.SimInspiralChooseFDWaveformSequence(*args), number=100) / 100
//! ```
use std::time::Instant;

use waveformrs::imrphenomx::{IMRPhenomXP, PrecessingBinary, PrecessionVersion};

fn time<F: FnMut()>(label: &str, iterations: u32, mut f: F) {
    f();
    let start = Instant::now();
    for _ in 0..iterations {
        f();
    }
    println!("{label:<55} {:>10.1?}", start.elapsed() / iterations);
}

fn main() {
    let binary = PrecessingBinary {
        mass_1: 40.0,
        mass_2: 20.0,
        spin_1: [0.4, -0.3, 0.2],
        spin_2: [0.1, 0.5, -0.3],
        luminosity_distance: 100.0,
        inclination: 0.7,
        phi_ref: 0.1,
        f_ref: 20.0,
    };
    let delta_f = 0.25;
    let frequencies: Vec<f64> = (0..5921).map(|i| 20.0 + i as f64 * delta_f).collect();
    println!("{} frequencies, M = 60 Msun, q = 0.5", frequencies.len());

    for (name, version) in [
        ("SpinTaylor (320)", PrecessionVersion::SpinTaylor),
        ("NNLO (102)", PrecessionVersion::Nnlo),
    ] {
        let model = IMRPhenomXP::with_version(&binary, version).unwrap();
        time(&format!("{name}: new"), 1000, || {
            std::hint::black_box(IMRPhenomXP::with_version(&binary, version).unwrap());
        });
        time(&format!("{name}: polarizations"), 200, || {
            std::hint::black_box(model.polarizations(&frequencies).unwrap());
        });
        time(&format!("{name}: new + polarizations"), 200, || {
            let model = IMRPhenomXP::with_version(&binary, version).unwrap();
            std::hint::black_box(model.polarizations(&frequencies).unwrap());
        });
        time(
            &format!("{name}: uniform grid (deltaF = 0.25, 20-1500 Hz)"),
            200,
            || {
                std::hint::black_box(model.polarizations_uniform(delta_f, 20.0, 1500.0).unwrap());
            },
        );
    }
}

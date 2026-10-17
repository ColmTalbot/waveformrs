//! Timing of IMRPhenomXPHM on a 0.25 Hz grid from 20 to 1500 Hz (the
//! configuration of `examples/imrxp.rs`: M = 60 Msun, q = 0.5, generic
//! spins), for the SpinTaylor (PrecVersion 320), MSA (223) and NNLO (102)
//! Euler angles.
//!
//! Compare with lalsimulation (single thread, `OMP_NUM_THREADS=1`):
//!
//! ```python
//! import lal, lalsimulation as ls, numpy as np, timeit
//! p = lal.CreateDict(); ls.SimInspiralWaveformParamsInsertPhenomXPrecVersion(p, 320)
//! ls.SimInspiralWaveformParamsInsertPhenomXHMThresholdMband(p, 0)
//! ls.SimInspiralWaveformParamsInsertPhenomXPHMThresholdMband(p, 0)
//! f = lal.CreateREAL8Vector(5921); f.data = 20 + 0.25 * np.arange(5921)
//! args = (0.1, 40 * lal.MSUN_SI, 20 * lal.MSUN_SI, 0.4, -0.3, 0.2, 0.1, 0.5, -0.3, 20.0,
//!         1e8 * lal.PC_SI, 0.7, p, ls.IMRPhenomXPHM, f)
//! timeit.timeit(lambda: ls.SimInspiralChooseFDWaveformSequence(*args), number=20) / 20
//! ```
use std::time::Instant;

use waveformrs::imrphenomx::{IMRPhenomXPHM, PrecessingBinary, PrecessionVersion};

fn time<F: FnMut()>(label: &str, iterations: u32, mut f: F) {
    f();
    let start = Instant::now();
    for _ in 0..iterations {
        f();
    }
    println!("{label:<58} {:>10.1?}", start.elapsed() / iterations);
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
        ("MSA (223)", PrecessionVersion::Msa),
        ("NNLO (102)", PrecessionVersion::Nnlo),
    ] {
        let model = IMRPhenomXPHM::with_version(&binary, version).unwrap();
        time(&format!("{name}: new"), 200, || {
            std::hint::black_box(IMRPhenomXPHM::with_version(&binary, version).unwrap());
        });
        time(&format!("{name}: polarizations"), 50, || {
            std::hint::black_box(model.polarizations(&frequencies).unwrap());
        });
        time(
            &format!("{name}: uniform grid (deltaF = 0.25, 20-1500 Hz)"),
            50,
            || {
                std::hint::black_box(model.polarizations_uniform(delta_f, 20.0, 1500.0).unwrap());
            },
        );
    }
}

use num_complex::Complex;
use std::time::Instant;

use waveformrs::{taylor::TaylorF2, waveform::Waveform};

fn main() {
    let total_mass = 90.0;
    let mass_ratio = 0.5;
    let chi_1 = 0.0;
    let chi_2 = 0.0;
    let luminosity_distance = 100.0;

    let delta_f: f64 = 0.25;
    let f_lower: f64 = 20.0;
    let f_max: f64 = 1500.0;
    let n_frequencies = (f_max / delta_f).ceil() as usize;

    let start = Instant::now();
    for _ in 0..1000 {
        let frequencies = (0..n_frequencies)
            .map(|i| i as f64 * delta_f)
            .collect::<Vec<f64>>();
        let taylor_f2 = TaylorF2::new(total_mass, mass_ratio, chi_1, chi_2, luminosity_distance);
        let _: Vec<Complex<f64>> = frequencies
            .iter()
            .map(|&f| {
                if f < f_lower || f > f_max {
                    return Complex::<f64>::new(0.0, 0.0);
                } else {
                    taylor_f2.waveform_single_frequency(f, 0.1)
                }
            })
            .collect();
    }
    let duration = start.elapsed();
    println!("Average evaluation time: {:?}", duration / 1000);
}

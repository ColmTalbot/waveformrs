use num_complex::Complex;
use std::f64::consts::PI;
use std::time::Instant;

use waveformrs::{imrphenomd::IMRPhenomD, waveform::Waveform};

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
        let imrd = IMRPhenomD::new(total_mass, mass_ratio, chi_1, chi_2, luminosity_distance);
        let _: Vec<Complex<f64>> = frequencies
            .iter()
            .map(|&f| {
                if f < f_lower || f > f_max {
                    return Complex::<f64>::new(0.0, 0.0);
                } else {
                    imrd.waveform_single_frequency(f, 0.1)
                }
            })
            .collect();
    }
    let duration = start.elapsed();
    println!("Average evaluation time: {:?}", duration / 1000);

    let imrd = IMRPhenomD::new(total_mass, mass_ratio, chi_1, chi_2, luminosity_distance);
    let value = imrd.waveform_single_frequency(20.0, 0.1);
    let phi_ref = value.im.atan2(value.re);
    println!("Inspiral");
    [19.5, 19.75, 20.0, 20.25, 20.5].iter().for_each(|&f| {
        let value = imrd.waveform_single_frequency(f, 0.1);
        println!(
            "{:.5}, {:.8e}",
            (value.im.atan2(value.re) - phi_ref + 2.0 * PI) % (2.0 * PI),
            value.norm()
        );
    });
    println!("Intermediate");
    let value = imrd.waveform_single_frequency(100.0, 0.1);
    let phi_ref = value.im.atan2(value.re);
    [99.5, 99.75, 100.0, 100.25, 100.5].iter().for_each(|&f| {
        let value = imrd.waveform_single_frequency(f, 0.1);
        println!(
            "{:.5}, {:.8e}",
            (value.im.atan2(value.re) - phi_ref + 2.0 * PI) % (2.0 * PI),
            value.norm()
        );
    });
    println!("Ringdown");
    let value = imrd.waveform_single_frequency(200.0, 0.1);
    let phi_ref = value.im.atan2(value.re);
    [199.5, 199.75, 200.0, 200.25, 200.5].iter().for_each(|&f| {
        let value = imrd.waveform_single_frequency(f, 0.1);
        println!(
            "{:.5}, {:.8e}",
            (value.im.atan2(value.re) - phi_ref + 2.0 * PI) % (2.0 * PI),
            value.norm()
        );
    });
}

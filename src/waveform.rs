use num_complex::Complex;

pub trait Waveform {
    fn orbital_speed(&self, frequency: f64) -> f64;
    fn phase(&self, v: f64, phi_c: f64) -> f64;
    fn amplitude(&self, v: f64) -> f64;

    fn waveform_single_frequency(&self, frequency: f64, phi_c: f64) -> Complex<f64> {
        let orbital_speed = self.orbital_speed(frequency);
        self.amplitude(orbital_speed)
            * Complex::<f64> {
                re: 0.0,
                im: -self.phase(orbital_speed, phi_c),
            }
            .exp()
    }
    fn waveform(&self, frequencies: &[f64], phi_c: f64) -> Vec<Complex<f64>> {
        frequencies
            .iter()
            .map(|&f| self.waveform_single_frequency(f, phi_c))
            .collect()
    }

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
}

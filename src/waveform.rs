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
        let waveform = self.waveform(frequencies, phi_c);
        let plus_factor = (1.0 + theta_jn.cos().powi(2)) / 2.0;
        let cross_factor = -Complex::I * theta_jn.cos();
        let plus = waveform.iter().map(|&w| w * plus_factor).collect();
        let cross = waveform.iter().map(|&w| w * cross_factor).collect();
        (plus, cross)
    }
}

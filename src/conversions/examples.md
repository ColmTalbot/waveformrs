# Examples

`SimInspiralTD` and `SimInspiralFD` for any approximant:

```
use waveformrs::conversions::simulate_td;
use waveformrs::{Approximant, Binary, WaveformOptions};

let binary = Binary {
    mass_1: 36.0,
    mass_2: 29.0,
    spin_1: [0.4, 0.2, 0.3],
    spin_2: [-0.2, 0.3, -0.1],
    luminosity_distance: 500.0,
    inclination: 0.7,
    phi_ref: 1.2,
    f_ref: 0.0,
};
// SimInspiralTD(..., deltaT = 1 / 4096, f_min = 20, f_ref = 0, IMRPhenomXPHM)
// with PhenomXPrecVersion = 320 and no multibanding
let options = WaveformOptions::default();
let h = simulate_td(Approximant::IMRPhenomXPHM, &binary, &options, 1.0 / 4096.0, 20.0)?;
println!(
    "{} samples from t = {} s",
    h.hplus.data.len(),
    h.hplus.epoch.seconds()
);
# Ok::<(), waveformrs::Error>(())
```

```no_run
use waveformrs::conversions::simulate_fd;
use waveformrs::{Approximant, Binary, WaveformOptions};

let binary = Binary {
    mass_1: 40.0,
    mass_2: 20.0,
    spin_1: [0.4, -0.3, 0.2],
    spin_2: [0.1, 0.5, -0.3],
    luminosity_distance: 100.0,
    inclination: 0.7,
    phi_ref: 0.1,
    f_ref: 20.0,
};
// SimInspiralFD(..., deltaF = 0, f_min = 20, f_max = 2048, f_ref = 20, NRSur7dq4)
// (needs NRSur7dq4_v1.0.h5 in LAL_DATA_PATH or NRSUR7DQ4_DATA)
let options = WaveformOptions::default();
let h = simulate_fd(Approximant::NRSur7dq4, &binary, &options, 0.0, 20.0, 2048.0)?;
println!("{} frequencies spaced by {} Hz", h.hplus.data.len(), h.hplus.delta_f);
# Ok::<(), waveformrs::Error>(())
```

Any frequency-domain generator, here a closure (a Gaussian pulse at
`t = 0`), through [`FdGeneratorFn`]:

```
use num_complex::Complex;
use waveformrs::conversions::{
    td_from_fd, ConditioningParameters, FdRequest, FdGeneratorFn,
    FrequencyDomainPolarizations, FrequencySeries, GpsTime,
};

let model = FdGeneratorFn {
    parameters: ConditioningParameters { mass_1: 30.0, mass_2: 30.0, spin_1z: 0.0, spin_2z: 0.0 },
    zero_reference_frequency_is_f_min: false,
    inspiral_only: false,
    generate: |request: &FdRequest| {
        let n = (request.f_max / request.delta_f) as usize + 1;
        let data: Vec<Complex<f64>> = (0..n)
            .map(|k| {
                let f = k as f64 * request.delta_f;
                Complex::new((-(f / 200.0).powi(2)).exp(), 0.0)
            })
            .collect();
        // LAL's frequency-domain models put t = 0 at the end of the segment
        let epoch = GpsTime::ZERO.add_seconds(-1.0 / request.delta_f);
        let series = FrequencySeries { epoch, delta_f: request.delta_f, data };
        Ok::<_, String>(FrequencyDomainPolarizations { hplus: series.clone(), hcross: series })
    },
};
let h = td_from_fd(&model, 1.0 / 2048.0, 20.0, 0.0).unwrap();
assert_eq!(h.hplus.delta_t, 1.0 / 2048.0);
```

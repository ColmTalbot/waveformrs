use std::f64::consts::PI;

use crate::constants::EGAMMA;

// Euler Gamma is currently only in nightly builds of the standard library
// see https://github.com/rust-lang/rust/issues/103883

pub(crate) fn chi_pn(params: &Params) -> f64 {
    let chi_s = (params.chi_1 + params.chi_2) / 2.0;
    let chi_a = (params.chi_1 - params.chi_2) / 2.0;

    return chi_s * (1.0 - params.eta * 76.0 / 113.0) + params.seta * chi_a;
}

#[derive(Clone, Copy, Debug)]
#[allow(dead_code)]
pub struct PNPhasing {
    pub v: [f64; 16],
    pub vlogv: [f64; 16],
    pub vlogvlogv: [f64; 16],
}

impl PNPhasing {
    pub fn new(params: &Params) -> Self {
        let scale = 3.0 / (128.0 * params.eta);
        let mut new = Self {
            v: [0.0; 16],
            vlogv: [0.0; 16],
            vlogvlogv: [0.0; 16],
        };
        new.v[0] = taylor_f2_phase_0(params) * scale;
        new.v[1] = taylor_f2_phase_1(params) * scale;
        new.v[2] = taylor_f2_phase_2(params) * scale;
        new.v[3] = taylor_f2_phase_3(params) * scale;
        new.v[4] = taylor_f2_phase_4(params) * scale;
        new.v[5] = taylor_f2_phase_5(params) * scale;
        new.v[6] = taylor_f2_phase_6(params) * scale;
        new.v[7] = taylor_f2_phase_7(params) * scale;
        new.v[10] = taylor_f2_phase_10(params) * scale;
        new.v[12] = taylor_f2_phase_12(params) * scale;
        new.v[13] = taylor_f2_phase_13(params) * scale;
        new.v[14] = taylor_f2_phase_14(params) * scale;
        new.v[15] = taylor_f2_phase_15(params) * scale;
        new.vlogv[5] = taylor_f2_phase_5l(params) * scale;
        new.vlogv[6] = taylor_f2_phase_6l(params) * scale;
        new
    }
}

#[derive(Clone, Copy, Debug)]
pub struct PNAmplitude {
    pub v: [f64; 10],
}

impl PNAmplitude {
    pub fn new(params: &Params) -> Self {
        let mut new = Self { v: [0.0; 10] };
        new.v[0] = taylor_f2_amplitude_0(params);
        new.v[1] = taylor_f2_amplitude_1(params);
        new.v[2] = taylor_f2_amplitude_2(params);
        new.v[3] = taylor_f2_amplitude_3(params);
        new.v[4] = taylor_f2_amplitude_4(params);
        new.v[5] = taylor_f2_amplitude_5(params);
        new.v[6] = taylor_f2_amplitude_6(params);
        new
    }
}

#[derive(Clone, Copy, Debug)]
pub struct Params {
    pub eta: f64,
    pub seta: f64,
    pub chi_1: f64,
    pub chi_2: f64,
    pub mass_ratio: f64,
    pub m1_on_m: f64,
    pub m2_on_m: f64,
    pub qm_def_1: f64,
    pub qm_def_2: f64,
    pub lambda_1: f64,
    pub lambda_2: f64,
}

impl Params {
    pub fn new(mass_ratio: f64, chi_1: f64, chi_2: f64) -> Self {
        let eta = mass_ratio / (1.0 + mass_ratio).powi(2);
        // seta = (m1 - m2) / (m1 + m2) = (1 - q) / (1 + q)
        let m1_on_m = 1.0 / (1.0 + mass_ratio);
        let m2_on_m = mass_ratio / (1.0 + mass_ratio);
        Self {
            eta,
            seta: m1_on_m - m2_on_m,
            chi_1,
            chi_2,
            mass_ratio,
            m1_on_m,
            m2_on_m,
            qm_def_1: 1.0,
            qm_def_2: 1.0,
            lambda_1: 0.0,
            lambda_2: 0.0,
        }
    }
}

pub fn _zero_function(_: &Params) -> f64 {
    0.0
}

pub fn taylor_f2_amplitude_0(_: &Params) -> f64 {
    1.0
}

pub fn taylor_f2_amplitude_1(_: &Params) -> f64 {
    0.0
}

pub fn taylor_f2_amplitude_2(args: &Params) -> f64 {
    -323.0 / 224.0 + 451.0 * args.eta / 168.0
}

pub fn taylor_f2_amplitude_3(args: &Params) -> f64 {
    args.chi_1 * (27.0 * args.mass_ratio / 16.0 - 11.0 * args.eta / 12.0 + 27.0 / 16.0)
        + args.chi_2 * (-27.0 * args.mass_ratio / 16.0 - 11.0 * args.eta / 12.0 + 27.0 / 16.0)
}

pub fn taylor_f2_amplitude_4(args: &Params) -> f64 {
    args.chi_1.powi(2) * (-81.0 * args.mass_ratio / 64.0 + 81.0 * args.eta / 32.0 - 81.0 / 64.0)
        + args.chi_2.powi(2)
            * (81.0 * args.mass_ratio / 64.0 + 81.0 * args.eta / 32.0 - 81.0 / 64.0)
        + (105271.0 / 24192.0 * args.eta.powi(2)
            - 1975055.0 / 338688.0 * args.eta
            - 27312085.0 / 8128512.0)
        - 47.0 / 16.0 * args.eta * args.chi_1 * args.chi_2
}

pub fn taylor_f2_amplitude_5(args: &Params) -> f64 {
    args.chi_1.powi(3)
        * (args.mass_ratio * (3.0 / 16.0 - 3.0 * args.eta / 16.0) - 9.0 * args.eta / 16.0
            + 3.0 / 16.0)
        + args.chi_1
            * (args.mass_ratio * (287213.0 / 32256.0 - 2083.0 * args.eta / 8064.0)
                - 2227.0 * args.eta.powi(2) / 2016.0
                - 15569.0 * args.eta / 1344.0
                + 287213.0 / 32256.0)
        + args.chi_2.powi(3)
            * (args.mass_ratio * (3.0 * args.eta / 16.0 - 3.0 / 16.0) - 9.0 * args.eta / 16.0
                + 3.0 / 16.0)
        + args.chi_2
            * (args.mass_ratio * (2083.0 * args.eta / 8064.0 - 287213.0 / 32256.0)
                - 2227.0 * args.eta.powi(2) / 2016.0
                - 15569.0 * args.eta / 1344.0
                + 287213.0 / 32256.0)
        - 85.0 * PI / 64.0
        + 85.0 * PI * args.eta / 16.0
}

pub fn taylor_f2_amplitude_6(args: &Params) -> f64 {
    args.chi_1
        * (-17.0 * PI * args.mass_ratio / 12.0 + 5.0 * PI * args.eta / 3.0 - 17.0 * PI / 12.0)
        + args.chi_2
            * (17.0 * PI * args.mass_ratio / 12.0 + 5.0 * PI * args.eta / 3.0 - 17.0 * PI / 12.0)
        + args.chi_1
            * args.chi_2
            * (-133249.0 * args.eta.powi(2) / 8064.0 - 319321.0 * args.eta / 32256.0)
        + args.chi_1.powi(2)
            * (args.mass_ratio * (-14139.0 * args.eta / 32256.0 - 49039.0 / 14336.0)
                + 163199.0 * args.eta.powi(2) / 16128.0
                + 158633.0 * args.eta / 64512.0
                - 49039.0 / 14336.0)
        + args.chi_2.powi(2)
            * (args.mass_ratio * (14139.0 * args.eta / 32256.0 + 49039.0 / 14336.0)
                + 163199.0 * args.eta.powi(2) / 16128.0
                + 158633.0 * args.eta / 64512.0
                - 49039.0 / 14336.0)
        - 177520268561.0 / 8583708672.0
        + (545384828789.0 / 5007163392.0 - 205.0 * PI.powi(2) / 48.0) * args.eta
        - 3248849057.0 * args.eta.powi(2) / 178827264.0
        + 34473079.0 * args.eta.powi(3) / 6386688.0
}

pub fn taylor_f2_phase_0(_: &Params) -> f64 {
    1.0
}

pub fn taylor_f2_phase_1(_: &Params) -> f64 {
    0.0
}

pub fn taylor_f2_phase_2(args: &Params) -> f64 {
    55.0 * args.eta / 9.0 + 3715.0 / 756.0
}

pub fn taylor_f2_phase_3(args: &Params) -> f64 {
    let mut phase = -16.0 * PI;
    for (m_on_m, chi) in vec![(args.m1_on_m, args.chi_1), (args.m2_on_m, args.chi_2)] {
        phase += m_on_m * (25.0 + 38.0 / 3.0 * m_on_m) * chi;
    }
    phase
}

pub fn taylor_f2_phase_4(args: &Params) -> f64 {
    let mut phase =
        15293365.0 / 508032.0 + 27145.0 / 504.0 * args.eta + 3085.0 / 72.0 * args.eta.powi(2);
    phase -= 395.0 / 4.0 * args.eta * args.chi_1 * args.chi_2;
    for (m_on_m, chi, qm_def) in vec![
        (args.m1_on_m, args.chi_1, args.qm_def_1),
        (args.m2_on_m, args.chi_2, args.qm_def_2),
    ] {
        phase -= (50.0 * qm_def + 5.0 / 8.0) * m_on_m.powi(2) * chi.powi(2);
    }
    phase
}

pub fn taylor_f2_phase_5(args: &Params) -> f64 {
    let mut phase = 5.0 / 9.0 * (7729.0 / 84.0 - 13.0 * args.eta) * PI;
    for (m_on_m, chi) in vec![(args.m1_on_m, args.chi_1), (args.m2_on_m, args.chi_2)] {
        phase -= chi
            * m_on_m
            * (13915.0 / 84.0 - m_on_m * (1.0 - m_on_m) * 10.0 / 3.0
                + m_on_m * (12760.0 / 81.0 + m_on_m * (1.0 - m_on_m) * 170.0 / 9.0));
    }
    phase
}

pub fn taylor_f2_phase_6(args: &Params) -> f64 {
    let mut phase =
        11583231236531.0 / 4694215680.0 - 640.0 / 3.0 * PI.powi(2) - 6848.0 / 21.0 * EGAMMA;
    phase += args.eta * (-15737765635.0 / 3048192.0 + 2255.0 / 12.0 * PI.powi(2));
    phase += args.eta.powi(2) * 76055.0 / 1728.0 - args.eta.powi(3) * 127825.0 / 1296.0;
    phase += taylor_f2_phase_6l(args) * 4.0f64.ln();
    phase += taylor_3pn_ss(args);
    for (m_on_m, chi) in vec![(args.m1_on_m, args.chi_1), (args.m2_on_m, args.chi_2)] {
        phase += PI * m_on_m * (1490.0 / 3.0 + m_on_m * 260.0) * chi;
    }
    phase
}

pub fn taylor_3pn_ss(args: &Params) -> f64 {
    let mut phase = (32675.0 / 112.0 + 5575.0 / 18.0 * args.eta) * args.chi_1 * args.chi_2;
    for (m_on_m, chi, qm_def) in vec![
        (args.m1_on_m, args.chi_1, args.qm_def_1),
        (args.m2_on_m, args.chi_2, args.qm_def_2),
    ] {
        phase += (47035.0 / 84.0 + 2935.0 / 6.0 * m_on_m - 120.0 * m_on_m.powi(2))
            * m_on_m.powi(2)
            * qm_def
            * chi.powi(2);
        phase += (-410825.0 / 672.0 - 1085.0 / 12.0 * m_on_m + 1255.0 / 36.0 * m_on_m.powi(2))
            * m_on_m.powi(2)
            * chi.powi(2);
    }
    phase
}

pub fn taylor_f2_phase_7(args: &Params) -> f64 {
    let mut phase = PI
        * (77096675.0 / 254016.0 + 378515.0 / 1512.0 * args.eta
            - 74045.0 / 756.0 * args.eta.powi(2));
    for (m_on_m, chi) in vec![(args.m1_on_m, args.chi_1), (args.m2_on_m, args.chi_2)] {
        phase += chi
            * m_on_m
            * (-170978035.0 / 48384.0
                + args.eta * 2876425.0 / 672.0
                + args.eta.powi(2) * 4735.0 / 144.0
                + m_on_m
                    * (-7189233785.0 / 1524096.0 + args.eta * 458555.0 / 3024.0
                        - args.eta.powi(2) * 5345.0 / 72.0));
    }
    phase
}

pub fn taylor_f2_phase_8(_: &Params) -> f64 {
    0.0
}

pub fn taylor_f2_phase_9(_: &Params) -> f64 {
    0.0
}

pub fn taylor_f2_phase_10(args: &Params) -> f64 {
    let mut phase = 0.0;
    for (lambda, m_on_m) in vec![(args.lambda_1, args.m1_on_m), (args.lambda_2, args.m2_on_m)] {
        phase += 24.0 * (-12.0 + 11.0 * m_on_m) * m_on_m.powi(4) * lambda;
    }
    phase
}

pub fn taylor_f2_phase_11(_: &Params) -> f64 {
    0.0
}

pub fn taylor_f2_phase_12(args: &Params) -> f64 {
    let mut phase = 0.0;
    for (lambda, m_on_m) in vec![(args.lambda_1, args.m1_on_m), (args.lambda_2, args.m2_on_m)] {
        phase += (-15895.0 / 28.0 + 4595.0 / 28.0 * m_on_m + 5715.0 / 14.0 * m_on_m.powi(2)
            - 325.0 / 7.0 * m_on_m.powi(3))
            * m_on_m.powi(4)
            * lambda;
    }
    phase
}

pub fn taylor_f2_phase_13(args: &Params) -> f64 {
    let mut phase = 0.0;
    for (lambda, m_on_m) in vec![(args.lambda_1, args.m1_on_m), (args.lambda_2, args.m2_on_m)] {
        phase += 24.0 * (12.0 - 11.0 * m_on_m) * PI * m_on_m.powi(4) * lambda;
    }
    phase
}

pub fn taylor_f2_phase_14(args: &Params) -> f64 {
    let mut phase = 0.0;
    for (lambda, m_on_m) in vec![(args.lambda_1, args.m1_on_m), (args.lambda_2, args.m2_on_m)] {
        phase += -(m_on_m.powi(4))
            * lambda
            * 5.0
            * (193986935.0 / 571536.0
                - 14415613.0 / 381024.0 * m_on_m
                - 57859.0 / 378.0 * m_on_m.powi(2)
                - 209495.0 / 1512.0 * m_on_m.powi(3)
                + 965.0 / 54.0 * m_on_m.powi(4)
                - 4.0 * m_on_m.powi(5));
    }
    phase
}

pub fn taylor_f2_phase_15(args: &Params) -> f64 {
    let mut phase = 0.0;
    for (lambda, m_on_m) in vec![(args.lambda_1, args.m1_on_m), (args.lambda_2, args.m2_on_m)] {
        phase += (m_on_m.powi(4)) * lambda * 1.0 / 28.0
            * PI
            * (27719.0 - 22415.0 * m_on_m + 7598.0 * m_on_m.powi(2) - 10520.0 * m_on_m.powi(3));
    }
    phase
}

pub fn taylor_f2_phase_0l(_: &Params) -> f64 {
    0.0
}

pub fn taylor_f2_phase_1l(_: &Params) -> f64 {
    0.0
}

pub fn taylor_f2_phase_2l(_: &Params) -> f64 {
    0.0
}

pub fn taylor_f2_phase_3l(_: &Params) -> f64 {
    0.0
}

pub fn taylor_f2_phase_4l(_: &Params) -> f64 {
    0.0
}

pub fn taylor_f2_phase_5l(args: &Params) -> f64 {
    taylor_f2_phase_5(args) * 3.0
}

pub fn taylor_f2_phase_6l(_: &Params) -> f64 {
    -6848.0 / 21.0
}

pub fn taylor_f2_phase_7l(_: &Params) -> f64 {
    0.0
}

pub fn taylor_f2_phase_8l(_: &Params) -> f64 {
    0.0
}

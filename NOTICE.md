# Notices

waveformrs is Copyright (C) 2024-2026 Colm Talbot.

waveformrs is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version. It is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
The full text is in [LICENSE](LICENSE).

## Code this work is derived from

Much of waveformrs is a translation into Rust of C code from the projects
below, modified (restructured, optimised and extended) in 2024-2026. The
copyright of the original code stays with its authors, listed as in the
original files. The licence of each project allows its use in a work
licensed under the GNU General Public License version 3 or later, which is
why waveformrs uses that licence.

### LALSuite

[LALSuite](https://git.ligo.org/lscsoft/lalsuite): lalsimulation 6.2.0 and
lal 7.7.0, and the changes of lalsuite merge request !2593 (LALSuite master).
Licensed under the GNU General Public License version 2 or (at your option)
any later version, used here under version 3 or later.

| waveformrs | LALSuite files | Copyright |
|---|---|---|
| `src/taylor.rs`, `src/pn.rs` | `LALSimInspiralTaylorF2.c`, `LALSimInspiralPNCoefficients.c` | (C) 2007 Jolien Creighton, B.S. Sathyaprakash, Thomas Cokelaer; (C) 2012 Leo Singer, Evan Ochsner, Les Wade, Alex Nitz; (C) 2011 Drew Keppel, 2012 Riccardo Sturani |
| `src/imrphenomd*` | `LALSimIMRPhenomD.c`, `LALSimIMRPhenomD_internals.c` | (C) 2015 Michael Puerrer, Sebastian Khan, Frank Ohme, Ofek Birnholtz, Lionel London, Francesco Pannarale |
| `src/imrphenomx*` (IMRPhenomXAS) | `LALSimIMRPhenomX.c`, `LALSimIMRPhenomX_internals.c`, `LALSimIMRPhenomXUtilities.c`, `LALSimIMRPhenomX_inspiral.c`, `LALSimIMRPhenomX_intermediate.c`, `LALSimIMRPhenomX_ringdown.c`, `LALSimIMRPhenomX_qnm.c` | (C) 2018 Geraint Pratten; (C) 2018 Sascha Husa, Geraint Pratten |
| `src/imrphenomx*` (IMRPhenomXHM) | `LALSimIMRPhenomXHM.c`, `LALSimIMRPhenomXHM_internals.c`, `LALSimIMRPhenomXHM_inspiral.c`, `LALSimIMRPhenomXHM_intermediate.c`, `LALSimIMRPhenomXHM_ringdown.c`, `LALSimIMRPhenomXHM_qnm.c` | (C) 2019 Marta Colleoni, Cecilio García Quirós |
| `src/imrphenomx*` (IMRPhenomXP, IMRPhenomXPHM) | `LALSimIMRPhenomXPHM.c`, `LALSimIMRPhenomX_precession.c` | (C) 2019 Cecilio García Quirós, Geraint Pratten; (C) 2018/2019 Geraint Pratten |
| `src/imrphenomx*` (IMRPhenomXPNR) | `LALSimIMRPhenomX_PNR.c`, `LALSimIMRPhenomX_PNR_alpha.c`, `LALSimIMRPhenomX_PNR_beta.c`, `LALSimIMRPhenomX_PNR_coefficients.c`, `LALSimIMRPhenomX_PNR_internals.c`, `LALSimIMRPhenomX_AntisymmetricWaveform.c`; `LALSimIMRPhenomX_PNR_deviations.c` | (C) 2022 Cardiff University; (C) 2021 The University of Amsterdam |
| `src/nrtidal.rs`, `src/imrphenomx/nrtidal.rs`, the NRTidal parts of `src/imrphenomd.rs` and `src/imrphenomx/` | `LALSimNRTunedTides.c`, `LALSimIMRPhenomD_NRTidal.c`, `LALSimUniversalRelations.c`, `LALSimInspiralEOS.c`, `LALSimIMRPhenomX.c`, `LALSimIMRPhenomX_internals.c` | (C) 2017 Tim Dietrich, Sebastiano Bernuzzi, Nathan Johnson-McDaniel, Shasvath J Kapadia, Francesco Pannarale and Sebastian Khan, Michael Puerrer, Adrian Abac; (C) 2015 Sebastian Khan, Michael Puerrer; (C) 2017 Andrea Taracchini; (C) 2012 Walter Del Pozzo, Tjonnie Li, Michalis Agathos; (C) 2018 Geraint Pratten |
| `src/spintaylor/`, `src/imrphenomx/spintaylor_angles.rs` | `LALSimInspiralSpinTaylor.c`; lal `LALAdaptiveRungeKuttaIntegrator.c` | (C) 2011 E. Ochsner, 2014 A. Klein, 2015 R. Sturani; (C) 2010 Michele Vallisneri, Will Farr, Evan Ochsner, 2014 A. Klein |
| `src/nrsur*` | `LALSimIMRPrecessingNRSur.c`, `LALSimNRSurrogateUtilities.c` | (C) 2017 Jonathan Blackman, Vijay Varma; (C) 2018 Jonathan Blackman |
| `src/conditioning.rs`, `src/conversions/`, `src/series.rs`, `src/fft.rs` | `LALSimInspiral.c`, `LALSimInspiralGeneratorConditioning.c`, `LALSimInspiralGeneratorLegacy.c`; lal `ButterworthTimeSeries.c`, `TimeFreqFFT.c` | (C) 2008 J. Creighton, S. Fairhurst, B. Krishnan, L. Santamaria, D. Keppel, Evan Ochsner, C. Pankow, 2014 A. Klein; (C) 2007 Jolien Creighton, Teviet Creighton; (C) 2007 Duncan Brown, Jolien Creighton, Patrick Brady, Reinhard Prix, Tania Regimbau, John Whelan (the conditioning files carry no notice of their own) |
| `src/imrphenomx/wigner.rs`, `src/model/modes.rs` | lal `SphericalHarmonics.c`; `LALSimSphHarmMode.c` | (C) 2007 S. Fairhurst, B. Krishnan, L. Santamaria, C. Robinson; (C) 2015 Jolien Creighton |
| `src/constants.rs` | lal `LALConstants.h` | (C) 2007 Jolien Creighton |

### GNU Scientific Library (GSL)

[GSL](https://www.gnu.org/software/gsl/) 2.8, licensed under the GNU General
Public License version 3 or (at your option) any later version.

| waveformrs | GSL files | Copyright |
|---|---|---|
| `src/imrphenomx/msa.rs` (Jacobi elliptic functions, elliptic integral of the first kind, `hypot`) | `specfunc/elljac.c`, `specfunc/ellint.c`, `sys/hypot.c` | (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman; (C) 1996, 1997, 1998, 1999, 2000, 2007 Gerard Jungman, Brian Gough, Patrick Alken |
| `src/interpolate.rs`, `src/nrsur/spline.rs` (natural cubic splines) | `interpolation/cspline.c` | (C) 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2004 Gerard Jungman |
| `src/imrphenomx/poly_roots.rs` (polynomial roots) | `poly/zsolve.c`, `poly/companion.c`, `poly/balance.c`, `poly/qr.c` | (C) 1996, 1997, 1998, 1999, 2000, 2007 Brian Gough |
| `src/imrphenomx/linalg.rs` (LU decomposition) | `linalg/lu.c` | (C) 1996, 1997, 1998, 1999, 2000, 2007, 2009 Gerard Jungman, Brian Gough; (C) 2019, 2021 Patrick Alken |
| `src/spintaylor/integrator.rs` (Runge-Kutta-Fehlberg integration) | `ode-initval/rkf45.c`, `ode-initval/evolve.c`, `ode-initval/cstd.c` | (C) 2001, 2004, 2007 Brian Gough; (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman |

### fdlibm

The sine, cosine and logarithm kernels and constants in `src/fastmath.rs`
follow [fdlibm](https://www.netlib.org/fdlibm/), under this notice:

> Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
>
> Developed at SunPro, a Sun Microsystems, Inc. business.
> Permission to use, copy, modify, and distribute this
> software is freely granted, provided that this notice
> is preserved.

## Dependencies

The Rust crates waveformrs depends on (among them `num-complex`, `once_cell`,
`rustfft`, and for the Python bindings `pyo3` and `numpy`) are licensed under
MIT, Apache-2.0, BSD-2-Clause and similar permissive licences, all compatible
with the GNU General Public License version 3. They are not included in this
repository. A binary distribution, such as a Python wheel, contains their
code and should also include their licence notices.

## Data

The reference data in `tests/data` is the output of LALSimulation, generated
by the scripts in `scripts/`. The NRSur7dq4 data file (`NRSur7dq4_v1.0.h5`)
is not part of this repository; it is distributed with LALSuite's data.

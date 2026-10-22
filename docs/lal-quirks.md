# Where waveformrs differs from LALSimulation, and the quirks it keeps

waveformrs ports the models of LALSimulation 6.2.0. This page lists every
place where the port deliberately does something other than what LAL's code
does, and the oddities of LAL's code that it still reproduces. The LAL
references are to lalsuite at `ba91eb79ca` (May 2025), which for these
functions is the same as the 6.2.0 release.

Each entry says what LAL does, what waveformrs does, and how much it
matters. "Fixed" means the port no longer copies LAL here; "kept" means it
does, with the reason.

Several entries concern the original release of the IMRPhenomXHM
amplitudes, `PhenomXHMReleaseVersion = 122019`. IMRPhenomXHM and
IMRPhenomXPHM use the recalibrated release 122022 (LAL's default), and
waveformrs offers no other. LAL builds IMRPhenomXPNR (XO4a) on the 122019
amplitudes ("XO4 uses previous version of XHM",
`LALSimInspiralGeneratorLegacy.c`), and so does waveformrs: there the
122019 code is used, quirks included, so that IMRPhenomXPNR stays LAL's
model.

## Summary

| | what | status | effect |
|---|---|---|---|
| 1 | elementary functions, fused multiply-adds, analytic derivatives | changed | ~1e-16 per operation; up to ~1e-4 where LAL amplifies rounding |
| 2 | IMRPhenomXPNR antisymmetric phase (parity) | fixed (lalsuite !2593) | large for the antisymmetric mode |
| 3 | (3,2) mixing coefficients in single precision | fixed | ~1e-7 in the (3,2) ringdown |
| 4 | `I/3` in single precision | fixed | ~3e-8 in one PN coefficient of the (2,1) amplitude (122019, IMRPhenomXPNR) |
| 5 | (2,1) high-spin intermediate phase without the 22 connection | fixed | where the collocation points straddle a 22 region boundary |
| 6 | NRTidalv3 `(4/3)` integer division | fixed | none on the waveforms |
| 7 | spline status: only the last evaluation checked | fixed | only near the ends of the splines |
| 8 | quintic zero-crossing test reads the wrong roots | fixed | 2% of binaries (122019, IMRPhenomXPNR) |
| 9 | PNR moving average in place, weighted | fixed (the paper's Eq. 17) | up to ~3e-3 on uniform grids, ~7e-3 on non-uniform ones |
| 10 | quartic zero-crossing test drops an imaginary part | kept | none seen (122019, IMRPhenomXPNR) |
| 11 | PNR inflection point: `acos` of a complex value | kept | none seen on the waveforms |
| 12 | `eta >= 0.0237954` against the comment "q > 40" | kept | 122019, IMRPhenomXPNR |
| 13 | `>` versus `>=` at the same 22 phase boundary | kept | a single frequency |
| 14 | smaller C-isms | kept | none |
| 15 | negative higher-mode amplitudes allowed | kept | 122019, IMRPhenomXPNR |
| 16 | IMRPhenomXPNR antisymmetric mode against its paper | noted | see entry 16 |

## Changes with respect to LAL

### 1. Elementary functions and fused multiply-adds

LAL's results depend on the C library and the compiler of each build:
glibc's `cbrt` on Linux and Apple's on macOS, and clang's or GCC's choice of
which products and sums to fuse. waveformrs used to copy one build (the
conda-forge osx-arm64 one) operation by operation. It now uses the best
available operations instead:

* `f64::cbrt` (correctly rounded on linux-gnu, Apple's on macOS), and
  `cbrt` instead of `pow(x, 1./3.)`;
* fused multiply-adds wherever the hardware has the instruction;
* the complex exponential, square root and principal cube root in Rust
  (`src/complex.rs`), the cube root from `cbrt(|z|)` rather than
  `cpow(z, 1./3.)`;
* `f64::sin_cos` in NRSur7dq4 rather than the calls that matched clang's
  use of `__sincos_stret`;
* the NRTidal quadrupole-monopole polynomial in Horner's form;
* the fractional powers of the IMRPhenomX frequencies (`UsefulPowers`)
  from `cbrt` and `sqrt` rather than LAL's repeated multiplication of
  `pow(x, 1/6)`, which accumulates several ulp in the higher powers (LAL's
  own inspiral phase is off by up to ~9 ulp of its largest term against an
  exact evaluation, see `tests/long_waveforms.rs`);
* one form of each IMRPhenomX ansatz (nested, with fused multiply-adds),
  instead of a copy of LAL's operations for the coefficients next to a
  faster one for the frequencies;
* exact constants: pi/2 (LAL's literal is 2 ulp off) and sqrt(10) (1 ulp
  off) in IMRPhenomXP/XPHM;
* fused multiply-adds in the NRSur7dq4 frame rotation (in the vectorised
  and per-sample paths alike), the empirical interpolants accumulated
  column by column without copying OpenBLAS's `dgemv` kernel, and the
  IMRPhenomX LU solves without copying OpenBLAS's `dger`/`dtrsv` kernels.

Kept on purpose: compiler-rt's scaled complex division (which avoids
overflow), complex rotations `exp(i k alpha)` built by repeated
multiplication (accurate to a few ulp for `k <= 4`, and much cheaper than
`sin_cos` per `k`), and the operation order of the GSL algorithms that
LAL calls (cubic splines, the RKF45 integrator, polynomial roots), which
are part of how the models are defined.

The results are mathematically the same. The differences from LAL are at
the rounding level except where the models amplify rounding: the (3,2)
mode of IMRPhenomXHM, which takes a second derivative by finite
differences with a step of 1e-7, the MSA angles and the SpinTaylor
integration. There the crate differs from LAL by about as much as LAL's
own builds differ from each other (up to ~2e-3 in the (3,2) mode).
`docs/linux.md` has the measurements.

LAL fixes several IMRPhenomXHM coefficients from derivatives taken by
finite differences: the second derivative of the rotated spheroidal
ringdown phase of the (3,2) mode (step 1e-7), and the amplitude
derivatives at the region boundaries (step 1e-9). waveformrs uses the
analytic derivatives of the same functions (`src/imrphenomx/jet.rs`). Its
own builds on macOS and Linux then agree to ~1e-12 in IMRPhenomXHM
instead of ~1e-5, and it differs from LAL by about as much as LAL's builds
differ from each other (CHANGELOG.md has the numbers).

### 2. IMRPhenomXPNR: the parity of the antisymmetric phase

LALSimulation 6.2.0 computes the phase of the antisymmetric (2,2) mode of
IMRPhenomXPNR in a frame that breaks the mirror symmetry of the binary.
lalsuite merge request !2593 (commit `64d2004def`, in lalsuite master but
not in a release) fixes this. waveformrs implements only the fixed model, so
its IMRPhenomXPNR differs from that of the 6.2.0 release. It is tested
against a liblalsimulation 6.2.0 built with the fix
(`scripts/build_patched_lalsimulation.sh`, the `imrphenomxpnr_parityfix*`
suites).

### 3. The (3,2) mixing coefficients in double precision

`IMRPhenomXHM_Initialize_MixingCoeffs` stores the spheroidal-to-spherical
mixing coefficients in `COMPLEX8 mixingCoeffs[4]`
(`LALSimIMRPhenomXHM_structs.h:140`), so the fits are rounded to single
precision. That is a relative error of ~6e-8 in each coefficient, and
~1e-7 in the (3,2) ringdown amplitude and phase. waveformrs keeps them in
double precision (`mixing_coefficients` in
`src/imrphenomx/hm_internals.rs`).

### 4. `I/3` in the (2,1) PN amplitude

`LALSimIMRPhenomXHM_internals.c:1098` writes `I/3*delta*(...)`. `I` is
`_Complex_I`, a `complex float`, and `3` an `int`, so `I/3` is computed in
single precision: `1/3` with a relative error of 3e-8 in one 1.5PN
coefficient (`x2`) of the (2,1) inspiral amplitude of the 122019 release.
waveformrs uses `1/3` in double precision.

### 5. The (2,1) intermediate phase at high spins

For the (2,1) mode with `STotR >= 0.8`, LAL extends the inspiral into the
intermediate region with differences of the 22 phase derivative,
`IMRPhenomX_dPhase_22`, at three collocation frequencies
(`LALSimIMRPhenomXHM_internals.c:2076`). It calls this before it computes
the 22 connection coefficients, which are still zero, so the derivative is
the uncorrected one of each 22 region and jumps at the 22 region
boundaries. If two collocation points lie in different 22 regions, the
difference includes the jump. waveformrs computes the 22 connection
coefficients first, so the derivative is continuous.

### 6. NRTidalv3: `(4/3)` is an integer division

`IMRPhenomX_TidalPhaseDerivative` (`LALSimIMRPhenomX_internals.c`, lines
3234–3235 in 6.2.0) writes

```c
+ (4/3)*pow(Mf, 1.0/3.0)*c_2A*powers_of_lalpi.four_thirds
```

for stars A and B. In C, `(4/3)` is 1, so the derivative of
`c_2 (pi Mf)^(4/3)` is too small by 3/4. Neighbouring terms use `2./3.` and
`5./3.`, so this is a typo. It came in with NRTidalv3 (lalsuite
`f1fa96d137`, December 2023).

It has no effect on LAL's waveforms. The derivative enters only through
the Planck taper that switches from the NRTidal phase to the PN phase
between 1.15 and 1.35 times the NRTidalv3 merger frequency. LAL evaluates
it only at `f_final = min(last frequency, f_merger)`, to set the time shift
of IMRPhenomXAS_NRTidalv3 and IMRPhenomXP_NRTidalv3. There the taper and its
derivative both vanish, so the term is multiplied by zero. It would matter
only if the derivative were evaluated above `1.15 f_merger`.

waveformrs uses `4/3`. The waveforms are unchanged.

### 7. Spline evaluations: every status is checked

The PNR and SpinTaylor `beta` connections evaluate the `cos(beta)` spline
at `f1 = 0.97 fmax`, `f2 = 0.98 fmax` and `fmax`, and fall back to a flat
`beta` if an evaluation fails. LAL keeps the status of the last evaluation
only (`LALSimIMRPhenomX_PNR_*.c`, `IMRPhenomX_SpinTaylorAnglesSplinesAll`),
so a point outside the spline at `f1` or `f2` goes unnoticed (and gives an
extrapolated value or NaN). waveformrs requires all three to be inside the
spline. This changes only binaries where the spline starts above
`0.97 fmax`.

### 8. The quintic zero-crossing test reads the wrong roots (fixed)

`CrossZeroP5` (`LALSimIMRPhenomXHM_intermediate.c:2150`) decides whether the
quintic intermediate amplitude polynomial of IMRPhenomXHM crosses zero
between the inspiral and the ringdown. If it does, `ChoosePolOrder` lowers
the order of the polynomial and tries again. `gsl_poly_complex_solve`
returns the five roots packed as `z[2k]` (real part) and `z[2k+1]`
(imaginary part), and LAL tests

| root | imaginary part tested | real part tested | correct |
|---|---|---|---|
| 1 | `z[1]` | `z[0]` | yes |
| 2 | `z[3]` | `z[2]` | yes |
| 3 | `z[4]` | `z[4]` | `z[5]`, `z[4]` |
| 4 | `z[5]` | `z[6]` | `z[7]`, `z[6]` |
| 5 | `z[7]` | `z[8]` | `z[9]`, `z[8]` |

so the third root is tested against its own real part as an imaginary
part, and the fourth and fifth with the imaginary parts of the root before
them. A real root in the interval can be missed, and the real part of a
complex root can be taken for a crossing. The code is unchanged since
IMRPhenomXHM was added (`f57c060a66`, February 2020). Its comments and
those of `ChoosePolOrder` describe checking all the roots, so this is an
indexing slip, not a choice.

It affects only the 122019 amplitudes, which waveformrs uses for
IMRPhenomXPNR. Measured on the aligned-spin IMRPhenomXHM modes of that
release over 20,000 binaries (mass ratio log-uniform in 1–1000, aligned
spins uniform in ±0.99), the corrected test decides differently in 840 of
71,061 calls, which changes 404 binaries (2%). The affected modes change by
3e-4 to 3e-3 relative to their peak, at most 2e-2.

waveformrs tests every root with its own real and imaginary parts
(`cross_zero_p5` in `src/imrphenomx/hm_amp_122019.rs`). Its IMRPhenomXPNR
therefore differs from LAL's in those binaries. Of the reference cases, one
mode differs from LAL by more than rounding (the (3,3) of `m1 = 148.5`,
`m2 = 3.5`, by 8e-4); the test lists it.

### 9. The PNR antisymmetric amplitude ratio: the moving average of the paper (fixed)

`IMRPhenomX_PNR_GenerateAntisymmetricAmpRatio`
(`LALSimIMRPhenomX_AntisymmetricWaveform.c`) smooths the ratio `kappa` of
the antisymmetric to the symmetric amplitude, which is the PN fit of Eq. 16
of Ghosh et al. 2024 (arXiv:2310.16980) below the ringdown frequency and
its ringdown value above it:

```c
for (size_t id = 0; id < kappa->length-width-1; id++) {
    ... smoothed_ratio += kappa->data[id+j]*df ...   /* j = 0 .. width = 80 */
    kappa->data[id+half_width] = smoothed_ratio / frequency_width;
}
```

The paper (Sec. IV, Eq. 17) describes this as a moving average over a
symmetric window of k = 40 points on each side, for 40 <= n <= N - 40,
"to avoid a sharp transition". It is written as a plain average of the
unsmoothed values. LAL's code differs from it in three ways:

* it writes into the array it reads, so each window's left half already
  holds smoothed values: the filter is recursive, and smooths more than
  the paper's average;
* it weights each sample by the frequency step (the same as the paper's
  average on a uniform grid);
* it leaves the first 40 and the last 41 samples unsmoothed (the paper:
  the first and last 40).

The paper's Eq. 17 (its sum runs from n - k + 1 to n + k + 1) is itself off
by one from the symmetric window its text describes.

Effect of a plain (non-recursive) average, over 300 IMRPhenomXPNR binaries
(mass ratio 1–20, total mass 20–200, isotropic spins up to 0.99), as
max|h_+ change| / max|h_+| (median, 90%, max) and white-noise mismatch:

| grid | relative change | mismatch |
|---|---|---|
| uniform, 0.125 Hz | 4e-5, 4e-4, 1.6e-3 | 2e-9, 6e-8, 4e-7 |
| uniform, 1 Hz | 2e-4, 9e-4, 2.9e-3 | 8e-8, 1.4e-6, 8e-6 |

The window is a number of samples, not a frequency width, in the paper
and in LAL, so the model depends on the frequency grid it is evaluated
on. The same binaries on the 1 Hz grid against the 0.125 Hz grid (at the
same frequencies) differ by 1e-3 (median), 8e-3 (90%) and 3.5e-2 (max)
with LAL's average, and by about as much (1e-3, 8e-3, 3.2e-2) with the
plain one. The recursion is therefore a small part of the model's
uncertainty; the grid dependence is the larger issue, and would need a
window of fixed width in frequency, which neither the paper nor LAL uses.

waveformrs uses the average of the paper (`smooth_amp_ratio` in
`src/imrphenomx/pnr.rs`): each sample `n` with `40 <= n < N - 40` becomes
the plain mean of the 81 unsmoothed samples `n - 40 ..= n + 40`, centred as
the text describes (not on `n + 1` as Eq. 17's indices are written), and
the first and last 40 samples are left as they are. On the reference
binaries this moves IMRPhenomXPNR from LAL by up to 2e-3 (mismatch 3e-6)
on uniform grids, and by up to 7e-3 (mismatch 3e-4) on the non-uniform
frequency sequences of the reference data. On such grids the paper's
average, a mean over samples, and LAL's, weighted by the frequency
steps, differ the most; the paper only considers the uniform grids of its
numerical-relativity data. The unit tests in `src/imrphenomx/xpnr` switch
to LAL's average and compare the rest of the model with LAL at the strict
tolerances.

## Quirks that are kept

### 10. The quartic zero-crossing test drops an imaginary part

`CrossZeroP4` (Ferrari's formula) computes
`bigdenom = 4*sqrt(q1 + frac1 + frac2)` with the real `sqrt` of a complex
value, which C converts to its real part, while the same expression a line
later is `csqrt(q1 + frac1 + frac2)`. This is a slip: `csqrt` was meant.
Over the same 20,000 binaries, `csqrt` changes none of the 7,890
decisions. Kept (122019, IMRPhenomXPNR).

### 11. The PNR inflection point: `acos` of a complex value

`IMRPhenomX_PNR_three_inflection_points` (`LALSimIMRPhenomX_PNR_beta.c`)
solves the cubic of Eq. 68 of arXiv:2107.08876 with the trigonometric
formula in `COMPLEX16` arithmetic, but calls `acos` and `cos`, not `cacos`
and `ccos`. C passes their real parts. When the cubic has a single real
root, the argument of `acos` is imaginary or outside [-1, 1], and the
complex formula would be needed. The selection that follows also tests
`cimag(f) < 1e-10` without `fabs`. Over 3,000 precessing binaries (mass
ratio 1–20, isotropic spins up to 0.99), the complex formula changes the
inflection frequency beyond rounding in 4 of 6,010 evaluations (2
binaries), and the waveform in none. Kept.

### 12. The (2,1) PN amplitude condition

For the 122019 release, `IMRPhenomXHM_SetHMWaveformVariables` sets
`useFAmpPN = 1`, the (2,1) inspiral amplitude ansatz without re-expansion,
if `eta >= 0.0237954`, with the comment "for EMR (q>40.)". `eta =
0.0237954` is `q = 40`, so the condition selects `q <= 40`, the opposite
of the comment. The code has been the same since 2020 and was kept in the
2023 recalibration (`9eaa372236`); the inspiral code says the (2,1) is
"better behaved without the reexpansion". Which of the two was intended
needs the model's authors or arXiv:2001.10914. Kept.

### 13. `>` versus `>=` at the same phase boundary

`IMRPhenomX_Phase_22` changes region at `Mf >= fPhaseMatchIM`, while the
per-frequency phase of IMRPhenomXAS changes at `Mf > fPhaseMatchIM`.
This matters at one frequency at most, where the two regions agree to the
connection's accuracy. Kept.

### 14. Smaller C-isms

These are copied as written and change nothing: `0.0 * x` terms from
complex multiplications by `I`, `(int)(tend - tinit) / deltat` for a buffer
size, and LAL's integer truncations of sample counts (part of LAL's
interface, so not quirks).

### 15. Negative higher-mode amplitudes in release 122019

The 122022 release replaces a negative mode amplitude by a tiny positive
value (`if (Amp < 0 && pWF->IMRPhenomXHMReleaseVersion != 122019) Amp =
FALSE_ZERO`, `LALSimIMRPhenomXHM_internals.c`); the 122019 release lets it
go negative, which flips the sign of the mode there. `ChoosePolOrder`
(entry 8) is meant to prevent this in the intermediate region. Kept for
IMRPhenomXPNR.

### 16. IMRPhenomXPNR's antisymmetric mode against its paper

A comparison of Ghosh et al. 2024 (arXiv:2310.16980, Secs. IV and V) with
LAL master (`LALSimIMRPhenomX_AntisymmetricWaveform.c`) and waveformrs.
Apart from the moving average (entry 9), waveformrs follows LAL master.

Where the code differs from the paper:

* **The reference phase `phi_A0`.** The paper sets it to `alpha(f_ref)`.
  LAL master (with lalsuite !2593) uses `-alpha(f_ref) - pi/2` and adds
  `phiJ_Sf` to the whole antisymmetric phase; lalsimulation 6.2 used
  `alpha(f_ref) - alpha0` and `zeta_polarization`, closer to the paper. The
  paper's sign looks wrong: it would double the `alpha` term at `f_ref`
  instead of removing it. The other terms are the code's frame convention.
* **Eq. 25 as printed is not continuous with Eq. 24.** The code uses the
  continuous `phi_B0 = alpha(fT) - phi_s(fT) / 2 + A0 fT + phi_A0`.
* **The window of Eqs. 21–22 (`fw = 0.005`) is not implemented.** The code
  switches from one branch of Eq. 24 to the other at `fT = 0.85 fRD`. With
  `A0` from Eq. 23 the phase and its derivative are continuous there, so
  the window would change little.
* **The cutoff at `Mf = 0.12`.** The paper gives the model below it; the
  code continues `phi_s + phi_B0` up to the end of the waveform.
* **Eq. 23 is evaluated at `fT`**, as written, although the text mentions
  `fRD`; `alpha'` is a central difference with `dMf = 0.0005`, which the
  paper does not specify.

Where the paper leaves things open and the code chooses: `theta_LS` of
Eq. 18 is in radians (degrees would make `b` ~ 1e3); `theta_q` of Eq. 19,
undefined in the paper, is `(q - 1) pi`; the single spin of Eq. 15 is the
mapping of arXiv:2107.08876 (`chi_par = M chi_eff / m1`,
`chi = sqrt(chi_par^2 + chi_perp^2)`, `theta = acos(chi_par / chi)`);
`fRD` is the ringdown frequency of the tuned co-precessing model.

Where they agree: `kappa_PN` of Eq. 15 term by term, the `1 + b v^5`
correction and the coefficients of Eq. 18, the plateau above `fRD`,
`A_a = kappa A_s`, the integrated phase of Eq. 24 below and above `fT`, and
the combination `h_22 = h_s + h_a`, `h*_{2,-2} = h_s - h_a` of Eqs. 6–7.

In LAL only: the IMRPhenomXP path with the antisymmetric mode
(`LALSimIMRPhenomX.c`) builds the antisymmetric phase without `phiJ_Sf`
and from the twisted-up `alpha` rather than the analytic one, which !2593
fixed only in the IMRPhenomXPHM path that IMRPhenomXPNR uses.
waveformrs has no counterpart of that path.

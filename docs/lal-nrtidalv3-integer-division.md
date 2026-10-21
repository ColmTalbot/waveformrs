# An integer division in LAL's NRTidalv3 phase derivative

`IMRPhenomX_TidalPhaseDerivative` in `LALSimIMRPhenomX_internals.c`
(lalsimulation 6.2.0 and lalsuite master as of May 2025, lines 3234–3235)
differentiates the 7.5PN tidal polynomial of NRTidalv3:

```c
REAL8 dPNpolyA = (2./3.)*pow(Mf, -1.0/3.0)*c_1A*powers_of_lalpi.two_thirds
               + c_3over2A*LAL_PI
               + (4/3)*pow(Mf, 1.0/3.0)*c_2A*powers_of_lalpi.four_thirds
               + (5./3.)*powers_of_Mf->two_thirds*c_5over2A*powers_of_lalpi.five_thirds;
```

and the same for star B. `(4/3)` is integer division in C, so it evaluates
to 1, not 4/3. The derivative of `c_2 (pi Mf)^(4/3)` is therefore too small
by a factor of 3/4. The neighbouring terms use `2./3.` and `5./3.`, so it
looks like a typo. It came in with the NRTidalv3 implementation (lalsuite
commit f1fa96d137, December 2023).

## Effect on the waveforms

None, as LAL uses the function at present. The wrong term is in the
derivative of the PN tidal phase, which enters only through the Planck
taper that switches from the NRTidal phase to the PN phase between 1.15
and 1.35 times the NRTidalv3 merger frequency. LAL evaluates the
derivative only at `f_final = min(last frequency, f_merger)`, with the same
NRTidalv3 merger frequency, to set the time shift of IMRPhenomXAS_NRTidalv3
and IMRPhenomXP_NRTidalv3 (`IMRPhenomXASGenerateFD`, `IMRPhenomXPGenerateFD`
and `IMRPhenomXHMMultiBandOneMode`). There the taper and its derivative
both vanish, so the PN term is multiplied by zero.

The term would matter if the derivative were evaluated above
`1.15 f_merger`, for example if the time shift were moved to a later
frequency.

## In waveformrs

The port keeps `(4/3)` as 1 so that it reproduces LAL exactly, with a
comment pointing here. Since the term does not change the output, fixing
it in LAL (or here) would change no waveform.

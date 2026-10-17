//! Orbit-averaged SpinTaylor evolution of precessing binaries.
//!
//! Ports of the adaptive Runge-Kutta-Fehlberg integrator of LAL
//! (`XLALAdaptiveRungeKutta4Hermite`) and GSL 2.8 (`gsl_odeiv_*`) used by
//! `XLALSimInspiralSpinTaylorPNEvolveOrbit`.

pub mod fp;
pub mod integrator;

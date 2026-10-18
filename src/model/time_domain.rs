//! [`TimeDomainModel`] for NRSur7dq4.

use crate::error::Result;
use crate::model::{
    Approximant, Binary, TimeDomainModel, TimeSeriesModes, TimeSeriesPolarizations,
};
use crate::nrsur::{NRSur7dq4, NRSur7dq4Data, NRSurWorkspace};

/// NRSur7dq4 as `SimInspiralChooseTDWaveform` / `SimInspiralChooseTDModes`.
/// [`TimeDomainModel::from_binary`] uses the process-wide surrogate data of
/// [`NRSur7dq4Data::global`] (found through `NRSUR7DQ4_DATA` or
/// `LAL_DATA_PATH`); use [`NRSur7dq4::new`] for other data.
impl TimeDomainModel for NRSur7dq4<'_> {
    type Workspace = NRSurWorkspace;

    fn from_binary(binary: &Binary) -> Result<Self> {
        Ok(NRSur7dq4::new(NRSur7dq4Data::global()?, binary)?)
    }

    fn approximant(&self) -> Approximant {
        Approximant::NRSur7dq4
    }

    fn polarizations_into(
        &self,
        workspace: &mut NRSurWorkspace,
        delta_t: f64,
        f_min: f64,
        out: &mut TimeSeriesPolarizations,
    ) -> Result<()> {
        Ok(self.polarizations_with(workspace, delta_t, f_min, out)?)
    }

    fn modes_into(
        &self,
        workspace: &mut NRSurWorkspace,
        delta_t: f64,
        f_min: f64,
        out: &mut TimeSeriesModes,
    ) -> Result<()> {
        Ok(self.modes_with(workspace, delta_t, f_min, out)?)
    }
}

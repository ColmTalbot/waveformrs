//! The surrogate data (`PrecessingNRSurData` in LALSimulation), read from the
//! HDF5 file with [`crate::hdf5`], plus tables precomputed once when loading.

use std::path::Path;

use crate::hdf5;
use crate::nrsur::fit::Fit;
use crate::nrsur::frames::WignerTables;
use crate::nrsur::spline::{FixedKnots, FixedPoints};
use crate::nrsur::NRSurError;

/// Highest `ell` of the coorbital modes (`NRSUR_LMAX`).
pub(crate) const LMAX: usize = 4;
/// Number of `(ell, m)` modes with `2 <= ell <= LMAX`.
pub(crate) const N_MODES: usize = LMAX * (LMAX + 2) - 3;

/// Index of the `(ell, m)` mode in the `MultiModalWaveform` ordering.
#[inline]
pub(crate) const fn mode_index(ell: usize, m: i32) -> usize {
    ((ell * (ell + 1)) as i32 - 4 + m) as usize
}

/// The fits at one dynamics node (`DynamicsNodeFitData`).
#[derive(Clone, Debug)]
pub(crate) struct DynamicsNode {
    /// Orbital angular frequency.
    pub(crate) omega: Fit,
    /// x, y components of `Omega^coorb`.
    pub(crate) omega_copr: [Fit; 2],
    /// Coorbital components of the time derivative of chiA in the
    /// coprecessing frame.
    pub(crate) chi_a_dot: [Fit; 3],
    /// The same for chiB.
    pub(crate) chi_b_dot: [Fit; 3],
}

/// One real component of the coorbital waveform (`WaveformDataPiece`).
#[derive(Clone, Debug)]
pub(crate) struct DataPiece {
    /// The empirical node times (indices into `t_coorb`).
    pub(crate) node_indices: Vec<usize>,
    /// The fit at each empirical node.
    pub(crate) fits: Vec<Fit>,
    /// The empirical interpolation basis, `n_nodes x n_coorb` row-major,
    /// with LAL's sign convention applied.
    pub(crate) basis: Vec<f64>,
}

/// The data pieces for fixed `ell` (`WaveformFixedEllModeData`).
#[derive(Clone, Debug)]
pub(crate) struct EllModes {
    pub(crate) m0_real: DataPiece,
    pub(crate) m0_imag: DataPiece,
    /// `Re[X_+]`, `Re[X_-]`, `Im[X_+]`, `Im[X_-]` for `m = 1..=ell`.
    pub(crate) x_real_plus: Vec<DataPiece>,
    pub(crate) x_real_minus: Vec<DataPiece>,
    pub(crate) x_imag_plus: Vec<DataPiece>,
    pub(crate) x_imag_minus: Vec<DataPiece>,
}

/// Which precessing surrogate the data belong to
/// (`PrecessingNRSurVersion`). Only NRSur7dq4 is implemented; NRSur7dq2
/// differs in the fit parametrisation, the vector fits and the limits.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SurrogateVersion {
    NRSur7dq4,
}

impl SurrogateVersion {
    /// The first dynamics time, nudged inside the first node
    /// (`NRSUR7DQ4_START_TIME`).
    pub(crate) fn start_time(self) -> f64 {
        match self {
            Self::NRSur7dq4 => -4299.99999999,
        }
    }

    /// Largest mass ratio without `unlimited_extrapolation`.
    pub fn q_max(self) -> f64 {
        match self {
            Self::NRSur7dq4 => 6.01,
        }
    }

    /// Mass ratio beyond which LAL warns about extrapolation.
    pub fn q_max_warn(self) -> f64 {
        match self {
            Self::NRSur7dq4 => 4.01,
        }
    }

    /// Largest spin magnitude without `unlimited_extrapolation`.
    pub fn chi_max(self) -> f64 {
        match self {
            Self::NRSur7dq4 => 1.0,
        }
    }

    /// Spin magnitude beyond which LAL warns about extrapolation.
    pub fn chi_max_warn(self) -> f64 {
        match self {
            Self::NRSur7dq4 => 0.81,
        }
    }
}

/// The NRSur7dq4 surrogate data, loaded once from `NRSur7dq4_v1.0.h5` and
/// shared (read-only) by any number of waveform evaluations.
///
/// Besides the fit coefficients and empirical-interpolation bases this holds
/// tables derived from them: the basis functions as indices into the power
/// table of the fit parameters, the factorised natural-spline systems for
/// the dynamics nodes and the coorbital times, the positions of the
/// coorbital times between the dynamics nodes, and the Wigner-D
/// coefficients.
#[derive(Clone, Debug)]
pub struct NRSur7dq4Data {
    pub(crate) version: SurrogateVersion,
    /// Dynamics node times `t_ds` (without the half nodes).
    pub(crate) t_ds: FixedKnots,
    /// The half nodes `t_1/2`, `t_3/2`, `t_5/2` used to start the RK4 steps.
    pub(crate) t_ds_half: [f64; 3],
    pub(crate) ds_nodes: Vec<DynamicsNode>,
    pub(crate) ds_half_nodes: Vec<DynamicsNode>,
    /// The coorbital (waveform) times `t_coorb`.
    pub(crate) t_coorb: FixedKnots,
    /// `t_coorb` located between the dynamics nodes.
    pub(crate) coorb_in_ds: FixedPoints,
    /// Data pieces for `ell = 2, 3, 4`.
    pub(crate) modes: Vec<EllModes>,
    pub(crate) wigner: WignerTables,
}

fn invalid(why: impl Into<String>) -> NRSurError {
    NRSurError::InvalidData(why.into())
}

fn read_f64(group: &hdf5::Group, name: &str) -> Result<Vec<f64>, NRSurError> {
    Ok(group.dataset(name)?.read_f64()?)
}

fn read_i64(group: &hdf5::Group, name: &str) -> Result<Vec<i64>, NRSurError> {
    Ok(group.dataset(name)?.read_i64()?)
}

fn load_fit(group: &hdf5::Group, name: &str) -> Result<Fit, NRSurError> {
    let coefs = read_f64(group, &format!("{name}_coefs"))?;
    let orders_set = group.dataset(&format!("{name}_bfOrders"))?;
    if orders_set.shape() != [coefs.len(), 7] {
        return Err(invalid(format!(
            "{}: shape {:?} for {} coefficients",
            orders_set.path(),
            orders_set.shape(),
            coefs.len()
        )));
    }
    let orders = orders_set.read_i64()?;
    Fit::new(coefs, &orders).ok_or_else(|| invalid(format!("{name}: invalid basis orders")))
}

fn load_dynamics_node(group: &hdf5::Group) -> Result<DynamicsNode, NRSurError> {
    Ok(DynamicsNode {
        omega: load_fit(group, "omega")?,
        omega_copr: [
            load_fit(group, "omega_orb_0")?,
            load_fit(group, "omega_orb_1")?,
        ],
        chi_a_dot: [
            load_fit(group, "chiA_0")?,
            load_fit(group, "chiA_1")?,
            load_fit(group, "chiA_2")?,
        ],
        chi_b_dot: [
            load_fit(group, "chiB_0")?,
            load_fit(group, "chiB_1")?,
            load_fit(group, "chiB_2")?,
        ],
    })
}

fn load_data_piece(
    group: &hdf5::Group,
    invert_sign: bool,
    n_coorb: usize,
) -> Result<DataPiece, NRSurError> {
    let basis_set = group.dataset("EIBasis")?;
    let node_indices = read_i64(group, "nodeIndices")?;
    let n_nodes = node_indices.len();
    if basis_set.shape() != [n_nodes, n_coorb] {
        return Err(invalid(format!(
            "{}: shape {:?}, expected [{n_nodes}, {n_coorb}]",
            basis_set.path(),
            basis_set.shape()
        )));
    }
    let mut basis = basis_set.read_f64()?;
    if invert_sign {
        basis.iter_mut().for_each(|b| *b *= -1.0);
    }
    let node_indices = node_indices
        .into_iter()
        .map(|i| usize::try_from(i).ok().filter(|&i| i < n_coorb))
        .collect::<Option<Vec<_>>>()
        .ok_or_else(|| invalid(format!("{}: node index out of range", group.path())))?;
    let modelers = group.group("nodeModelers")?;
    let mut fits = Vec::with_capacity(n_nodes);
    for i in 0..n_nodes {
        let coefs = read_f64(&modelers, &format!("coefs_{i}"))?;
        let orders_set = modelers.dataset(&format!("bfOrders_{i}"))?;
        if orders_set.shape() != [coefs.len(), 7] {
            return Err(invalid(format!("{}: unexpected shape", orders_set.path())));
        }
        let orders = orders_set.read_i64()?;
        fits.push(
            Fit::new(coefs, &orders)
                .ok_or_else(|| invalid(format!("{}: invalid basis orders", orders_set.path())))?,
        );
    }
    Ok(DataPiece {
        node_indices,
        fits,
        basis,
    })
}

impl NRSur7dq4Data {
    /// Load the surrogate data from `NRSur7dq4_v1.0.h5`
    /// (`PrecessingNRSur_Init`).
    pub fn load(path: impl AsRef<Path>) -> Result<Self, NRSurError> {
        let file = hdf5::File::open(path)?;
        Self::from_file(&file)
    }

    /// Load the surrogate data from an open HDF5 file.
    pub fn from_file(file: &hdf5::File) -> Result<Self, NRSurError> {
        let root = file.root()?;
        let t_with_halves = read_f64(&root, "t_ds")?;
        if t_with_halves.len() < 10 {
            return Err(invalid("too few dynamics nodes"));
        }
        let n_ds = t_with_halves.len() - 3;
        let mut t_ds = Vec::with_capacity(n_ds);
        let mut t_ds_half = [0.0; 3];
        for i in 0..3 {
            t_ds.push(t_with_halves[2 * i]);
            t_ds_half[i] = t_with_halves[2 * i + 1];
        }
        t_ds.extend_from_slice(&t_with_halves[6..]);
        if t_ds.windows(2).any(|w| w[1] <= w[0]) {
            return Err(invalid("t_ds is not increasing"));
        }

        // Nodes are stored as ds_node_0 (t_0), ds_node_1 (t_1/2), ds_node_2
        // (t_1), ds_node_3 (t_3/2), ds_node_4 (t_2), ds_node_5 (t_5/2),
        // ds_node_6 (t_3), ds_node_{n+3} (t_n) for n >= 3.
        let mut ds_nodes = Vec::with_capacity(n_ds);
        let mut ds_half_nodes = Vec::with_capacity(3);
        for i in 0..n_ds {
            let j = if i < 3 { 2 * i } else { i + 3 };
            ds_nodes.push(load_dynamics_node(&root.group(&format!("ds_node_{j}"))?)?);
            if i < 3 {
                ds_half_nodes.push(load_dynamics_node(
                    &root.group(&format!("ds_node_{}", j + 1))?,
                )?);
            }
        }

        let t_coorb = read_f64(&root, "t_coorb")?;
        let n_coorb = t_coorb.len();
        if n_coorb < 4 || t_coorb.windows(2).any(|w| w[1] <= w[0]) {
            return Err(invalid("t_coorb is not increasing"));
        }
        if t_coorb[0] < t_ds[0] || t_coorb[n_coorb - 1] > t_ds[n_ds - 1] {
            return Err(invalid("t_coorb extends beyond the dynamics nodes"));
        }

        let mut modes = Vec::with_capacity(LMAX - 1);
        for ell in 2..=LMAX {
            let piece = |name: String, invert: bool| -> Result<DataPiece, NRSurError> {
                load_data_piece(&root.group(&name)?, invert, n_coorb)
            };
            // The file stores Y_pm = (h^{l,-m} pm h^{l,m}*) / 2; X_pm = pm Y_pm*
            // gives the sign flips of Re[X_-] and Im[X_+].
            let mut x = [Vec::new(), Vec::new(), Vec::new(), Vec::new()];
            for m in 1..=ell {
                x[0].push(piece(format!("hCoorb_{ell}_{m}_Re+"), false)?);
                x[1].push(piece(format!("hCoorb_{ell}_{m}_Re-"), true)?);
                x[2].push(piece(format!("hCoorb_{ell}_{m}_Im+"), true)?);
                x[3].push(piece(format!("hCoorb_{ell}_{m}_Im-"), false)?);
            }
            let [x_real_plus, x_real_minus, x_imag_plus, x_imag_minus] = x;
            modes.push(EllModes {
                m0_real: piece(format!("hCoorb_{ell}_0_real"), false)?,
                m0_imag: piece(format!("hCoorb_{ell}_0_imag"), false)?,
                x_real_plus,
                x_real_minus,
                x_imag_plus,
                x_imag_minus,
            });
        }

        let t_ds = FixedKnots::new(&t_ds);
        let coorb_in_ds = FixedPoints::new(&t_ds, &t_coorb);
        Ok(Self {
            version: SurrogateVersion::NRSur7dq4,
            t_ds,
            t_ds_half,
            ds_nodes,
            ds_half_nodes,
            t_coorb: FixedKnots::new(&t_coorb),
            coorb_in_ds,
            modes,
            wigner: WignerTables::new(LMAX),
        })
    }

    /// The surrogate version of the data.
    pub fn version(&self) -> SurrogateVersion {
        self.version
    }

    /// The dynamics node times (in units of the total mass, with the peak
    /// of the waveform at zero).
    pub fn dynamics_times(&self) -> &[f64] {
        &self.t_ds.x
    }

    /// The half nodes `t_1/2`, `t_3/2`, `t_5/2` whose fits are used for the
    /// substeps of the first three RK4 steps.
    pub fn half_node_times(&self) -> [f64; 3] {
        self.t_ds_half
    }

    /// The times of the coorbital waveform samples.
    pub fn coorbital_times(&self) -> &[f64] {
        &self.t_coorb.x
    }

    /// The total number of fit coefficients (a measure of the cost of an
    /// evaluation).
    pub fn n_coefficients(&self) -> usize {
        let node = |n: &DynamicsNode| {
            n.omega.len()
                + n.omega_copr.iter().map(Fit::len).sum::<usize>()
                + n.chi_a_dot.iter().map(Fit::len).sum::<usize>()
                + n.chi_b_dot.iter().map(Fit::len).sum::<usize>()
        };
        let piece = |p: &DataPiece| p.fits.iter().map(Fit::len).sum::<usize>();
        self.ds_nodes
            .iter()
            .chain(&self.ds_half_nodes)
            .map(node)
            .sum::<usize>()
            + self
                .modes
                .iter()
                .map(|m| {
                    piece(&m.m0_real)
                        + piece(&m.m0_imag)
                        + m.x_real_plus
                            .iter()
                            .chain(&m.x_real_minus)
                            .chain(&m.x_imag_plus)
                            .chain(&m.x_imag_minus)
                            .map(piece)
                            .sum::<usize>()
                })
                .sum::<usize>()
    }
}

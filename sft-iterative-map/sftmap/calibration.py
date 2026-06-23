"""Calibrate the iterative-map parameters from SFT simulations of AR emergence.

This realises Objectives 1 and 2 of the proposal: it runs Monte-Carlo ensembles
of observation-constrained AR emergence (:mod:`sftmap.ar_emergence`) through the
surface-flux-transport step (:mod:`sftmap.sft`) across a grid of prescribed
cycle amplitudes ``SN``, builds the relationship ``dDM(SN)`` -- both its mean
trend and its stochastic scatter -- and fits it to the J20 / Wang et al. form

    dDM = k1 * erf( SN / quench ) * (1 + stoch * X)

to recover the saturation level ``k1``, the saturation scale ``quench`` and the
multiplicative stochasticity ``stoch``.  The polar-precursor coefficient
``k0`` is taken from observations (Jiang et al. 2018; default 58.7) and is not
calibrated here.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np
from scipy.optimize import curve_fit
from scipy.special import erf

from .ar_emergence import EmergenceParameters, generate_cycle
from .iterative_map import MapParameters
from .sft import DEFAULT_DIPOLE_NORM, DEFAULT_LAMBDA_R, SFTModel, cycle_delta_dm

__all__ = ["CalibrationResult", "simulate_delta_dm", "calibrate"]


def _erf_model(sn, k1, quench):
    return k1 * erf(sn / quench)


@dataclass
class CalibrationResult:
    """Outcome of a calibration run."""

    sn_grid: np.ndarray
    dm_mean: np.ndarray             # mean dDM at each SN
    dm_std: np.ndarray              # std of dDM at each SN
    k1: float
    quench: float
    stoch: float
    k0: float = 58.7
    constrained_min: float = 107.1  # weakest observationally constrained cycle
    samples: Optional[np.ndarray] = field(default=None, repr=False)

    def parameters(self) -> MapParameters:
        """Return the calibrated :class:`~sftmap.iterative_map.MapParameters`."""
        return MapParameters(k0=self.k0, k1=self.k1, quench=self.quench, stoch=self.stoch)

    def summary(self) -> str:
        return (
            "Calibrated iterative-map parameters\n"
            "-----------------------------------\n"
            f"  k1     (saturation level) = {self.k1:.3f}\n"
            f"  quench (saturation scale) = {self.quench:.2f}\n"
            f"  stoch  (multiplicative noise) = {self.stoch:.3f}\n"
            f"  k0     (fixed, observed)  = {self.k0:.2f}\n"
            f"  => k0*k1 = {self.k0 * self.k1:.1f}"
        )


def simulate_delta_dm(
    sn: float,
    n_realizations: int = 200,
    emergence: EmergenceParameters = EmergenceParameters(),
    lambda_R: float = DEFAULT_LAMBDA_R,
    dipole_norm: float = DEFAULT_DIPOLE_NORM,
    rng: Optional[np.random.Generator] = None,
    use_pde: bool = False,
    pde_kwargs: Optional[dict] = None,
) -> np.ndarray:
    """Monte-Carlo samples of ``dDM`` for a single prescribed cycle amplitude.

    By default the fast semi-analytic surrogate is used.  Set ``use_pde=True``
    to evolve each realisation with the full 1-D PDE solver (far slower).
    """
    if rng is None:
        rng = np.random.default_rng()
    out = np.empty(n_realizations, dtype=float)
    pde_kwargs = pde_kwargs or {}
    for i in range(n_realizations):
        ars = generate_cycle(sn, emergence, rng=rng)
        if use_pde:
            model = SFTModel(**pde_kwargs)
            out[i] = model.cycle_delta_dm(ars)
        else:
            out[i] = cycle_delta_dm(ars, lambda_R=lambda_R, dipole_norm=dipole_norm)
    return out


def calibrate(
    sn_grid: Optional[np.ndarray] = None,
    n_realizations: int = 300,
    emergence: EmergenceParameters = EmergenceParameters(),
    lambda_R: float = DEFAULT_LAMBDA_R,
    dipole_norm: float = DEFAULT_DIPOLE_NORM,
    k0: float = 58.7,
    constrained_min: float = 107.1,
    seed: Optional[int] = 0,
    use_pde: bool = False,
    keep_samples: bool = False,
) -> CalibrationResult:
    """Run the calibration pipeline and fit ``(k1, quench, stoch)``.

    Parameters
    ----------
    sn_grid : ndarray, optional
        Grid of cycle amplitudes to simulate (default: 20..360 in steps of 20).
    n_realizations : int
        Monte-Carlo realisations per grid point.
    emergence : EmergenceParameters
        Observation-constrained AR-emergence model.
    lambda_R, dipole_norm : float
        Surrogate SFT parameters.
    k0 : float
        Fixed polar-precursor coefficient.
    constrained_min : float
        Amplitudes below this are not used to estimate ``stoch`` (they are not
        observationally constrained and have noisy relative scatter).
    seed : int, optional
        Base RNG seed.
    use_pde : bool
        Use the full PDE solver instead of the surrogate.
    keep_samples : bool
        Store every ``dDM`` sample in the result (memory heavy).

    Returns
    -------
    CalibrationResult
    """
    if sn_grid is None:
        sn_grid = np.arange(20.0, 361.0, 20.0)
    sn_grid = np.asarray(sn_grid, dtype=float)

    rng = np.random.default_rng(seed)
    dm_mean = np.empty_like(sn_grid)
    dm_std = np.empty_like(sn_grid)
    all_samples = [] if keep_samples else None

    for i, sn in enumerate(sn_grid):
        samples = simulate_delta_dm(
            sn, n_realizations=n_realizations, emergence=emergence,
            lambda_R=lambda_R, dipole_norm=dipole_norm, rng=rng, use_pde=use_pde,
        )
        dm_mean[i] = np.mean(samples)
        dm_std[i] = np.std(samples)
        if keep_samples:
            all_samples.append(samples)

    # Fit the mean trend to the erf saturation law.
    p0 = (max(dm_mean.max(), 1e-3), max(sn_grid.mean(), 1.0))
    popt, _ = curve_fit(_erf_model, sn_grid, dm_mean, p0=p0, maxfev=20000)
    k1, quench = float(popt[0]), float(abs(popt[1]))

    # Estimate the multiplicative stochasticity from the relative scatter over
    # the observationally constrained range.
    mask = (sn_grid >= constrained_min) & (dm_mean > 0)
    rel = dm_std[mask] / dm_mean[mask]
    stoch = float(np.median(rel)) if rel.size else float("nan")

    return CalibrationResult(
        sn_grid=sn_grid,
        dm_mean=dm_mean,
        dm_std=dm_std,
        k1=k1,
        quench=quench,
        stoch=stoch,
        k0=k0,
        constrained_min=constrained_min,
        samples=(np.array(all_samples) if keep_samples else None),
    )

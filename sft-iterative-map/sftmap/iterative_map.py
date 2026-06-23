"""One-dimensional iterative map for solar-cycle amplitudes.

This module implements the observation-based iterative map of Wang, Jiang &
Wang (2025, ApJ 984, 183, "Paper I"; RAA 25, 125013, "Paper II"), following the
Babcock-Leighton (B-L) quantification of Jiang (2020, "J20").

The map relates the amplitude ``SN(n)`` of cycle ``n`` (maximum of the 13-month
smoothed monthly sunspot number, SSN v2) to the amplitude of cycle ``n+1``:

    SN(n+1) = k0 * k1 * erf( SN(n) / quench ) * (1 + stoch * X) - SN(n)         (1)

with a reflecting boundary at ``SN = 0`` (a negative draw is reflected back to
its absolute value so the iteration can continue).

The two physical steps behind Eq. (1) are:

* The linear (and weakly noisy) Omega-effect / polar precursor:
      SN(n+1) = k0 * DM(n),                                                     (2)
  where ``DM(n)`` is the axial dipole moment at the cycle n/n+1 minimum and
  ``k0`` is the polar-precursor coefficient.

* The nonlinear, stochastic B-L generation of poloidal field, which produces
  the *increment* of the dipole moment contributed by cycle ``n``:
      DM(n) = dDM(n) - DM(n-1)      (Hale reversal: the new field cancels the
                                     residual old field),                       (3)
      dDM   = k1 * erf( SN / quench ) * (1 + stoch * X).                        (4)

  ``k1`` is the saturated (maximum) dipole contribution of all active regions
  in a strong cycle, ``quench`` sets the saturation scale in sunspot number,
  and ``stoch`` is the standard deviation of the multiplicative noise ``X``
  (a standard normal variable).  Equation (4) is the quantity this package
  *calibrates* from surface-flux-transport simulations (see
  :mod:`sftmap.calibration`).

Combining (2)-(4) gives the closed first-order recursion, Eq. (1).
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Optional

import numpy as np
from scipy.optimize import brentq
from scipy.special import erf

__all__ = [
    "MapParameters",
    "STANDARD",
    "USOSKIN_MATCHED",
    "delta_dm",
    "recurrence",
    "generate_series",
    "fixed_point",
    "map_derivative",
    "lyapunov_exponent",
    "amplitude_pdf",
    "predict_next",
    "usoskin_average_to_max",
]


@dataclass(frozen=True)
class MapParameters:
    """Parameters of the solar-cycle iterative map, Eq. (1).

    Attributes
    ----------
    k0 : float
        Polar-precursor coefficient (linear Omega-effect), relating the axial
        dipole moment at minimum to the next cycle amplitude.
    k1 : float
        Saturated dipole-moment contribution of all active regions in a strong
        cycle (sets the saturation level of ``dDM``).
    quench : float
        Saturation scale in sunspot number (controls how fast ``dDM`` saturates
        with ``SN``).  Encodes the combined effect of tilt and latitude
        quenching.
    stoch : float
        Standard deviation of the multiplicative B-L stochasticity ``X``.
    """

    k0: float = 58.7
    k1: float = 6.94
    quench: float = 75.85
    stoch: float = 0.17

    def with_(self, **kwargs) -> "MapParameters":
        """Return a copy with selected fields replaced (e.g. ``p.with_(k1=5.0)``)."""
        return replace(self, **kwargs)


#: J20 / Wang et al. (2025) "standard" parameter set.
STANDARD = MapParameters(k0=58.7, k1=6.94, quench=75.85, stoch=0.17)

#: Illustrative set closer to the Usoskin et al. (2014) reconstructed amplitude
#: PDF (Paper I, Fig. 4d): the B-L dynamo saturates at higher amplitude (larger
#: ``quench``) but with a smaller saturated dipole (smaller ``k1``) and slightly
#: weaker stochasticity.  These values are obtained by trial, not optimisation,
#: and roughly reproduce the Paper I cycle-26 prediction of ~194 +/- 44.
USOSKIN_MATCHED = MapParameters(k0=58.7, k1=6.1, quench=100.0, stoch=0.13)


def delta_dm(
    sn: np.ndarray | float,
    p: MapParameters = STANDARD,
    rng: Optional[np.random.Generator] = None,
    stochastic: bool = True,
) -> np.ndarray | float:
    """Dipole-moment contribution ``dDM`` of a cycle of amplitude ``sn``, Eq. (4).

    Parameters
    ----------
    sn : float or ndarray
        Cycle amplitude(s).
    p : MapParameters
        Map parameters.
    rng : numpy.random.Generator, optional
        Random generator used for the stochastic term.  A fresh default
        generator is used if omitted.
    stochastic : bool
        If False, return the deterministic mean trend ``k1 * erf(sn/quench)``.

    Returns
    -------
    float or ndarray
        The (possibly stochastic) dipole contribution ``dDM``.
    """
    sn = np.asarray(sn, dtype=float)
    mean = p.k1 * erf(sn / p.quench)
    if not stochastic or p.stoch == 0.0:
        return mean
    if rng is None:
        rng = np.random.default_rng()
    x = rng.standard_normal(size=sn.shape if sn.ndim else None)
    return mean * (1.0 + p.stoch * x)


def recurrence(
    sn: float,
    p: MapParameters = STANDARD,
    rng: Optional[np.random.Generator] = None,
    stochastic: bool = True,
) -> float:
    """Advance the map one step: return ``SN(n+1)`` given ``SN(n) = sn``, Eq. (1).

    A reflecting boundary is applied at ``SN = 0`` (the absolute value is
    returned), as in Paper II.
    """
    nxt = p.k0 * delta_dm(sn, p, rng=rng, stochastic=stochastic) - sn
    return abs(float(nxt))


def generate_series(
    n: int,
    sn0: float = 150.0,
    p: MapParameters = STANDARD,
    seed: Optional[int] = None,
    stochastic: bool = True,
    burn_in: int = 0,
) -> np.ndarray:
    """Generate a series of ``n`` cycle amplitudes from the iterative map.

    Parameters
    ----------
    n : int
        Number of amplitudes to return.
    sn0 : float
        Initial cycle amplitude.
    p : MapParameters
        Map parameters.
    seed : int, optional
        Seed for reproducibility.
    stochastic : bool
        Use the stochastic map (Eq. 1) if True, else the deterministic map.
    burn_in : int
        Number of leading iterations to discard before recording (useful to
        remove dependence on the initial condition).

    Returns
    -------
    ndarray
        Array of ``n`` cycle amplitudes.
    """
    rng = np.random.default_rng(seed)
    sn = float(sn0)
    for _ in range(burn_in):
        sn = recurrence(sn, p, rng=rng, stochastic=stochastic)
    out = np.empty(n, dtype=float)
    for i in range(n):
        out[i] = sn
        sn = recurrence(sn, p, rng=rng, stochastic=stochastic)
    return out


def fixed_point(p: MapParameters = STANDARD) -> float:
    """Nonzero fixed point ``SN*`` of the deterministic map.

    Solves ``SN* = k0 * k1 * erf(SN*/quench) - SN*``, i.e.
    ``2 SN* = k0 k1 erf(SN*/quench)``.
    """
    def g(sn: float) -> float:
        return p.k0 * p.k1 * erf(sn / p.quench) - 2.0 * sn

    # The map saturates at SN = k0 k1; the fixed point lies below it.
    hi = p.k0 * p.k1
    # g(0) = 0 (trivial root); search just above 0 for the nonzero root.
    lo = 1e-6
    if g(lo) <= 0:
        # No nonzero fixed point (origin already attracting).
        return 0.0
    return float(brentq(g, lo, hi))


def map_derivative(sn: np.ndarray | float, p: MapParameters = STANDARD) -> np.ndarray | float:
    """Derivative of the deterministic recursion, Paper I, Eq. (6).

    ``d SN(n+1) / d SN(n) = k0 k1 * (2 / (quench sqrt(pi))) *
    exp(-(SN/quench)^2) - 1``.
    """
    sn = np.asarray(sn, dtype=float)
    return p.k0 * p.k1 * (2.0 / (p.quench * np.sqrt(np.pi))) * np.exp(-((sn / p.quench) ** 2)) - 1.0


def lyapunov_exponent(
    p: MapParameters = STANDARD,
    sn0: float = 150.0,
    n: int = 5000,
) -> float:
    """Largest Lyapunov exponent of the *deterministic* map.

    Computed as the average of ``log|f'(SN_n)|`` along a deterministic orbit of
    length ``n`` (Paper I, Fig. 2).  A negative value indicates convergence to
    the stable fixed point (no deterministic chaos).
    """
    sn = float(sn0)
    acc = 0.0
    count = 0
    for _ in range(n):
        d = abs(float(map_derivative(sn, p)))
        if d > 0:
            acc += np.log(d)
            count += 1
        sn = recurrence(sn, p, stochastic=False)
    return acc / count if count else float("nan")


def amplitude_pdf(
    series: np.ndarray,
    bins: int = 200,
    range_: Optional[tuple[float, float]] = None,
):
    """Probability density function of a cycle-amplitude series.

    Returns
    -------
    centers, density : ndarray, ndarray
        Bin centres and the normalised probability density.
    """
    density, edges = np.histogram(series, bins=bins, range=range_, density=True)
    centers = 0.5 * (edges[:-1] + edges[1:])
    return centers, density


def predict_next(
    sn_now: float,
    p: MapParameters = STANDARD,
    n_realizations: int = 100_000,
    seed: Optional[int] = None,
):
    """Probabilistic prediction of the next cycle amplitude.

    Draws ``n_realizations`` stochastic realisations of one map step from the
    current amplitude and returns summary statistics.

    Returns
    -------
    dict
        ``mean``, ``std`` and the 5th/50th/95th percentiles of ``SN(n+1)``.
    """
    rng = np.random.default_rng(seed)
    x = rng.standard_normal(n_realizations)
    nxt = np.abs(p.k0 * p.k1 * erf(sn_now / p.quench) * (1.0 + p.stoch * x) - sn_now)
    return {
        "mean": float(np.mean(nxt)),
        "std": float(np.std(nxt)),
        "p05": float(np.percentile(nxt, 5)),
        "p50": float(np.percentile(nxt, 50)),
        "p95": float(np.percentile(nxt, 95)),
    }


def usoskin_average_to_max(sn_average: np.ndarray | float) -> np.ndarray | float:
    """Convert a cycle-average SSN to the SSN-v2 cycle maximum convention.

    Implements Paper I, Eq. (8):
    ``SN = [1.95 * SN_ave + 24] / 0.6`` (the 1/0.6 factor converts the
    Group/Wolf sunspot number used by Usoskin et al. (2014) to SSN v2).
    """
    sn_average = np.asarray(sn_average, dtype=float)
    return (1.95 * sn_average + 24.0) / 0.6

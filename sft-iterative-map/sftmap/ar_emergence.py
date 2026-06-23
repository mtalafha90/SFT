"""Observation-constrained synthetic active-region (AR) emergence.

For a prescribed cycle amplitude ``SN`` this module generates a synthetic
population of active regions with distributions of flux, emergence latitude and
tilt angle consistent with empirical emergence laws, and with the two nonlinear
B-L feedbacks required by the proposal:

* **Tilt quenching** -- the mean Joy's-law tilt coefficient decreases with
  increasing cycle strength (Dasi-Espuig et al. 2010; Jiao et al. 2021);
* **Latitude quenching** -- the mean emergence latitude increases with cycle
  strength (Li et al. 2003; Solanki et al. 2008; Jiang 2020),

together with the **stochasticity** of AR properties (large cycle-to-cycle
scatter in tilt and latitude; Jiang et al. 2011, 2014) that is the physical
origin of the multiplicative noise in the iterative map.

The output feeds :mod:`sftmap.sft`, which converts each AR into its
contribution to the axial dipole moment at the following cycle minimum.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

__all__ = ["EmergenceParameters", "ActiveRegions", "mean_tilt_coefficient",
           "mean_latitude", "generate_cycle"]


@dataclass(frozen=True)
class EmergenceParameters:
    """Parameters of the synthetic AR-emergence model.

    Defaults are representative of cycle-19-class activity and the empirical
    laws cited above; they are deliberately exposed so the whole emergence
    model can be tuned to observations.

    Attributes
    ----------
    rate : float
        Mean number of (dipole-relevant) ARs emerging per cycle, per unit
        sunspot number.  The expected AR count for a cycle is ``rate * SN``.
    sn_ref : float
        Reference cycle amplitude at which the un-quenched laws are defined.
    lat0 : float
        Mean unsigned emergence latitude (deg) at ``SN = sn_ref``.
    q_lat : float
        Latitude-quenching slope (deg per unit SN): the mean latitude is
        ``lat0 + q_lat * (SN - sn_ref)``.
    lat_spread : float
        1-sigma scatter of the unsigned emergence latitude (deg).
    lat_min, lat_max : float
        Clamps on the unsigned emergence latitude (deg).
    joy0 : float
        Joy's-law tilt coefficient at ``SN = sn_ref`` (tilt = joy * latitude).
    q_tilt : float
        Tilt-quenching strength (per unit SN).  The Joy's-law coefficient is
        the *saturating* form ``joy0 / (1 + q_tilt * (SN - sn_ref))``, so that
        for strong cycles the tilt coefficient decreases roughly as ``1/SN``.
        Combined with the ``~SN`` growth of the AR count this produces a clean
        increase-then-saturate ``dDM(SN)`` (Dasi-Espuig et al. 2010; J20),
        rather than the increase-then-decrease of a purely linear law.
    tilt_scatter : float
        1-sigma scatter of the tilt angle about Joy's law (deg) -- the dominant
        source of B-L stochasticity.
    flux_mu, flux_sigma : float
        Log-normal parameters of the AR (unsigned) flux in units of 1e22 Mx
        (``mean = exp(flux_mu + flux_sigma**2 / 2)``).
    """

    rate: float = 1.6
    sn_ref: float = 180.0
    lat0: float = 15.0
    q_lat: float = 0.004
    lat_spread: float = 6.0
    lat_min: float = 2.0
    lat_max: float = 40.0
    joy0: float = 0.5
    q_tilt: float = 0.005
    tilt_scatter: float = 15.0
    flux_mu: float = 0.0
    flux_sigma: float = 0.5


@dataclass
class ActiveRegions:
    """A synthetic AR population for one cycle (arrays are per-AR)."""

    latitude: np.ndarray   # signed emergence latitude, deg (+N / -S)
    tilt: np.ndarray       # signed tilt angle, deg (Joy's law convention)
    flux: np.ndarray       # unsigned flux, 1e22 Mx
    sn: float              # the prescribed cycle amplitude

    def __len__(self) -> int:
        return int(self.latitude.size)


def mean_tilt_coefficient(sn: float, p: EmergenceParameters) -> float:
    """Tilt-quenched Joy's-law coefficient for a cycle of amplitude ``sn``.

    Saturating (hyperbolic) form ``joy0 / (1 + q_tilt (sn - sn_ref))``, clamped
    so the denominator stays positive and the coefficient bounded.
    """
    denom = max(1.0 + p.q_tilt * (sn - p.sn_ref), 0.2)
    return float(np.clip(p.joy0 / denom, 0.05 * p.joy0, 3.0 * p.joy0))


def mean_latitude(sn: float, p: EmergenceParameters) -> float:
    """Latitude-quenched mean (unsigned) emergence latitude, deg."""
    return float(np.clip(p.lat0 + p.q_lat * (sn - p.sn_ref), p.lat_min, p.lat_max))


def generate_cycle(
    sn: float,
    p: EmergenceParameters = EmergenceParameters(),
    rng: Optional[np.random.Generator] = None,
) -> ActiveRegions:
    """Generate a synthetic AR population for a cycle of amplitude ``sn``.

    The number of ARs is Poisson-distributed with mean ``rate * sn``.  Each AR
    is assigned a hemisphere, an unsigned emergence latitude (Gaussian about the
    latitude-quenched mean), a tilt (Joy's law with tilt quenching plus Gaussian
    scatter), and a log-normal flux.  Joy's-law sign convention: the tilt is
    such that the leading polarity lies equatorward of the following polarity,
    so the per-hemisphere dipole contributions add coherently.
    """
    if rng is None:
        rng = np.random.default_rng()

    n = rng.poisson(p.rate * max(sn, 0.0))
    if n == 0:
        empty = np.array([], dtype=float)
        return ActiveRegions(empty, empty.copy(), empty.copy(), float(sn))

    hemi = rng.choice((-1.0, 1.0), size=n)

    abslat = rng.normal(mean_latitude(sn, p), p.lat_spread, size=n)
    abslat = np.clip(abslat, p.lat_min, p.lat_max)
    latitude = hemi * abslat

    joy = mean_tilt_coefficient(sn, p)
    mean_tilt = joy * abslat                       # Joy's law (magnitude)
    tilt_mag = mean_tilt + rng.normal(0.0, p.tilt_scatter, size=n)
    tilt = hemi * tilt_mag                          # signed by hemisphere

    flux = rng.lognormal(p.flux_mu, p.flux_sigma, size=n)

    return ActiveRegions(latitude=latitude, tilt=tilt, flux=flux, sn=float(sn))

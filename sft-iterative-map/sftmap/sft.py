"""Surface flux transport (SFT): from active regions to the axial dipole moment.

Two routes from a synthetic AR population (:mod:`sftmap.ar_emergence`) to the
cycle's contribution to the axial dipole moment ``dDM`` are provided:

1. A fast **semi-analytic surrogate** (:func:`bmr_dipole_contribution`,
   :func:`cycle_delta_dm`).  The asymptotic contribution of a single bipolar
   magnetic region (BMR) to the global axial dipole, after surface flux
   transport has run to completion, is well approximated by

       dD_inf  ~  c * Phi * sin(tilt) * cos(lat) * exp(-lat^2 / (2 lambda_R^2)),

   where ``lambda_R`` is the characteristic "dynamo effective latitude" set by
   the competition between the meridional flow and supergranular diffusion
   (Petrovay et al. 2020; Iijima et al. 2017; Wang & Sheeley 1991; Jiang 2020).
   Low-latitude, high-tilt regions dominate the net dipole; high-latitude
   regions are strongly cancelled (the Gaussian factor).  This makes the
   Monte-Carlo ensembles required for calibration tractable.

2. A **full 1-D PDE solver** (:class:`SFTModel`) that integrates the
   axisymmetric flux-transport equation (meridional advection, supergranular
   diffusion, optional radial-decay term) on a colatitude grid, refactored and
   modernised from the original ``transp.py`` of this project.  Use it to
   validate the surrogate or to evolve an individual cycle with full physics.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

from .ar_emergence import ActiveRegions

__all__ = [
    "bmr_dipole_contribution",
    "cycle_delta_dm",
    "MeridionalFlow",
    "SFTModel",
]

#: Characteristic dipole (dynamo effective) latitude in degrees (Petrovay et
#: al. 2020): low-latitude regions contribute to the global dipole, higher
#: latitude ones are increasingly cancelled.
DEFAULT_LAMBDA_R = 22.0

#: Overall normalisation of the surrogate per-BMR dipole contribution.  Its
#: absolute value is convention-dependent (it absorbs the flux unit and BMR
#: polarity separation); the default is chosen so the saturated ``dDM`` of the
#: default emergence model matches the J20 value (k1 ~ 6.94).
DEFAULT_DIPOLE_NORM = 0.24


def bmr_dipole_contribution(
    latitude: np.ndarray,
    tilt: np.ndarray,
    flux: np.ndarray,
    lambda_R: float = DEFAULT_LAMBDA_R,
    dipole_norm: float = DEFAULT_DIPOLE_NORM,
) -> np.ndarray:
    """Asymptotic axial-dipole contribution of each BMR (surrogate).

    Parameters
    ----------
    latitude, tilt : ndarray
        Signed emergence latitude and tilt of each AR, in degrees.
    flux : ndarray
        Unsigned AR flux (arbitrary units; absorbed into ``dipole_norm``).
    lambda_R : float
        Dynamo effective latitude (deg).
    dipole_norm : float
        Overall normalisation constant.

    Returns
    -------
    ndarray
        Signed per-AR contribution to the axial dipole moment.  The sign is set
        by the tilt and hemisphere so that Joy's-law regions add coherently and
        anti-Joy (large negative-tilt fluctuation) regions subtract.
    """
    lat_rad = np.deg2rad(latitude)
    tilt_rad = np.deg2rad(tilt)
    geom = np.cos(lat_rad) * np.exp(-(latitude ** 2) / (2.0 * lambda_R ** 2))
    return dipole_norm * flux * np.sin(tilt_rad) * np.sign(lat_rad) * geom


def cycle_delta_dm(
    ars: ActiveRegions,
    lambda_R: float = DEFAULT_LAMBDA_R,
    dipole_norm: float = DEFAULT_DIPOLE_NORM,
) -> float:
    """Net dipole-moment contribution ``dDM`` of a whole AR population."""
    if len(ars) == 0:
        return 0.0
    contrib = bmr_dipole_contribution(
        ars.latitude, ars.tilt, ars.flux, lambda_R=lambda_R, dipole_norm=dipole_norm
    )
    return float(np.sum(contrib))


# --------------------------------------------------------------------------- #
#  Full 1-D PDE surface-flux-transport solver (modernised from transp.py)     #
# --------------------------------------------------------------------------- #

@dataclass
class MeridionalFlow:
    """Meridional-flow profile ``u(latitude)`` in m/s.

    Default is the van Ballegooijen (1998) / Jiang et al. (2014) profile
    ``u0 * sin(2.4 * lat)`` poleward of which the flow is set to zero.
    """

    u0: float = 11.0           # peak speed, m/s
    cutoff_deg: float = 75.0   # flow set to zero poleward of this latitude
    factor: float = 2.4        # argument scaling of the sine profile

    def __call__(self, latitude_deg: np.ndarray) -> np.ndarray:
        u = self.u0 * np.sin(self.factor * np.deg2rad(latitude_deg))
        return np.where(np.abs(latitude_deg) > self.cutoff_deg, 0.0, u)


class SFTModel:
    """Axisymmetric 1-D surface flux transport on a colatitude grid.

    Solves, in units of Mm / day / Gauss,

        dB/dt = -1/(R sin(t)) d/dt[ sin(t) u B ]
                + eta/(R^2 sin(t)) d/dt[ sin(t) dB/dt ]
                - B / tau + S,

    with ``t`` the colatitude, ``u`` the meridional flow, ``eta`` the
    supergranular diffusivity and ``tau`` an optional radial-decay timescale.
    The conservative annular-flux variable ``W = R sin(t) B`` is stepped with
    centred differences, as in the original ``transp.py``.
    """

    R_SUN = 695.7  # Mm

    def __init__(
        self,
        n_grid: int = 181,
        eta: float = 250.0,                 # km^2/s
        tau_years: Optional[float] = None,  # radial decay; None disables it
        flow: Optional[MeridionalFlow] = None,
    ):
        self.N = n_grid
        self.theta = np.linspace(0.0, np.pi, n_grid)
        self.dx = np.pi / (n_grid - 1)
        self.latitude = 90.0 - self.theta * 180.0 / np.pi
        # convert to Mm/day: 1 km^2/s = 8.64e-2 Mm^2/day ; 1 m/s = 8.64e-2 Mm/day
        self.eta = eta * 8.64e-2
        self.tau = None if tau_years is None else tau_years * 365.25
        self.flow = flow if flow is not None else MeridionalFlow()
        self.uc = self.flow(self.latitude) * 8.64e-2  # Mm/day
        self.B = np.zeros(n_grid)

    # -- diagnostics -------------------------------------------------------- #
    def axial_dipole_moment(self, B: Optional[np.ndarray] = None) -> float:
        """Axial dipole moment ``(3/4) * integral B sin(2 theta) d theta``."""
        B = self.B if B is None else B
        return float(0.75 * np.trapezoid(B * np.sin(2.0 * self.theta), self.theta))

    # -- sources ------------------------------------------------------------ #
    def add_bipole(self, lat0_deg: float, tilt_deg: float, amplitude: float,
                   sep_deg: float = 6.0, fwhm_deg: float = 4.0) -> None:
        """Add a bipolar magnetic region centred at ``lat0_deg`` to the field.

        The two polarities are separated in latitude by ``sep_deg * sin(tilt)``;
        ``amplitude`` scales the (signed) peak field of each polarity.
        """
        dlat = sep_deg * np.sin(np.deg2rad(tilt_deg))
        sigma = fwhm_deg / 2.3548
        lead = amplitude * np.exp(-((self.latitude - (lat0_deg - dlat / 2)) ** 2) / (2 * sigma ** 2))
        follow = -amplitude * np.exp(-((self.latitude - (lat0_deg + dlat / 2)) ** 2) / (2 * sigma ** 2))
        self.B += lead + follow

    # -- time stepping ------------------------------------------------------ #
    def step(self, dt: float, source: Optional[np.ndarray] = None) -> None:
        """Advance the field by ``dt`` days (one explicit FTCS step)."""
        x, dx, R = self.theta, self.dx, self.R_SUN
        W = R * np.sin(x) * self.B
        # advection (centred)
        Wfladv = (W / R) * self.uc
        dWflux = (np.roll(Wfladv, -1) - np.roll(Wfladv, 1)) / 2.0
        # diffusion (centred)
        diff_r = np.sin(x + dx / 2) * (np.roll(self.B, -1) - self.B) / dx
        diff_l = np.sin(x - dx / 2) * (self.B - np.roll(self.B, 1)) / dx
        dWflux += self.eta / R * (diff_r - diff_l)
        dW = dWflux / dx
        if self.tau is not None:
            dW -= W / self.tau
        if source is not None:
            dW += source * R * np.sin(x)
        W = W + dW * dt
        # recover B, regularise the poles (assume third derivative ~ 0)
        self.B[1:self.N - 1] = W[1:self.N - 1] / R / np.sin(x[1:self.N - 1])
        self.B[0] = self.B[2] + 0.5 * (self.B[1] - self.B[3])
        self.B[self.N - 1] = self.B[self.N - 3] + 0.5 * (self.B[self.N - 2] - self.B[self.N - 4])

    def run(self, n_steps: int, dt: float = 1.0) -> float:
        """Integrate freely (no source) for ``n_steps`` and return the dipole."""
        for _ in range(n_steps):
            self.step(dt)
        return self.axial_dipole_moment()

    def cycle_delta_dm(self, ars: ActiveRegions, years: float = 11.0,
                       dt: float = 1.0, amp_per_flux: float = 1.0) -> float:
        """Evolve one cycle of emergence with full physics and return ``dDM``.

        Active regions are deposited as bipoles spread uniformly over the cycle;
        the field is relaxed for one extra cycle so the dipole reaches its
        asymptotic value.  Returns the change in axial dipole moment.
        """
        dm0 = self.axial_dipole_moment()
        n_steps = int(years * 365.25 / dt)
        if len(ars) > 0:
            # Hale's law: the leading-polarity sign is opposite in the two
            # hemispheres, so the per-hemisphere dipole contributions add
            # coherently.  Encode it via the sign of the deposited amplitude
            # (the overall sign is the convention that makes dDM > 0 for a
            # Joy's-law population, matching the surrogate).
            hale = -np.sign(ars.latitude)
            emerge_step = np.random.default_rng(0).integers(0, n_steps, size=len(ars))
            order = np.argsort(emerge_step)
            idx = 0
            for s in range(n_steps):
                while idx < len(ars) and emerge_step[order[idx]] == s:
                    j = order[idx]
                    self.add_bipole(ars.latitude[j], ars.tilt[j],
                                    amp_per_flux * ars.flux[j] * hale[j])
                    idx += 1
                self.step(dt)
        else:
            for _ in range(n_steps):
                self.step(dt)
        # relax to the asymptotic dipole
        for _ in range(n_steps):
            self.step(dt)
        return self.axial_dipole_moment() - dm0

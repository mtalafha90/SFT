"""sftmap -- Observation-constrained iterative-map calibration via SFT.

A reusable Python workflow that connects the parameters of the one-dimensional
solar-cycle iterative map (Wang, Jiang & Wang 2025, Papers I & II) to the
physics of the Babcock-Leighton dynamo, through surface-flux-transport (SFT)
simulations of observation-constrained active-region emergence.

Pipeline
--------
``ar_emergence``  synthetic AR populations with tilt + latitude quenching
``sft``           AR -> axial-dipole contribution (surrogate or full 1-D PDE)
``calibration``   Monte-Carlo ensembles -> fit ``dDM(SN)`` -> (k1, quench, stoch)
``iterative_map`` the recursion SN(n+1)=k0 k1 erf(SN/quench)(1+stoch X)-SN(n)
``go_rule``       Gnevyshev-Ohl statistics
``validation``    amplitude PDF, grand maxima, prediction, G-O benchmarks
"""

from __future__ import annotations

from . import ar_emergence, calibration, go_rule, iterative_map, sft, validation
from .iterative_map import (
    MapParameters,
    STANDARD,
    USOSKIN_MATCHED,
    generate_series,
    predict_next,
)

__version__ = "0.1.0"

__all__ = [
    "ar_emergence",
    "sft",
    "calibration",
    "iterative_map",
    "go_rule",
    "validation",
    "MapParameters",
    "STANDARD",
    "USOSKIN_MATCHED",
    "generate_series",
    "predict_next",
]

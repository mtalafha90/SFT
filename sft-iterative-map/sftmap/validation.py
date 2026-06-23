"""Statistical validation of the calibrated map against observed signatures.

Implements Objective 3 of the proposal: reproduce the observed amplitude
probability density function (PDF) and the Gnevyshev-Ohl (G-O) statistics, and
bundle them with the key benchmark numbers reported by Wang et al. (2025).
"""

from __future__ import annotations

import numpy as np

from . import go_rule
from .iterative_map import (
    MapParameters,
    STANDARD,
    amplitude_pdf,
    fixed_point,
    generate_series,
    lyapunov_exponent,
    predict_next,
)

__all__ = ["grand_maxima_fraction", "peak_amplitude", "validate"]


def grand_maxima_fraction(
    series: np.ndarray, threshold: float = 258.0, block: int = 1000
) -> dict:
    """Fraction of cycles in the grand-maxima tail (amplitude > ``threshold``).

    Computed per block of ``block`` cycles; the mean and standard deviation over
    blocks are returned (Paper I reports 0.25 +/- 0.01 for the standard set,
    0.037 +/- 0.005 for the Usoskin-matched set).
    """
    series = np.asarray(series, dtype=float)
    n_blocks = len(series) // block
    if n_blocks == 0:
        return {"mean": float(np.mean(series > threshold)), "std": float("nan")}
    trimmed = series[: n_blocks * block].reshape(n_blocks, block)
    frac = np.mean(trimmed > threshold, axis=1)
    return {"mean": float(np.mean(frac)), "std": float(np.std(frac))}


def peak_amplitude(series: np.ndarray, bins: int = 200) -> float:
    """Most probable cycle amplitude (peak of the PDF)."""
    centers, density = amplitude_pdf(series, bins=bins)
    return float(centers[int(np.argmax(density))])


def validate(
    p: MapParameters = STANDARD,
    n: int = 1_000_000,
    seed: int = 0,
    sn_now: float = 150.0,
) -> dict:
    """Compute the full set of validation diagnostics for a parameter set.

    Returns a dictionary with the amplitude-PDF peak, grand-maxima fraction,
    deterministic fixed point and Lyapunov exponent, the next-cycle prediction
    from ``sn_now``, and the G-O statistics.
    """
    series = generate_series(n, sn0=sn_now, p=p, seed=seed, stochastic=True, burn_in=100)
    return {
        "parameters": p,
        "fixed_point": fixed_point(p),
        "lyapunov_exponent": lyapunov_exponent(p),
        "pdf_peak_amplitude": peak_amplitude(series),
        "grand_maxima_fraction": grand_maxima_fraction(series),
        "prediction": predict_next(sn_now, p, seed=seed),
        "go": go_rule.go_statistics(series),
        "series": series,
    }

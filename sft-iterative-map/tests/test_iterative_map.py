"""Tests for the iterative map against Wang et al. (2025) benchmarks."""

import numpy as np

from sftmap import iterative_map as im


def test_no_deterministic_chaos():
    # The deterministic map has a stable nonzero fixed point (|f'|<1) and a
    # negative Lyapunov exponent -- no chaos (Paper I).
    p = im.STANDARD
    sn_star = im.fixed_point(p)
    assert 150 < sn_star < 260
    assert abs(im.map_derivative(sn_star, p)) < 1.0
    assert im.lyapunov_exponent(p) < 0.0


def test_deterministic_converges_to_fixed_point():
    p = im.STANDARD
    series = im.generate_series(2000, sn0=50.0, p=p, stochastic=False)
    assert abs(series[-1] - im.fixed_point(p)) < 1.0


def test_prediction_standard_set():
    # Paper I: SN=150 predicts cycle 26 at 255 +/- 69 (standard set).
    r = im.predict_next(150.0, im.STANDARD, n_realizations=400_000, seed=0)
    assert abs(r["mean"] - 255) < 6
    assert abs(r["std"] - 69) < 6


def test_prediction_usoskin_set():
    # Paper I: the Usoskin-matched set predicts ~194 +/- 44.
    r = im.predict_next(150.0, im.USOSKIN_MATCHED, n_realizations=400_000, seed=0)
    assert abs(r["mean"] - 194) < 12
    assert abs(r["std"] - 44) < 8


def test_reflecting_boundary_keeps_positive():
    p = im.STANDARD.with_(stoch=0.5)  # strong noise -> would go negative
    series = im.generate_series(10_000, sn0=150.0, p=p, seed=1)
    assert np.all(series >= 0.0)


def test_pdf_peak_near_200():
    series = im.generate_series(500_000, sn0=150.0, p=im.STANDARD, seed=0, burn_in=100)
    centers, density = im.amplitude_pdf(series, bins=200, range_=(0, 500))
    peak = centers[np.argmax(density)]
    assert 180 < peak < 240

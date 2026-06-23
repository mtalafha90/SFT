"""Tests for the Gnevyshev-Ohl statistics (Wang et al. 2025, Paper II)."""

import numpy as np

from sftmap import go_rule, iterative_map as im


def _series(n=500_000, seed=0):
    return im.generate_series(n, sn0=150.0, p=im.STANDARD, seed=seed, burn_in=100)


def test_e_to_a_ratio():
    # Paper II: 0.4555 (the following cycle is more often the stronger one).
    assert abs(go_rule.e_to_a_ratio(_series()) - 0.4555) < 0.01


def test_delta_sn_distribution():
    # Paper II: mean ~ 0, median ~ 19, std ~ 157.
    stats = go_rule.delta_sn_stats(_series())
    assert abs(stats["mean"]) < 5
    assert 10 < stats["median"] < 28
    assert 130 < stats["std"] < 185


def test_correlations_negative_and_equal():
    # Paper II: E-O and O-E correlations are both about -0.42.
    c = go_rule.eo_oe_correlations(_series())
    assert abs(c["EO"] - (-0.42)) < 0.05
    assert abs(c["EO"] - c["OE"]) < 0.05


def test_block_lengths_exponential():
    # Block lengths are positive integers and skew short (exponential).
    s = _series()
    go = go_rule.go_block_lengths(s)
    alt = go_rule.alternation_block_lengths(s)
    assert go.min() >= 1 and alt.min() >= 1
    assert np.median(go) <= np.mean(go)  # right-skewed

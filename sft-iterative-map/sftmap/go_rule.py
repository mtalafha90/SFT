"""Gnevyshev-Ohl (G-O) rule statistics from a cycle-amplitude series.

Implements the diagnostics of Wang, Jiang & Wang (2025, RAA 25, 125013,
"Paper II"):

* the **E-to-A ratio** (cycle-strength definition): the fraction of G-O pairs
  in which the leading (even) cycle is stronger than the following (odd) one;
* the distribution of **DeltaSN = SN(2n+1) - SN(2n)** within pairs;
* the **E-O / O-E correlations** (correlation definition);
* the length distribution of **G-O blocks** (runs of same-sign pairs) and of
  **cycle-alternation blocks** (runs of the zig-zag pattern), both of which are
  exponential for a memoryless stochastic map.

For the standard parameter set and 1,000,000 cycles, Paper II reports an E-to-A
ratio of 0.4555 +/- 0.0003 (so the *following* cycle is more often the stronger
one), a DeltaSN distribution with mean ~ 0 but median ~ 19, and E-O/O-E
correlations of about -0.42.
"""

from __future__ import annotations

import numpy as np

__all__ = [
    "e_to_a_ratio",
    "delta_sn",
    "delta_sn_stats",
    "eo_oe_correlations",
    "go_block_lengths",
    "alternation_block_lengths",
    "go_statistics",
]


def _pairs(series: np.ndarray, start: int = 0):
    """Return (even, odd) arrays for G-O pairs (2n, 2n+1), pairing from ``start``."""
    s = np.asarray(series, dtype=float)[start:]
    npair = len(s) // 2
    s = s[: 2 * npair]
    return s[0::2], s[1::2]


def e_to_a_ratio(series: np.ndarray, start: int = 0) -> float:
    """Fraction of G-O pairs whose leading (even) cycle is the stronger one.

    A value below 0.5 means the *following* cycle is more often stronger
    (the sense of the G-O rule found by Paper II).
    """
    even, odd = _pairs(series, start)
    if len(even) == 0:
        return float("nan")
    return float(np.mean(even > odd))


def delta_sn(series: np.ndarray, start: int = 0) -> np.ndarray:
    """Within-pair amplitude difference ``DeltaSN = SN(2n+1) - SN(2n)``."""
    even, odd = _pairs(series, start)
    return odd - even


def delta_sn_stats(series: np.ndarray, start: int = 0) -> dict:
    """Mean, median and standard deviation of ``DeltaSN``."""
    d = delta_sn(series, start)
    return {
        "mean": float(np.mean(d)),
        "median": float(np.median(d)),
        "std": float(np.std(d)),
    }


def eo_oe_correlations(series: np.ndarray) -> dict:
    """E-O and O-E Pearson correlations.

    ``EO`` correlates each even cycle with the *following* odd cycle; ``OE``
    correlates each odd cycle with the *following* even cycle.  For the
    memoryless map both equal the generic lag-1 correlation.
    """
    even, odd = _pairs(series, start=0)
    n = min(len(even), len(odd))
    eo = float(np.corrcoef(even[:n], odd[:n])[0, 1]) if n > 1 else float("nan")

    even2, odd2 = _pairs(series, start=1)  # pairs (1-2), (3-4), ...
    n2 = min(len(even2), len(odd2))
    # here "even2" holds odd-indexed cycles and "odd2" the following even ones
    oe = float(np.corrcoef(even2[:n2], odd2[:n2])[0, 1]) if n2 > 1 else float("nan")
    return {"EO": eo, "OE": oe}


def _run_lengths(boolean_pattern: np.ndarray) -> np.ndarray:
    """Lengths of maximal runs of identical values in a boolean/int array."""
    a = np.asarray(boolean_pattern)
    if a.size == 0:
        return np.array([], dtype=int)
    change = np.nonzero(np.diff(a.astype(int)))[0]
    bounds = np.concatenate(([-1], change, [a.size - 1]))
    return np.diff(bounds)


def go_block_lengths(series: np.ndarray, start: int = 0) -> np.ndarray:
    """Lengths (in pairs) of G-O blocks.

    A G-O block is a maximal run of consecutive pairs that all have the same
    sign of ``SN(even) - SN(odd)``.  Block lengths are exponentially
    distributed for the stochastic map (Paper II, Fig. 3a).
    """
    even, odd = _pairs(series, start)
    sign = even > odd
    return _run_lengths(sign)


def alternation_block_lengths(series: np.ndarray) -> np.ndarray:
    """Lengths of cycle-alternation (zig-zag) blocks.

    An alternation block is a maximal run over which the sign of
    ``SN(n+1) - SN(n)`` keeps flipping (high-low-high-...).  Identified directly
    from adjacent cycles, as in Paper II (the map has no memory beyond one
    cycle).  Returns the lengths, in cycles, of the alternating runs.
    """
    s = np.asarray(series, dtype=float)
    if s.size < 2:
        return np.array([], dtype=int)
    up = np.diff(s) > 0  # True where the cycle rises
    if up.size < 2:
        return np.array([up.size + 1], dtype=int)
    # The zig-zag is broken wherever two consecutive steps have the *same* sign.
    broken = up[1:] == up[:-1]
    # Each "broken" point ends an alternation block.  Convert the boolean
    # break-pattern into run lengths of the difference series, then +1 to count
    # cycles rather than steps.
    bounds = np.concatenate(([-1], np.nonzero(broken)[0], [up.size - 1]))
    step_runs = np.diff(bounds)
    return step_runs + 1


def go_statistics(series: np.ndarray) -> dict:
    """Convenience bundle of the main G-O diagnostics for a series."""
    return {
        "e_to_a_ratio": e_to_a_ratio(series),
        "delta_sn": delta_sn_stats(series),
        "correlations": eo_oe_correlations(series),
        "mean_go_block": float(np.mean(go_block_lengths(series))),
        "mean_alternation_block": float(np.mean(alternation_block_lengths(series))),
    }

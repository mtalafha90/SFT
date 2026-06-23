#!/usr/bin/env python3
"""Analyse the iterative map: cobweb diagram, amplitude series and PDF.

Reproduces the core figures of Wang et al. (2025, Paper I): the recursion
function with cobweb iterations (deterministic and stochastic), a synthesized
amplitude series, and the amplitude probability density function.  Also prints
the fixed point, Lyapunov exponent and a next-cycle prediction.

Example
-------
    python scripts/run_map_analysis.py --param standard --out map_analysis.png
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from sftmap import iterative_map as im
from sftmap import validation


def cobweb_points(series, p):
    """Return staircase (x, y) arrays tracing the cobweb iteration."""
    xs, ys = [], []
    for i in range(len(series) - 1):
        xs += [series[i], series[i]]
        ys += [series[i], series[i + 1]]
    return np.array(xs), np.array(ys)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--param", choices=["standard", "usoskin"], default="standard")
    ap.add_argument("--n", type=int, default=500_000, help="cycles for the PDF")
    ap.add_argument("--sn0", type=float, default=150.0, help="initial amplitude")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default="map_analysis.png")
    args = ap.parse_args()

    p = im.STANDARD if args.param == "standard" else im.USOSKIN_MATCHED

    print(f"Parameter set: {args.param}  {p}")
    print(f"  fixed point        : {im.fixed_point(p):.2f}")
    print(f"  Lyapunov exponent  : {im.lyapunov_exponent(p):.4f} (<0 => no chaos)")
    pred = im.predict_next(args.sn0, p, seed=args.seed)
    print(f"  predict {args.sn0:g} -> {pred['mean']:.1f} +/- {pred['std']:.1f} "
          f"(90% CI [{pred['p05']:.0f}, {pred['p95']:.0f}])")
    gm = validation.grand_maxima_fraction(
        im.generate_series(args.n, args.sn0, p, seed=args.seed, burn_in=100))
    print(f"  grand-maxima frac  : {gm['mean']:.3f} +/- {gm['std']:.3f}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.5))

    # (a) recursion function + cobweb (deterministic, 3 initial conditions)
    grid = np.linspace(0, 2 * im.fixed_point(p) + 100, 400)
    rec = np.array([im.recurrence(g, p, stochastic=False) for g in grid])
    axes[0].plot(grid, rec, "k", lw=2, label="recursion (deterministic)")
    axes[0].plot(grid, grid, "k--", lw=0.8)
    for sn0, c in zip((50, 150, 250), ("C2", "C3", "C0")):
        det = im.generate_series(40, sn0, p, stochastic=False)
        cx, cy = cobweb_points(det, p)
        axes[0].plot(cx, cy, c, lw=0.8, alpha=0.8, label=f"$S_N(0)$={sn0}")
    axes[0].set(xlabel="$S_N(n)$", ylabel="$S_N(n+1)$", title="(a) Cobweb diagram")
    axes[0].legend(fontsize=8)

    # (b) stochastic amplitude series
    series = im.generate_series(60, args.sn0, p, seed=args.seed)
    axes[1].plot(np.arange(len(series)), series, "-o", ms=3, color="C0")
    axes[1].axhline(im.fixed_point(p), ls="--", color="grey", lw=1, label="fixed point")
    axes[1].set(xlabel="cycle number", ylabel="$S_N$",
                title="(b) Synthesized amplitude series")
    axes[1].legend(fontsize=8)

    # (c) amplitude PDF
    big = im.generate_series(args.n, args.sn0, p, seed=args.seed, burn_in=100)
    centers, density = im.amplitude_pdf(big, bins=200, range_=(0, 500))
    axes[2].plot(centers, density, color="C3", lw=2)
    axes[2].axvline(107.1, ls="--", color="grey", lw=1, label="constrained limit")
    axes[2].set(xlabel="cycle amplitude $S_N$", ylabel="probability density",
                title=f"(c) Amplitude PDF (peak={validation.peak_amplitude(big):.0f})")
    axes[2].legend(fontsize=8)

    fig.tight_layout()
    fig.savefig(args.out, dpi=150)
    print(f"\nSaved figure to {args.out}")


if __name__ == "__main__":
    main()

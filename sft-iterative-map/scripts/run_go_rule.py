#!/usr/bin/env python3
"""Gnevyshev-Ohl rule statistics from the iterative map (Wang et al. Paper II).

Generates a long amplitude series and computes the E-to-A ratio, the DeltaSN
distribution, the E-O/O-E correlations, and the exponential block-length
distributions for G-O blocks and cycle-alternation blocks.  Saves a figure.

Example
-------
    python scripts/run_go_rule.py --n 1000000 --out go_rule.png
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from sftmap import go_rule, iterative_map as im


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--param", choices=["standard", "usoskin"], default="standard")
    ap.add_argument("--n", type=int, default=1_000_000)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default="go_rule.png")
    args = ap.parse_args()

    p = im.STANDARD if args.param == "standard" else im.USOSKIN_MATCHED
    series = im.generate_series(args.n, sn0=150.0, p=p, seed=args.seed, burn_in=100)

    stats = go_rule.go_statistics(series)
    dsn = go_rule.delta_sn_stats(series)
    print(f"Parameter set: {args.param}")
    print(f"  E-to-A ratio          : {stats['e_to_a_ratio']:.4f}  (Paper II: 0.4555)")
    print(f"  DeltaSN  mean/median/std : {dsn['mean']:.1f} / {dsn['median']:.1f} / {dsn['std']:.1f}"
          f"   (Paper II: ~0 / 19 / 157)")
    print(f"  E-O / O-E correlations: {stats['correlations']['EO']:.3f} / "
          f"{stats['correlations']['OE']:.3f}  (Paper II: ~ -0.42)")
    print(f"  mean G-O block length : {stats['mean_go_block']:.2f} pairs")
    print(f"  mean alternation block: {stats['mean_alternation_block']:.2f} cycles")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.5))

    d = go_rule.delta_sn(series)
    axes[0].hist(d, bins=200, range=(-600, 600), density=True, color="C0")
    axes[0].axvline(0, color="k", lw=0.8)
    axes[0].axvline(dsn["median"], color="r", ls="--", lw=1,
                    label=f"median={dsn['median']:.0f}")
    axes[0].set(xlabel=r"$\Delta S_N = S_N(2n{+}1)-S_N(2n)$", ylabel="density",
                title="(a) Within-pair difference")
    axes[0].legend(fontsize=8)

    for ax, lengths, title in (
        (axes[1], go_rule.go_block_lengths(series), "(b) G-O block lengths"),
        (axes[2], go_rule.alternation_block_lengths(series), "(c) Alternation block lengths"),
    ):
        maxlen = int(lengths.max())
        bins = np.arange(0.5, maxlen + 1.5)
        ax.hist(lengths, bins=bins, density=True, color="C2", alpha=0.8)
        ax.set_yscale("log")
        ax.set(xlabel="block length", ylabel="probability",
               title=f"{title} (exponential)")

    fig.tight_layout()
    fig.savefig(args.out, dpi=150)
    print(f"\nSaved figure to {args.out}")


if __name__ == "__main__":
    main()

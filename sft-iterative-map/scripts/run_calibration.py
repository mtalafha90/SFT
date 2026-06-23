#!/usr/bin/env python3
"""Calibrate the iterative-map parameters from SFT simulations of AR emergence.

Runs Monte-Carlo ensembles of observation-constrained active-region emergence
through the surface-flux-transport step over a grid of cycle amplitudes, fits
the resulting dDM(SN) relation to the erf saturation law, and reports the
calibrated (k1, quench, stoch).  Saves a figure of the dDM(SN) relation with the
erf fit and the 1-sigma scatter band.

Examples
--------
    python scripts/run_calibration.py --realizations 500 --out calibration.png
    python scripts/run_calibration.py --use-pde --realizations 30
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from sftmap import calibration as cal
from sftmap.iterative_map import STANDARD


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--realizations", type=int, default=300,
                    help="Monte-Carlo realisations per grid point")
    ap.add_argument("--sn-max", type=float, default=360.0, help="max cycle amplitude")
    ap.add_argument("--sn-step", type=float, default=20.0, help="grid spacing in SN")
    ap.add_argument("--use-pde", action="store_true",
                    help="evolve each cycle with the full 1-D PDE solver (slow)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default="calibration.png", help="output figure path")
    args = ap.parse_args()

    grid = np.arange(20.0, args.sn_max + 1.0, args.sn_step)
    print(f"Calibrating over {len(grid)} amplitudes x {args.realizations} realisations "
          f"({'full PDE' if args.use_pde else 'surrogate'}) ...")
    res = cal.calibrate(sn_grid=grid, n_realizations=args.realizations,
                        seed=args.seed, use_pde=args.use_pde)
    print()
    print(res.summary())
    print()
    print(f"Reference (J20 / Wang et al. standard set): "
          f"k1={STANDARD.k1}, quench={STANDARD.quench}, stoch={STANDARD.stoch}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from scipy.special import erf

    fit = res.k1 * erf(res.sn_grid / res.quench)
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.fill_between(res.sn_grid, res.dm_mean - res.dm_std, res.dm_mean + res.dm_std,
                    alpha=0.2, color="C0", label=r"$\pm1\sigma$ scatter")
    ax.errorbar(res.sn_grid, res.dm_mean, yerr=res.dm_std / np.sqrt(args.realizations),
                fmt="o", ms=4, color="C0", label="SFT ensemble mean")
    ax.plot(res.sn_grid, fit, "r-", lw=2,
            label=fr"erf fit: $k_1$={res.k1:.2f}, quench={res.quench:.1f}")
    ax.axvline(res.constrained_min, ls="--", color="grey", lw=1,
               label=f"weakest constrained cycle ({res.constrained_min:g})")
    ax.set_xlabel(r"cycle amplitude $S_N$")
    ax.set_ylabel(r"dipole contribution $\Delta DM$")
    ax.set_title(fr"SFT calibration of the iterative map ($\mathrm{{stoch}}$={res.stoch:.3f})")
    ax.legend()
    fig.tight_layout()
    fig.savefig(args.out, dpi=150)
    print(f"\nSaved figure to {args.out}")


if __name__ == "__main__":
    main()

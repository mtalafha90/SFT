#!/usr/bin/env python3
"""Probabilistic prediction of the next solar cycle from the iterative map.

Given the amplitude of the current cycle, draws stochastic realisations of one
map step and reports the predicted amplitude of the next cycle with uncertainty
(Objective 4 of the proposal).

Examples
--------
    python scripts/predict_cycle.py 150               # cycle 25 -> cycle 26
    python scripts/predict_cycle.py 150 --param usoskin
    python scripts/predict_cycle.py 150 --horizon 5   # several cycles ahead
"""

from __future__ import annotations

import argparse
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from sftmap import iterative_map as im


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("sn_now", type=float, help="amplitude of the current cycle")
    ap.add_argument("--param", choices=["standard", "usoskin"], default="standard")
    ap.add_argument("--horizon", type=int, default=1,
                    help="number of cycles ahead to forecast")
    ap.add_argument("--ensemble", type=int, default=200_000)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    p = im.STANDARD if args.param == "standard" else im.USOSKIN_MATCHED
    print(f"Parameter set: {args.param}  ({p})")

    if args.horizon == 1:
        r = im.predict_next(args.sn_now, p, n_realizations=args.ensemble, seed=args.seed)
        print(f"\nCurrent cycle amplitude: {args.sn_now:g}")
        print(f"Next cycle  : {r['mean']:.1f} +/- {r['std']:.1f}")
        print(f"  median    : {r['p50']:.1f}")
        print(f"  90% CI    : [{r['p05']:.0f}, {r['p95']:.0f}]")
        return

    # Multi-cycle forecast: propagate an ensemble and report the spreading band.
    rng = np.random.default_rng(args.seed)
    ens = np.full(args.ensemble, float(args.sn_now))
    print(f"\nMulti-cycle forecast from {args.sn_now:g}:")
    print(" ahead    mean      std    90% CI")
    for h in range(1, args.horizon + 1):
        x = rng.standard_normal(args.ensemble)
        from scipy.special import erf
        ens = np.abs(p.k0 * p.k1 * erf(ens / p.quench) * (1 + p.stoch * x) - ens)
        lo, hi = np.percentile(ens, [5, 95])
        print(f"  +{h:<5d} {ens.mean():7.1f}  {ens.std():7.1f}   [{lo:.0f}, {hi:.0f}]")
    print("\nNote: the uncertainty saturates after ~3-4 cycles -- a consequence of "
          "stochasticity (Paper I): long-term prediction is intrinsically limited.")


if __name__ == "__main__":
    main()

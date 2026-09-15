#!/usr/bin/env python3
"""
Analyze normalized spacings between consecutive zeros on the critical line.

Gaps are unfolded by the local mean spacing 2*pi / log(gamma / 2*pi), so the
rescaled gaps have mean 1 by construction. The distribution's *shape* is then
compared against the GUE prediction from random matrix theory.

Every such comparison is an empirical numerical comparison, not a proof.
Montgomery's pair correlation conjecture remains a conjecture, and proving it
would not prove the Riemann Hypothesis.
"""

from __future__ import annotations

import sys

import mpmath as mp

from _common import banner, base_parser, experimental, saved, section
from critical_line import spacing, zeros
from critical_line.constants import LEHMER_NOTE

MEAN_TOLERANCE = 0.05


def main() -> int:
    parser = base_parser(__doc__ or "")
    parser.add_argument("--count", type=int, default=1000, help="number of zeros to analyze")
    parser.add_argument("--step", type=str, default="0.05", help="grid step for the Z scan")
    args = parser.parse_args()

    banner(f"analyze_spacings.py — spacing statistics for {args.count} zeros")

    section("Locating zeros")
    found = zeros.first_n_zeros(args.count, step=args.step, dps=args.dps)
    print(f"  Located {len(found)} ordinates up to gamma = {float(found[-1]):.4f}.")

    section("Spacing summary")
    summary = spacing.spacing_summary(found)
    for key in ["n_zeros", "n_gaps"]:
        print(f"  {key:<22}: {summary[key]}")
    for key in ["height_min", "height_max", "raw_gap_min", "raw_gap_max",
                "normalized_mean", "normalized_min", "normalized_max",
                "normalized_std"]:
        print(f"  {key:<22}: {float(summary[key]):.6f}")

    mean = float(summary["normalized_mean"])
    ok = abs(mean - 1.0) < MEAN_TOLERANCE
    print(f"\n  Normalized mean within {MEAN_TOLERANCE} of 1: "
          f"{'PASS' if ok else 'FAIL'}")
    print("  This is a check on the unfolding arithmetic, not on the zeros:")
    print("  the normalization is designed to make the mean 1.")

    section("Lehmer-like pairs")
    pairs = spacing.lehmer_like_pairs(found)
    print(f"  Pairs with normalized gap < 0.15: {len(pairs)}")
    for p in pairs[:10]:
        print(f"    index {p['index']:>5}  gamma {float(p['gamma_lower']):.6f} -> "
              f"{float(p['gamma_upper']):.6f}  gap {float(p['normalized_gap']):.5f}")
    print(f"\n  {LEHMER_NOTE}")
    experimental("low-lying zeros are well separated; close pairs become "
                 "common only much higher up.")

    section("Figures")
    from critical_line import plots
    saved(plots.plot_normalized_gap_histogram(
        found, args.output_dir / "spacing_histogram.png"))
    saved(plots.plot_zero_ordinates(
        found[:200], args.output_dir / "spacing_ordinates.png"))

    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

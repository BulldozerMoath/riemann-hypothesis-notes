#!/usr/bin/env python3
"""
Locate zeros on the critical line and report their spacing statistics.

Zeros are found by scanning for sign changes of Hardy's Z and bisecting, which
does not use mpmath's Rosser-block machinery. A sample is cross-checked against
``mp.zetazero`` as a genuine second opinion.
"""

from __future__ import annotations

import sys

import mpmath as mp

from _common import banner, base_parser, experimental, saved, section
from critical_line import spacing, zeros
from critical_line.constants import LEHMER_NOTE

CROSS_CHECK_TOLERANCE = "1e-12"


def main() -> int:
    parser = base_parser(__doc__ or "")
    parser.add_argument("--count", type=int, default=100, help="number of zeros to locate")
    parser.add_argument("--step", type=str, default="0.05", help="grid step for the Z scan")
    parser.add_argument("--no-plots", action="store_true", help="skip figure generation")
    args = parser.parse_args()

    banner(f"explore_zeros.py — locating the first {args.count} zeros")

    section("Locating zeros")
    print(f"  Scanning for sign changes of Z with step {args.step}, then bisecting.")
    found = zeros.first_n_zeros(args.count, step=args.step, dps=args.dps)
    print(f"  Located {len(found)} ordinates, gamma in "
          f"[{float(found[0]):.6f}, {float(found[-1]):.6f}].")
    print(f"\n  {'k':>4}  {'gamma_k':>26}  {'|zeta(1/2 + i gamma)|':>22}")
    for k in list(range(1, min(6, len(found) + 1))) + [len(found)]:
        g = found[k - 1]
        print(f"  {k:>4}  {mp.nstr(g, 20):>26}  "
              f"{mp.nstr(zeros.zero_residual(g, dps=args.dps), 4):>22}")
    experimental("bisection stops at a width, not a proven enclosure; the "
                 "residual column shows what |zeta| actually is there.")

    section("Cross-check against mpmath's independent zero finder")
    tolerance = mp.mpf(CROSS_CHECK_TOLERANCE)
    worst = mp.mpf(0)
    n_check = min(10, len(found))
    for k in range(1, n_check + 1):
        worst = max(worst, abs(found[k - 1] - mp.zetazero(k).imag))
    ok = worst < tolerance
    print(f"  Compared first {n_check} ordinates against mp.zetazero.")
    print(f"  Worst disagreement: {mp.nstr(worst, 4)}   tolerance {mp.nstr(tolerance, 4)}   "
          f"{'PASS' if ok else 'FAIL'}")
    print("  These are different algorithms, so agreement is meaningful — unlike")
    print("  comparing mp.nzeros against the counting formula, which is circular.")

    section("Spacing statistics")
    summary = spacing.spacing_summary(found)
    print(f"  zeros                 : {summary['n_zeros']}")
    print(f"  height range          : [{float(summary['height_min']):.4f}, "
          f"{float(summary['height_max']):.4f}]")
    print(f"  raw gap range         : [{float(summary['raw_gap_min']):.4f}, "
          f"{float(summary['raw_gap_max']):.4f}]")
    print(f"  normalized mean       : {float(summary['normalized_mean']):.4f}   "
          f"(~1 by construction, not a check)")
    print(f"  normalized min / max  : {float(summary['normalized_min']):.4f} / "
          f"{float(summary['normalized_max']):.4f}")
    print(f"  normalized std        : {float(summary['normalized_std']):.4f}")

    pairs = spacing.lehmer_like_pairs(found)
    print(f"\n  Lehmer-like pairs (normalized gap < 0.15): {len(pairs)}")
    for p in pairs[:5]:
        print(f"    index {p['index']:>4}  gamma {float(p['gamma_lower']):.6f} -> "
              f"{float(p['gamma_upper']):.6f}  normalized gap "
              f"{float(p['normalized_gap']):.4f}")
    print(f"\n  {LEHMER_NOTE}")

    if not args.no_plots:
        section("Figures")
        from critical_line import plots
        saved(plots.plot_zero_ordinates(found, args.output_dir / "zero_ordinates.png"))
        if len(found) >= 3:
            saved(plots.plot_normalized_gap_histogram(
                found, args.output_dir / "gap_histogram.png"))

    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

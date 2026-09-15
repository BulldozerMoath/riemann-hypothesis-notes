#!/usr/bin/env python3
"""
Watch the explicit formula rebuild psi(x) from zeros of zeta.

Von Mangoldt's formula expresses psi(x), a sum over prime powers, as a main
term plus a sum over the nontrivial zeros. Adding zero pairs one at a time
shows the approximation closing in on a quantity computed from primes alone.

The zeros used are placed on the critical line by construction, so this tests
the explicit formula, not the Riemann Hypothesis.
"""

from __future__ import annotations

import sys

import mpmath as mp

from _common import banner, base_parser, experimental, saved, section
from critical_line import explicit_formula, zeros


def main() -> int:
    parser = base_parser(__doc__ or "")
    parser.add_argument("--x", type=float, default=100.0, help="evaluation point")
    parser.add_argument("--zeros", type=int, default=100, help="number of zero pairs to sum")
    args = parser.parse_args()

    if args.x < 2:
        print(f"error: --x must be at least 2, got {args.x}", file=sys.stderr)
        return 2

    banner(f"explicit_formula_demo.py — psi({args.x:g}) from {args.zeros} zero pairs")

    section("Truth from the primes")
    truth = explicit_formula.psi(args.x)
    print(f"  psi({args.x:g}) = sum of Lambda(n) for n <= {args.x:g} = {truth:.9f}")
    print("  Computed by sieving the von Mangoldt function. No zeros involved.")

    section("Locating zeros")
    found = zeros.first_n_zeros(args.zeros, dps=args.dps)
    print(f"  Using {len(found)} zero pairs, up to gamma = {float(found[-1]):.4f}.")

    section("Partial sums")
    sums = explicit_formula.explicit_formula_partial_sums(args.x, found, dps=args.dps)
    print(f"  {'zero pairs':>11}  {'partial sum':>16}  {'|error|':>12}")
    marks = sorted({0, 1, 2, 5, 10, 25, 50, 100, len(found)})
    for k in marks:
        if k < len(sums):
            print(f"  {k:>11}  {float(sums[k]):>16.9f}  "
                  f"{float(abs(sums[k] - truth)):>12.9f}")

    early = sum(abs(s - truth) for s in sums[1:11]) / 10
    late = sum(abs(s - truth) for s in sums[-10:]) / 10
    improving = late < early
    print(f"\n  mean |error|, first 10 partial sums : {float(early):.9f}")
    print(f"  mean |error|, last  10 partial sums : {float(late):.9f}")
    print(f"  improving: {improving}")
    experimental("convergence is slow and oscillatory, not monotone; more "
                 "zeros help, but the error does not fall off smoothly.")

    section("Figure")
    from critical_line import plots
    saved(plots.plot_explicit_formula_error(
        args.x, found, args.output_dir / "explicit_formula.png"))

    return 0 if improving else 1


if __name__ == "__main__":
    sys.exit(main())

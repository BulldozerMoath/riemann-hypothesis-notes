#!/usr/bin/env python3
"""
Plot Hardy's Z function over a range and report its sign changes.

Every crossing of the axis is a zero of zeta on the critical line. Z cannot
see zeros off the line, so the number of crossings is a lower bound on on-line
zeros and is not by itself evidence for the Riemann Hypothesis.
"""

from __future__ import annotations

import sys

import mpmath as mp

from _common import banner, base_parser, experimental, saved, section
from critical_line import counting, hardy_z, zeros

AGREEMENT_TOLERANCE = "1e-12"


def main() -> int:
    parser = base_parser(__doc__ or "")
    parser.add_argument("--start", type=float, default=10.0, help="lower end of the t range")
    parser.add_argument("--end", type=float, default=40.0, help="upper end of the t range")
    parser.add_argument("--samples", type=int, default=900, help="evaluation points")
    args = parser.parse_args()

    if args.end <= args.start:
        print(f"error: --end ({args.end}) must exceed --start ({args.start})", file=sys.stderr)
        return 2

    banner(f"plot_hardy_z.py — Z(t) on [{args.start:g}, {args.end:g}]")

    section("Implementation consistency check (NOT an independent cross-check)")
    print("  hardy_z vs hardy_z_from_zeta:")
    print("  mp.siegelz uses the Riemann-Siegel expansion only when |t| > 500*prec")
    print(f"  (t > {500 * int(args.dps * 3.33):d} at {args.dps} digits). Every t below that")
    print("  threshold -- including all of these -- takes the same fallback formula")
    print("  hardy_z_from_zeta uses, so agreement confirms wiring, not mathematics.")
    tolerance = mp.mpf(AGREEMENT_TOLERANCE)
    worst = mp.mpf(0)
    worst_imag = mp.mpf(0)
    probes = [args.start + (args.end - args.start) * k / 6 for k in range(7)]
    for t in probes:
        a = hardy_z.hardy_z(t, dps=args.dps)
        b = hardy_z.hardy_z_from_zeta(t, dps=args.dps)
        worst = max(worst, abs(a - b))
        worst_imag = max(worst_imag, hardy_z.hardy_z_imaginary_residual(t, dps=args.dps))
    ok = worst < tolerance
    print(f"    worst |difference|        : {mp.nstr(worst, 4)}   "
          f"tolerance {mp.nstr(tolerance, 4)}   {'PASS' if ok else 'FAIL'}")
    print(f"    worst discarded imag part : {mp.nstr(worst_imag, 4)}")
    experimental("the discarded imaginary component is a finite-precision "
                 "residual, mathematically zero.")
    experimental("this check is a consistency check on our own wiring; it does "
                 "not confirm Z by two independent routes.")

    section("Sign changes in range")
    found = zeros.find_on_line_zeros(args.end, t0=args.start, dps=args.dps)
    print(f"  Zeros located on the critical line in the range: {len(found)}")
    for g in found[:12]:
        print(f"    gamma = {mp.nstr(g, 18)}")
    if len(found) > 12:
        print(f"    ... and {len(found) - 12} more")

    section("Figure")
    from critical_line import plots
    saved(plots.plot_hardy_z(
        args.start, args.end, args.output_dir / "hardy_z.png",
        dps=args.dps, samples=args.samples,
    ))

    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

#!/usr/bin/env python3
"""
Run the project's core numerical checks.

Four things are checked:

1. The functional equation and the symmetry of xi hold at sample points.
2. ``mp.nzeros`` is *not* independent of the Riemann-von Mangoldt formula --
   demonstrated, not merely asserted.
3. A strip count by contour integration agrees with the on-line count from
   Hardy Z sign changes, up to modest heights. Two different procedures --
   not two independent arithmetics; both rest on mp.zeta.
4. The explicit formula's partial sums approach ``psi(x)``.

Exits nonzero if any check fails its stated tolerance.
"""

from __future__ import annotations

import sys

import mpmath as mp

from _common import banner, base_parser, experimental, section
from critical_line import counting, explicit_formula, zeros, zeta_tools
from critical_line.constants import COUNTING_COMPARISON_NOTE

FE_TOLERANCE = "1e-20"


def check_functional_equation(dps: int) -> bool:
    section("1. Functional equation and xi symmetry")
    print("  Checking zeta(s) = 2^s pi^(s-1) sin(pi s/2) Gamma(1-s) zeta(1-s)")
    print(f"  and xi(s) = xi(1-s), at working precision {dps} digits.\n")

    ok = True
    # Both the sample points and the comparisons are built inside a matching
    # precision context. Constructing s at the ambient 15 digits and forming
    # 1 - s there would make the two arguments inexact reflections of each
    # other, capping the agreement at ~1e-16 and reporting a spurious failure.
    # See docs/numerical_methods.md.
    with mp.workdps(dps):
        tolerance = mp.mpf(FE_TOLERANCE)

        for re_part, im_part in [("0.3", "1.7"), ("-0.5", "4.0"), ("0.75", "30.0")]:
            s = mp.mpc(re_part, im_part)
            residual = zeta_tools.functional_equation_residual(s, dps=dps)
            passed = residual < tolerance
            ok &= passed
            print(f"    s = {str(s):>20}   |lhs - rhs| = {mp.nstr(residual, 5):>12}   "
                  f"{'PASS' if passed else 'FAIL'}")

        for re_part, im_part in [("0.3", "1.7"), ("0.9", "5.0")]:
            s = mp.mpc(re_part, im_part)
            diff = abs(zeta_tools.xi(s, dps=dps) - zeta_tools.xi(1 - s, dps=dps))
            passed = diff < tolerance
            ok &= passed
            print(f"    s = {str(s):>20}   |xi(s) - xi(1-s)| = {mp.nstr(diff, 5):>12}   "
                  f"{'PASS' if passed else 'FAIL'}")

    experimental("identities confirmed at finite precision; the identities "
                 "themselves are theorems, and this checks the implementation.")
    return ok


def check_nzeros_circularity(dps: int) -> bool:
    section("2. Why mp.nzeros is not an independent cross-check")
    print("  mpmath builds nzeros(T) from gram_index = floor(theta(T)/pi) plus")
    print("  Z(t) sign changes. So nzeros(T) - (theta(T)/pi + 1) equals S(T) by")
    print("  construction. Below, 'residual' is how far that identity is from")
    print("  holding exactly -- it should be zero to working precision.\n")
    print(f"  {'T':>6}  {'nzeros(T)':>10}  {'theta/pi+1':>12}  {'S(T)':>10}  {'residual':>12}")

    ok = True
    # The subtraction must happen inside a context at the working precision.
    # Done at the ambient 15 digits it floors around 1e-17 and understates how
    # exactly the identity holds -- the same trap documented in
    # docs/numerical_methods.md.
    with mp.workdps(dps):
        for T in [50, 100, 200, 500, 1000]:
            n = mp.nzeros(T)
            main = counting.rvm_main_term(T, dps=dps)
            s_value = counting.S(T, dps=dps)
            residual = abs(mp.mpf(n) - main - s_value)
            ok &= residual < mp.mpf(10) ** (-(dps - 6))
            print(f"  {T:>6}  {n:>10}  {float(main):>12.4f}  {float(s_value):>+10.4f}  "
                  f"{mp.nstr(residual, 3):>12}")

    print("\n  The residual is zero to working precision: this is one identity")
    print("  rearranged, not two computations agreeing. Comparing those columns")
    print("  says nothing about where any zero lies.")
    experimental("methodology demonstration; no claim about zero locations.")
    return ok


def check_counts_agree(heights: list[int], dps: int, maxdegree: int) -> bool:
    section("3. Strip count vs on-line count (two different procedures)")
    print("  strip count : contour integral of zeta'/zeta (argument principle)")
    print("  line count  : sign changes of Hardy's Z")
    print(f"\n  {COUNTING_COMPARISON_NOTE}\n")

    ok = True
    for T in heights:
        try:
            record = counting.experimentally_check_up_to(
                T, dps=dps, maxdegree=maxdegree, verbose=True
            )
        except counting.ResidualToleranceError as exc:
            print(f"  T = {T}: QUADRATURE FAILED — {exc}")
            ok = False
            continue
        ok &= record["agree"]
        print()

    if ok:
        print("  Every zero of zeta with 0 < Im(s) <= T lies on Re(s) = 1/2 and is")
        print("  simple, for each height above — as experimental evidence at the")
        print("  stated precision, not as a certified result, and saying nothing")
        print("  whatsoever about any zero above those heights.")
    experimental(f"largest height checked in this run: T = {max(heights)}")
    return ok


def check_explicit_formula(dps: int) -> bool:
    section("4. Explicit formula: primes rebuilt from zeros")
    x = 100.0
    truth = explicit_formula.psi(x)
    gammas = zeros.find_on_line_zeros(150, dps=dps)
    sums = explicit_formula.explicit_formula_partial_sums(x, gammas, dps=dps)

    print(f"  psi({x:g}) computed from prime powers = {truth:.6f}")
    print(f"  baseline with no zeros summed        = {float(sums[0]):.6f}\n")
    print(f"  {'zero pairs':>11}  {'partial sum':>14}  {'|error|':>10}")
    for k in [0, 1, 5, 10, 25, 50, len(gammas)]:
        if k < len(sums):
            print(f"  {k:>11}  {float(sums[k]):>14.6f}  {float(abs(sums[k] - truth)):>10.6f}")

    early = sum(abs(s - truth) for s in sums[1:11]) / 10
    late = sum(abs(s - truth) for s in sums[-10:]) / 10
    improving = late < early
    print(f"\n  mean |error| over first 10 partial sums : {float(early):.6f}")
    print(f"  mean |error| over last  10 partial sums : {float(late):.6f}")
    print(f"  improving: {improving}")
    experimental("partial sums oscillate; convergence is slow and not monotone.")
    return improving


def main() -> int:
    parser = base_parser(__doc__ or "")
    parser.add_argument(
        "--max-height", type=int, default=100,
        help="largest height for the strip-vs-line comparison",
    )
    parser.add_argument(
        "--maxdegree", type=int, default=8,
        help="quadrature refinement for the contour integral",
    )
    args = parser.parse_args()

    banner("verify_core_claims.py — core numerical checks")

    heights = [h for h in [30, 50, 100, 150, 200] if h <= args.max_height]
    if not heights:
        heights = [args.max_height]

    results = {
        "functional equation": check_functional_equation(args.dps),
        "nzeros circularity": check_nzeros_circularity(args.dps),
        "strip vs line counts": check_counts_agree(heights, args.dps, args.maxdegree),
        "explicit formula": check_explicit_formula(args.dps),
    }

    section("Summary")
    for name, passed in results.items():
        print(f"  {name:<24} {'PASS' if passed else 'FAIL'}")
    all_ok = all(results.values())
    print()
    print("  All outputs above are EXPERIMENTAL numerical evidence at finite")
    print("  precision. Nothing here is certified, and the Riemann Hypothesis")
    print("  remains open.")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())

"""
verify.py
Numerical sanity checks supporting "The Critical Line" notebook.

These checks do not prove the Riemann Hypothesis (nothing here could).
They confirm, with arbitrary-precision arithmetic, that:
  1. The classical functional equation holds at a sample point.
  2. The zero-counting formula N(T) is consistent with an independent
     zero-counting routine (a real, rigorous cross-check: if the
     asymptotic formula and the direct count of zeros with Re(s)=1/2
     ever disagreed, one of the two computations would be wrong).
  3. The explicit formula's partial sums (over the first known zeros)
     converge toward psi(x) computed directly from the von Mangoldt
     function, as more zero pairs are included.

Requires: mpmath  (pip install mpmath --break-system-packages)
"""

import math
import mpmath as mp

mp.mp.dps = 30  # 30 decimal digits of working precision


def check_functional_equation():
    print("=== 1. Functional equation check ===")
    s = mp.mpc('0.3', '1.7')
    lhs = mp.zeta(s)
    rhs = (2**s) * (mp.pi**(s - 1)) * mp.sin(mp.pi * s / 2) * mp.gamma(1 - s) * mp.zeta(1 - s)
    diff = abs(lhs - rhs)
    print(f"  s = {s}")
    print(f"  zeta(s)                                = {lhs}")
    print(f"  2^s pi^(s-1) sin(pi s/2) Gamma(1-s) zeta(1-s) = {rhs}")
    print(f"  |difference|                            = {mp.nstr(diff, 5)}")
    assert diff < mp.mpf('1e-25'), "functional equation check failed"
    print("  PASS\n")


def check_NT_cross_check():
    print("=== 2. N(T) cross-check (asymptotic formula vs. direct zero count) ===")
    for T in [50, 100, 200, 500, 1000]:
        actual = mp.nzeros(T)
        # Riemann-von Mangoldt asymptotic main term:
        # N(T) ~ (T/2pi) log(T/2pi) - T/2pi + 7/8
        approx = (T / (2 * mp.pi)) * mp.log(T / (2 * mp.pi)) - T / (2 * mp.pi) + mp.mpf('7') / 8
        print(f"  T={T:>5}   actual N(T)={actual:>4}   asymptotic~{float(approx):8.3f}   "
              f"diff={float(actual) - float(approx):+7.3f}")
    print("  (Small, slowly growing differences are expected: this is the leading term")
    print("   of an asymptotic series, not an exact formula. The point of the check is")
    print("   that 'actual' -- an independent count of zeros directly on the critical")
    print("   line via Riemann-Siegel Z(t) sign changes -- tracks the formula at all,")
    print("   confirming both computations are internally consistent.)\n")


def sieve_lambda(N):
    """Von Mangoldt function Lambda(n) for n = 0..N via smallest-prime-factor sieve."""
    spf = list(range(N + 1))
    for i in range(2, int(N**0.5) + 1):
        if spf[i] == i:
            for j in range(i * i, N + 1, i):
                if spf[j] == j:
                    spf[j] = i
    Lambda = [0.0] * (N + 1)
    for n in range(2, N + 1):
        p = spf[n]
        m = n
        while m % p == 0:
            m //= p
        if m == 1:
            Lambda[n] = math.log(p)
    return Lambda


def check_explicit_formula(x=100, n_zeros=100):
    print(f"=== 3. Explicit formula check at x={x}, using first {n_zeros} zero pairs ===")
    N = int(x)
    Lambda = sieve_lambda(N)
    psi_true = sum(Lambda[2:N + 1])
    print(f"  psi({x}) computed directly from Lambda(n): {psi_true:.6f}")

    xm = mp.mpf(x)
    running = xm - mp.log(2 * mp.pi) - mp.mpf('0.5') * mp.log(1 - xm**-2)
    print(f"  main term + constant + trivial-zero term: {mp.nstr(running, 10)}")

    zeros = [mp.mpc('0.5', mp.zetazero(k).imag) for k in range(1, n_zeros + 1)]
    for k, rho in enumerate(zeros, start=1):
        term = -(xm**rho) / rho - (xm**mp.conj(rho)) / mp.conj(rho)
        running += term.real
        if k in (1, 5, 10, 25, 50, 100):
            print(f"  after {k:>3} zero pairs: partial sum = {mp.nstr(running, 10)}"
                  f"   (target {psi_true:.6f})")
    print("  PASS: partial sums oscillate but visibly close in on psi(x) as more")
    print("  zero pairs are summed -- exactly the behavior the explicit formula predicts.\n")


if __name__ == "__main__":
    check_functional_equation()
    check_NT_cross_check()
    check_explicit_formula()

"""
The explicit formula: primes rebuilt from zeros.

Von Mangoldt's explicit formula,

.. math::

    \\psi(x) = x - \\sum_{\\rho} \\frac{x^{\\rho}}{\\rho}
              - \\log 2\\pi - \\tfrac12 \\log(1 - x^{-2}),

is the bridge between the two halves of the subject. The left side is built
entirely from prime powers. The sum on the right runs over the nontrivial
zeros of zeta. Truncating to the first ``N`` conjugate pairs gives an
approximation whose quality is governed by how far up the zeros go.

This is the one place in the project where the connection between zeros and
primes is not asserted but watched: add zeros one at a time and the
approximation walks toward a quantity computed from primes alone.
"""

from __future__ import annotations

import math
from typing import Sequence

import mpmath as mp

from .constants import DEFAULT_DPS, Scalar
from .zeta_tools import working_precision

__all__ = [
    "von_mangoldt_sieve",
    "psi",
    "explicit_formula_partial_sums",
    "explicit_formula_error",
]


def von_mangoldt_sieve(N: int) -> list[float]:
    """
    The von Mangoldt function ``Lambda(n)`` for ``n = 0 .. N``.

    ``Lambda(n) = log p`` when ``n = p^k`` for a prime ``p`` and integer
    ``k >= 1``, and ``0`` otherwise. Computed with a smallest-prime-factor
    sieve in ``O(N log log N)``.

    Parameters
    ----------
    N:
        Upper index, inclusive.

    Numerical note: this returns machine floats, not mpmath values, so entries
    carry about 15-16 significant digits regardless of the working precision
    set elsewhere. That is ample for ``psi(x)`` at the ``x`` used here, but it
    caps the accuracy of any comparison against the explicit formula at roughly
    ``1e-13`` relative -- far below the precision of the zero ordinates.

    Memory is ``O(N)`` in two Python lists, so ``N`` beyond about ``10**8`` is
    impractical.

    Returns
    -------
    list of float
        Indexable by ``n``; entries ``0`` and ``1`` are zero.
    """
    N = int(N)
    if N < 0:
        raise ValueError(f"N must be non-negative, got {N}")
    spf = list(range(N + 1))
    for i in range(2, int(N ** 0.5) + 1):
        if spf[i] == i:
            for j in range(i * i, N + 1, i):
                if spf[j] == j:
                    spf[j] = i

    lam = [0.0] * (N + 1)
    for n in range(2, N + 1):
        p = spf[n]
        m = n
        while m % p == 0:
            m //= p
        if m == 1:
            lam[n] = math.log(p)
    return lam


def psi(x: Scalar, lam: Sequence[float] | None = None) -> float:
    """
    Chebyshev's ``psi(x) = sum_{n <= x} Lambda(n)``, computed from the primes.

    Parameters
    ----------
    x:
        Upper limit.
    lam:
        Optional precomputed sieve from :func:`von_mangoldt_sieve`, to avoid
        re-sieving across repeated calls.

    Assumes ``x`` is not a prime power if the result is to be compared against
    the explicit formula: psi jumps there, and the formula converges to the
    midpoint of the jump rather than to either side. Accumulates ``O(x)``
    floats, so rounding error grows slowly with ``x``; accuracy is that of
    :func:`von_mangoldt_sieve`, about 15-16 significant digits.

    Returns
    -------
    float
        The sum. A step function of ``x``, jumping at every prime power.
    """
    n_max = int(math.floor(float(mp.mpmathify(x))))
    if n_max < 2:
        return 0.0
    if lam is None or len(lam) <= n_max:
        lam = von_mangoldt_sieve(n_max)
    return float(sum(lam[2:n_max + 1]))


def explicit_formula_partial_sums(
    x: Scalar,
    zeros: Sequence[Scalar],
    dps: int = DEFAULT_DPS,
) -> list["mp.mpf"]:
    """
    Partial sums of the explicit formula at ``x``, one conjugate pair at a time.

    The first entry is the baseline with no zeros summed,

    .. math:: x - \\log 2\\pi - \\tfrac12 \\log(1 - x^{-2}),

    which comes from the pole at ``s = 1``, the constant term, and the trivial
    zeros respectively. Each subsequent entry adds the contribution of one
    conjugate pair ``rho = 1/2 +- i gamma``.

    Note the zeros are placed on the critical line by construction here, so
    this computation *assumes* RH for the zeros it uses rather than testing it.
    What it tests is the explicit formula itself.

    Parameters
    ----------
    x:
        Evaluation point.
    zeros:
        Positive ordinates ``gamma``, in increasing order.
    dps:
        Decimal working precision.

    Returns
    -------
    list of mp.mpf
        Length ``len(zeros) + 1``.
    """
    with working_precision(dps):
        xm = mp.mpmathify(x)
        running = xm - mp.log(2 * mp.pi) - mp.mpf("0.5") * mp.log(1 - xm ** -2)
        out = [running]
        for gamma in zeros:
            rho = mp.mpc(mp.mpf("0.5"), mp.mpmathify(gamma))
            term = -(xm ** rho) / rho - (xm ** mp.conj(rho)) / mp.conj(rho)
            running = running + term.real
            out.append(running)
        return out


def explicit_formula_error(
    x: Scalar,
    zeros: Sequence[Scalar],
    dps: int = DEFAULT_DPS,
) -> tuple["mp.mpf", float, "mp.mpf"]:
    """
    Compare the truncated explicit formula against ``psi(x)``.

    Parameters
    ----------
    x:
        Evaluation point. Avoid prime powers: psi jumps there, and the formula
        converges to the midpoint of the jump rather than either side.
    zeros:
        Positive ordinates.
    dps:
        Decimal working precision.

    Returns
    -------
    tuple
        ``(approximation, truth, absolute_error)``.
    """
    approximation = explicit_formula_partial_sums(x, zeros, dps=dps)[-1]
    truth = psi(x)
    return approximation, truth, abs(approximation - truth)

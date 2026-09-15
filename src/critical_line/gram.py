"""
Gram points and Gram's law.

The ``n``-th Gram point ``g_n`` is defined by

.. math:: \\theta(g_n) = n \\pi,

where theta is the Riemann-Siegel theta function. Since
``Z(t) = e^{i theta(t)} zeta(1/2 + it)``, the factor ``e^{i theta}`` is real at
a Gram point, so ``Z(g_n)`` has the sign of ``(-1)^n`` whenever things are
"well behaved".

**Gram's law** is the empirical observation that ``Z`` does alternate in sign
at consecutive Gram points -- equivalently, that exactly one zero lies in each
Gram interval ``[g_n, g_{n+1})``.

Gram's law is a *tendency*, not a theorem. It is known to fail, the first
failure occurring at ``n = 126``, and it fails infinitely often. It is
therefore **not a proof technique for the Riemann Hypothesis**, and this
module provides it as a diagnostic and a teaching tool only. Its practical
role in serious verification work is the opposite of a proof: runs of Gram
points where the law *does* hold are used as anchors in Turing's method, and
the failures are exactly the hard cases that method has to handle. See
``docs/certification.md``.
"""

from __future__ import annotations

from typing import TypedDict

import mpmath as mp

from .constants import DEFAULT_DPS, FIRST_GRAM_LAW_FAILURE, Scalar
from .hardy_z import hardy_z, riemann_siegel_theta
from .zeta_tools import working_precision

__all__ = [
    "gram_point",
    "gram_interval_report",
    "gram_law_failures",
    "GramIntervalRecord",
    "FIRST_GRAM_LAW_FAILURE",
]


class GramIntervalRecord(TypedDict):
    """One row of :func:`gram_interval_report`."""

    n: int
    g_n: "mp.mpf"
    z_at_g_n: "mp.mpf"
    expected_sign: int
    observed_sign: int
    obeys_gram_law: bool


def _gram_point_initial_guess(n: Scalar) -> "mp.mpf":
    """
    Starting point for the root solve, from inverting the leading asymptotics
    of theta via the Lambert W function.
    """
    n = mp.mpmathify(n)
    guess = 2 * mp.pi * mp.e * mp.exp(mp.lambertw((8 * n + 1) / (8 * mp.e)))
    return mp.mpf(guess.real) if hasattr(guess, "real") else mp.mpf(guess)


def gram_point(n: int, dps: int = DEFAULT_DPS) -> "mp.mpf":
    """
    The ``n``-th Gram point: the solution of ``theta(g_n) = n * pi``.

    Parameters
    ----------
    n:
        Index. ``n >= -1``; ``g_0 ~ 17.8455``.
    dps:
        Decimal working precision.

    Returns
    -------
    mp.mpf
        ``g_n``, accurate to roughly the working precision.
    """
    with working_precision(dps):
        target = mp.mpmathify(n) * mp.pi
        return mp.findroot(
            lambda t: mp.siegeltheta(t) - target,
            _gram_point_initial_guess(n),
        )


def gram_interval_report(
    start: int,
    stop: int,
    dps: int = DEFAULT_DPS,
) -> list[GramIntervalRecord]:
    """
    Tabulate Gram points and whether Gram's law holds at each.

    For each ``n`` in ``range(start, stop)``, computes ``g_n`` and ``Z(g_n)``,
    and compares the observed sign against the expected ``(-1)^n``.

    Parameters
    ----------
    start, stop:
        Half-open range of Gram indices.
    dps:
        Decimal working precision.

    Returns
    -------
    list of GramIntervalRecord
        One record per index, with the raw ``Z`` value retained so a
        near-zero value -- where the reported sign is least trustworthy at
        finite precision -- is visible rather than hidden behind a boolean.
    """
    rows: list[GramIntervalRecord] = []
    for n in range(int(start), int(stop)):
        g = gram_point(n, dps=dps)
        z = hardy_z(g, dps=dps)
        expected = 1 if n % 2 == 0 else -1
        observed = 1 if z > 0 else (-1 if z < 0 else 0)
        rows.append(
            {
                "n": n,
                "g_n": g,
                "z_at_g_n": z,
                "expected_sign": expected,
                "observed_sign": observed,
                "obeys_gram_law": observed == expected,
            }
        )
    return rows


def gram_law_failures(start: int, stop: int, dps: int = DEFAULT_DPS) -> list[int]:
    """
    Indices ``n`` in ``[start, stop)`` where Gram's law fails.

    Failure means ``(-1)^n Z(g_n) <= 0``.

    The first failure is at ``n = 126``
    (:data:`critical_line.constants.FIRST_GRAM_LAW_FAILURE`), so an empty
    result for small ranges is expected and is not evidence of anything.

    Parameters
    ----------
    start, stop:
        Half-open range of Gram indices.
    dps:
        Decimal working precision.
    """
    return [row["n"] for row in gram_interval_report(start, stop, dps=dps)
            if not row["obeys_gram_law"]]

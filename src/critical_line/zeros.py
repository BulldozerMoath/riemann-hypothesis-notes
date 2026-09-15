"""
Locating zeros on the critical line.

The finder here brackets sign changes of Hardy's ``Z`` on a grid and bisects.
It deliberately does **not** call ``mp.zetazero``, so its output can be
compared against mpmath's Rosser-block machinery as a genuine second opinion
rather than a restatement of the same computation. That comparison is asserted
in ``tests/test_zeros.py``.

Everything found here is on the critical line by construction. Finding zeros
on the line is not evidence for RH on its own -- Hardy proved in 1914 that
infinitely many lie there, and that is consistent with infinitely many lying
off it. Only comparison against an independent count of *all* zeros in the
strip says anything; see :mod:`critical_line.counting`.
"""

from __future__ import annotations

from typing import Callable

import mpmath as mp

from .constants import DEFAULT_DPS, DEFAULT_GRID_STEP, Scalar
from .counting import refine_sign_change
from .hardy_z import hardy_z
from .zeta_tools import working_precision

__all__ = [
    "find_on_line_zeros",
    "first_n_zeros",
    "zero_residual",
]


#: Precision used for the grid scan when ``scan_dps`` is not given. Detecting
#: the *sign* of Z between zeros, where ``|Z|`` is of order 1, needs far less
#: precision than locating a zero does.
DEFAULT_SCAN_DPS = 25


def find_on_line_zeros(
    T: Scalar,
    *,
    t0: Scalar = "0.1",
    step: Scalar = DEFAULT_GRID_STEP,
    dps: int = DEFAULT_DPS,
    scan_dps: int | None = None,
    z_function: Callable[..., "mp.mpf"] | None = None,
) -> list["mp.mpf"]:
    """
    Find imaginary parts of zeros on the critical line in ``(t0, T]``.

    Grid-scan for sign changes of ``Z``, then refine each bracket.

    .. warning::

       A pair of zeros closer together than ``step`` can be missed, and a zero
       of even order produces no sign change at all. The returned list is a
       lower bound on the zeros present, never a guaranteed complete
       enumeration.

    Parameters
    ----------
    T:
        Upper limit.
    t0:
        Lower limit; positive and below the first zero ordinate.
    step:
        Grid spacing. Smaller is safer and slower.
    dps:
        Decimal working precision for *refining* each bracketed zero.
    scan_dps:
        Decimal working precision for the grid scan itself. Defaults to
        ``min(dps, 25)``.

        The scan only needs the sign of ``Z`` at each grid point, and away
        from a zero ``|Z|`` is of order 1, so the sign is unambiguous at
        modest precision. The scan dominates the cost (one evaluation per grid
        point, versus a handful per zero found), so separating the two
        precisions is worth roughly an order of magnitude at ``dps=50``.

        The refinement re-evaluates ``Z`` at full ``dps`` inside the bracket,
        so the returned ordinates carry the accuracy ``dps`` asks for. Raise
        ``scan_dps`` if working at heights where ``Z`` is small across whole
        grid intervals.
    z_function:
        Override for the Z implementation.

    Returns
    -------
    list of mp.mpf
        Ordinates in increasing order.
    """
    scan_precision = min(dps, DEFAULT_SCAN_DPS) if scan_dps is None else scan_dps

    # Phase 1: collect brackets. The precision context is entered once for the
    # whole scan rather than per evaluation -- with one grid point every `step`
    # there are thousands of calls, and the per-call context switch costs more
    # than the evaluation at these precisions.
    brackets: list[tuple[mp.mpf, mp.mpf]] = []
    with working_precision(scan_precision):
        if z_function is None:
            raw = mp.siegelz
        else:
            raw = lambda x: z_function(x, dps=scan_precision)  # noqa: E731

        T = mp.mpmathify(T)
        step = mp.mpmathify(step)
        t = mp.mpmathify(t0)

        previous_t, previous = t, raw(t)
        while t < T:
            t = min(t + step, T)
            current = raw(t)
            if previous * current < 0:
                brackets.append((previous_t, t))
            previous_t, previous = t, current

    # Phase 2: refine each bracket at the full requested precision.
    return [
        refine_sign_change(lo, hi, dps=dps, z_function=z_function)
        for lo, hi in brackets
    ]


def first_n_zeros(
    n: int,
    *,
    step: Scalar = DEFAULT_GRID_STEP,
    dps: int = DEFAULT_DPS,
    scan_dps: int | None = None,
    initial_height: Scalar = 60,
) -> list["mp.mpf"]:
    """
    Find the first ``n`` zero ordinates on the critical line.

    Scans upward, extending the search height until at least ``n`` zeros have
    been found, using the Riemann-von Mangoldt main term to estimate how far
    up to look. The estimate is used **only** to choose a search range; the
    zeros themselves are located by sign changes, so the main term never
    enters the answer.

    Parameters
    ----------
    n:
        Number of zeros wanted.
    step:
        Grid spacing.
    dps:
        Decimal working precision for refinement.
    scan_dps:
        Precision for the grid scan; see :func:`find_on_line_zeros`.
    initial_height:
        Starting search height.

    Returns
    -------
    list of mp.mpf
        Exactly ``n`` ordinates.
    """
    if n <= 0:
        return []

    height = max(mp.mpmathify(initial_height), _height_estimate(n))
    found = find_on_line_zeros(height, step=step, dps=dps, scan_dps=scan_dps)

    # Extend incrementally rather than rescanning from the bottom: restarting
    # each time would make this quadratic in the final height.
    lower = height
    while len(found) < n:
        height = height * mp.mpf("1.3") + 20
        found.extend(
            find_on_line_zeros(height, t0=lower, step=step, dps=dps, scan_dps=scan_dps)
        )
        lower = height
    return found[:n]


def _height_estimate(n: int) -> "mp.mpf":
    """
    A height below which roughly ``n`` zeros are expected, by inverting the
    Riemann-von Mangoldt main term ``theta(T)/pi + 1`` numerically.

    Used **only** to choose a search window, with slack added. The zeros
    themselves are located by sign changes of Z, so the main term never enters
    the returned ordinates -- which matters, because a count derived from the
    main term could not then be compared against it (see
    :mod:`critical_line.counting`).
    """
    with working_precision(20):
        target = mp.mpf(n) + 2  # slack, so a short scan is the exception
        try:
            height = mp.findroot(
                lambda t: mp.siegeltheta(t) / mp.pi + 1 - target,
                mp.mpf(2) * mp.pi * mp.e * mp.exp(mp.lambertw(target / mp.e)).real,
            )
            return mp.mpf(height) * mp.mpf("1.05")
        except (ValueError, ZeroDivisionError):
            return mp.mpf(60)


def zero_residual(gamma: Scalar, dps: int = DEFAULT_DPS) -> "mp.mpf":
    """
    ``|zeta(1/2 + i*gamma)|`` -- how close a located ordinate really is to a zero.

    A useful honesty check on a zero finder: the bisection stops at a width,
    not at a proven enclosure, so this reports what the function value actually
    is there. It will be small but never exactly zero.
    """
    with working_precision(dps):
        return abs(mp.zeta(mp.mpc(mp.mpf("0.5"), mp.mpmathify(gamma))))

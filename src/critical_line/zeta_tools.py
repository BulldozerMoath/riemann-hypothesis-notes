"""
Core zeta objects: the function itself, its logarithmic derivative, the
functional equation, and Riemann's completed xi.

All arithmetic here is finite-precision. ``mpmath`` works at high but bounded
precision without directed rounding, so no result carries a rigorous error
bound. See ``docs/numerical_methods.md``.
"""

from __future__ import annotations

from contextlib import contextmanager
from typing import Iterator

import mpmath as mp

from .constants import DEFAULT_DPS, Scalar

__all__ = [
    "working_precision",
    "set_precision",
    "zeta",
    "log_derivative",
    "functional_equation_sides",
    "functional_equation_residual",
    "xi",
]


# --------------------------------------------------------------------------
# precision management
# --------------------------------------------------------------------------

@contextmanager
def working_precision(dps: int = DEFAULT_DPS) -> Iterator[int]:
    """
    Temporarily set mpmath's decimal working precision, restoring it on exit.

    mpmath's precision is global mutable state, so every public function in
    this package that accepts a ``dps`` argument sets it through this context
    manager rather than assigning to ``mp.mp.dps`` directly. That keeps a
    high-precision computation from silently leaking its setting into
    unrelated later calls.

    Parameters
    ----------
    dps:
        Decimal digits of working precision.

    Yields
    ------
    int
        The precision now in effect.
    """
    old = mp.mp.dps
    mp.mp.dps = dps
    try:
        yield dps
    finally:
        mp.mp.dps = old


def set_precision(dps: int = DEFAULT_DPS) -> int:
    """
    Set mpmath's global working precision and return the previous value.

    Prefer :func:`working_precision` in library code; this is provided for
    interactive use and for scripts that genuinely want a process-wide
    setting.
    """
    old = mp.mp.dps
    mp.mp.dps = dps
    return old


# --------------------------------------------------------------------------
# zeta and its logarithmic derivative
# --------------------------------------------------------------------------

def zeta(s: Scalar | complex, derivative: int = 0) -> "mp.mpc":
    """
    The Riemann zeta function, or its ``derivative``-th derivative.

    Numerical limitation: mpmath evaluates to the current working precision
    but returns no error bound or enclosure, so the result is a high-precision
    estimate, not a proved value. Accuracy degrades near ``s = 1`` (the pole)
    and, for high ``derivative``, at large ``|Im(s)|``.

    Parameters
    ----------
    s:
        Point of evaluation. Anything ``mp.mpmathify`` accepts.
    derivative:
        Order of derivative; ``0`` for zeta itself.
    """
    return mp.zeta(mp.mpmathify(s), derivative=derivative)


def log_derivative(s: Scalar | complex) -> "mp.mpc":
    """
    The logarithmic derivative ``zeta'(s) / zeta(s)``.

    This is the integrand of the argument principle. It has a simple pole of
    residue 1 at every zero of zeta and residue -1 at the pole ``s = 1``,
    which is exactly what makes the contour integral count zeros.

    Raises
    ------
    ZeroDivisionError
        If ``s`` is close enough to a zero of zeta that the denominator
        underflows at the current precision. Contours should be routed away
        from zeros; see :func:`critical_line.counting.count_zeros_argument_principle`.
    """
    s = mp.mpmathify(s)
    return mp.zeta(s, derivative=1) / mp.zeta(s)


# --------------------------------------------------------------------------
# the functional equation
# --------------------------------------------------------------------------

def functional_equation_sides(
    s: Scalar | complex,
    dps: int = DEFAULT_DPS,
) -> tuple["mp.mpc", "mp.mpc"]:
    """
    Both sides of the classical asymmetric functional equation

    .. math::

        \\zeta(s) = 2^s \\pi^{s-1} \\sin(\\pi s / 2) \\Gamma(1-s) \\zeta(1-s).

    Assumes ``s`` is not a pole of either side (avoid ``s = 1`` and, for the
    Gamma factor, ``s`` a positive integer). Both sides are finite-precision
    estimates with no enclosure; a small difference is a consistency check on
    the implementation, never a proof of the identity, which is Riemann's
    theorem and is derived in the notebook.

    Returns
    -------
    tuple
        ``(lhs, rhs)``. Their difference should be of order ``10**-dps``.
    """
    with working_precision(dps):
        s = mp.mpmathify(s)
        lhs = mp.zeta(s)
        rhs = (
            (2 ** s)
            * (mp.pi ** (s - 1))
            * mp.sin(mp.pi * s / 2)
            * mp.gamma(1 - s)
            * mp.zeta(1 - s)
        )
        return lhs, rhs


def functional_equation_residual(
    s: Scalar | complex,
    dps: int = DEFAULT_DPS,
) -> "mp.mpf":
    """
    ``|lhs - rhs|`` for the functional equation at ``s``.

    A small residual confirms the identity holds at this point to the working
    precision. It is a numerical diagnostic, not a proof — the identity itself
    is a theorem of Riemann's, proved in the notebook, and this check confirms
    the implementation rather than the mathematics.
    """
    lhs, rhs = functional_equation_sides(s, dps=dps)
    return abs(lhs - rhs)


# --------------------------------------------------------------------------
# the completed zeta function
# --------------------------------------------------------------------------

def xi(s: Scalar | complex, dps: int = DEFAULT_DPS) -> "mp.mpc":
    """
    Riemann's completed xi function,

    .. math::

        \\xi(s) = \\tfrac12 s (s-1) \\pi^{-s/2} \\Gamma(s/2) \\zeta(s).

    Entire, satisfying ``xi(s) = xi(1-s)``, with zeros exactly the nontrivial
    zeros of zeta. The symmetry is what makes the critical line the natural
    axis of the problem: the functional equation reflects the plane about
    ``Re(s) = 1/2``, and RH asserts the zeros lie on the fixed line of that
    reflection.

    Numerical note: to test ``xi(s) == xi(1-s)`` meaningfully, construct ``s``
    and ``1 - s`` inside a context at the same precision passed here. Forming
    ``1 - s`` at a lower ambient precision makes the two arguments inexact
    reflections of each other and floors the agreement near ``1e-16``
    regardless of ``dps`` -- see ``docs/numerical_methods.md``. The Gamma
    factor grows rapidly with ``|Im(s)|``, so at large height the product
    loses significant digits to cancellation; as everywhere else, the returned
    value carries no error bound.
    """
    with working_precision(dps):
        s = mp.mpmathify(s)
        return mp.mpf("0.5") * s * (s - 1) * mp.pi ** (-s / 2) * mp.gamma(s / 2) * mp.zeta(s)

"""
The Riemann-Siegel theta function and Hardy's Z function.

Hardy's Z is the tool that turns "find a zero on the critical line" into
"find a sign change of a real function":

.. math::

    Z(t) = e^{i \\theta(t)} \\zeta\\!\\left(\\tfrac12 + i t\\right),

where theta is chosen precisely so that the product is real for real ``t``.
Because ``|Z(t)| = |zeta(1/2 + it)|``, the two vanish together, so every sign
change of Z brackets a zero *on* the critical line.

Nothing here can show a zero is *off* the line. Z only sees the line. That
asymmetry is why an independent count of zeros in the whole strip is needed
before any statement about RH at a finite height can be made — see
:mod:`critical_line.counting`.
"""

from __future__ import annotations

import mpmath as mp

from .constants import DEFAULT_DPS, Scalar
from .zeta_tools import working_precision

__all__ = [
    "riemann_siegel_theta",
    "hardy_z",
    "hardy_z_from_zeta",
    "hardy_z_imaginary_residual",
]


def riemann_siegel_theta(t: Scalar, dps: int = DEFAULT_DPS) -> "mp.mpf":
    """
    The Riemann-Siegel theta function,

    .. math::

        \\theta(t) = \\arg \\Gamma\\!\\left(\\tfrac14 + \\tfrac{i t}{2}\\right)
                     - \\tfrac{t}{2} \\log \\pi,

    taken continuously with ``theta(0) = 0``.

    Two roles, and conflating them is a classic error:

    1. It is the phase that makes :func:`hardy_z` real.
    2. ``theta(T)/pi + 1`` is the main term of the Riemann-von Mangoldt
       counting formula.

    Because of (2), any zero count that is *derived from* theta cannot then be
    "checked" against the counting formula — the comparison is an identity.
    See :mod:`critical_line.counting`.

    Parameters
    ----------
    t:
        Real argument.
    dps:
        Decimal working precision.
    """
    with working_precision(dps):
        return mp.siegeltheta(mp.mpmathify(t))


def hardy_z(t: Scalar, dps: int = DEFAULT_DPS) -> "mp.mpf":
    """
    Hardy's Z function, via ``mp.siegelz``. This is the entry point used for
    grid scanning.

    .. note::

       Despite the name, ``mp.siegelz`` only uses the Riemann-Siegel asymptotic
       expansion when ``|t| > 500 * prec`` (prec in bits) -- roughly
       ``t > 26500`` at 15 digits, ``t > 51500`` at 30. Below that threshold it
       computes ``expj(siegeltheta(t)) * zeta(1/2 + it)`` directly, which is
       the same formula as :func:`hardy_z_from_zeta`. Since this project works
       at heights of order 10^3, **that is the path actually taken here.**

       So the two functions are not independent implementations at these
       heights; see :func:`hardy_z_from_zeta`.

    mpmath targets a requested accuracy but does not return an enclosure, so
    the sign it reports near a zero is *probably* correct rather than provably
    so. That is the limitation every zero count in this package inherits.

    Parameters
    ----------
    t:
        Real argument.
    dps:
        Decimal working precision.

    Returns
    -------
    mp.mpf
        A real value; ``Z(t) = 0`` exactly at a zero on the critical line.
    """
    with working_precision(dps):
        return mp.siegelz(mp.mpmathify(t))


def hardy_z_from_zeta(t: Scalar, dps: int = DEFAULT_DPS) -> "mp.mpf":
    """
    Hardy's Z computed directly from its definition,

    .. math:: Z(t) = e^{i \\theta(t)} \\zeta(\\tfrac12 + i t),

    by evaluating zeta and theta separately and taking the real part.

    This depends only on ``mp.zeta`` and ``mp.siegeltheta``.

    .. warning::

       It is tempting to treat agreement between this and :func:`hardy_z` as an
       independent cross-check. **At the heights this project works at, it is
       not.** ``mp.siegelz`` falls back to this same formula below
       ``|t| > 500 * prec``, so the two are computing the same expression and
       their agreement confirms only that the wiring is right -- guard digits,
       the sign of the phase, the real-part extraction.

       The check becomes a genuine comparison of two algorithms only above that
       threshold, which is far beyond anything computed here. The test suite
       asserts the agreement and says plainly that this is what it means.

       This is the same trap as the ``mp.nzeros`` circularity documented in
       :mod:`critical_line.counting`: a comparison that looks independent
       because the two sides have different names.

    .. note::

       In exact arithmetic the product above is real. At finite precision it
       is not: the computed value carries a small imaginary component, which
       is discarded here. That discarded part is a **finite-precision
       residual**, not a mathematical quantity, and its size is a useful
       diagnostic of how much precision the evaluation actually retained.
       :func:`hardy_z_imaginary_residual` returns it.

    Parameters
    ----------
    t:
        Real argument.
    dps:
        Decimal working precision.
    """
    with working_precision(dps):
        t = mp.mpmathify(t)
        value = mp.e ** (1j * mp.siegeltheta(t)) * mp.zeta(mp.mpc(mp.mpf("0.5"), t))
        return mp.mpf(value.real)


def hardy_z_imaginary_residual(t: Scalar, dps: int = DEFAULT_DPS) -> "mp.mpf":
    """
    The magnitude of the imaginary part discarded by :func:`hardy_z_from_zeta`.

    Mathematically this is zero. Numerically it is not, and its size bounds
    how much confidence the computed real part deserves: a residual comparable
    to ``|Z(t)|`` itself means the evaluation has lost essentially all of its
    significant digits and the reported sign cannot be trusted.
    """
    with working_precision(dps):
        t = mp.mpmathify(t)
        value = mp.e ** (1j * mp.siegeltheta(t)) * mp.zeta(mp.mpc(mp.mpf("0.5"), t))
        return abs(mp.mpf(value.imag))

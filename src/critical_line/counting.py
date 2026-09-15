"""
Counting zeros three different ways, and keeping straight which comparisons
carry information.

Three quantities are in play. Conflating any two of them is the standard way
to believe you have verified more than you have.

**(A) N(T)** -- the number of zeros of zeta in the critical *strip* with
``0 < Im(s) <= T``. Says nothing about where in the strip they sit.

**(B)** The number of sign changes of Hardy's ``Z`` on ``(0, T]``. Each one
brackets a zero *on* the critical line, so this is a lower bound for the
number of on-line zeros.

**(C)** ``theta(T)/pi + 1`` -- the Riemann-von Mangoldt main term, a smooth
explicit function with no zero-finding in it at all.

The Riemann Hypothesis restricted to height ``T`` is the statement
**(A) == (B)**. That comparison means something only if (A) is computed
*without* already assuming where the zeros are.

Why ``mp.nzeros`` cannot serve as (A)
-------------------------------------

Reading mpmath's source (``mpmath/functions/zetazeros.py``), ``nzeros(t)`` is
built from ``gram_index(t) = floor(siegeltheta(t)/pi)`` plus Rosser-block
bookkeeping of ``Z(t)`` sign changes. Empirically, at every height tested::

    mp.nzeros(T) - (theta(T)/pi + 1) - backlunds(T) = 0   to working precision

(measured worst case over T in {50, 100, 200, 500, 1000}: 1.7e-24 at 25
digits, 7.3e-29 at 30, 6.6e-49 at 50 -- it tracks the precision, as an
identity must)

That is not two computations agreeing. Since ``S(T)`` is *defined* as
``N(T) - theta(T)/pi - 1``, the "discrepancy" such a comparison reports is
``S(T)`` by construction -- an identity rearranged, carrying no information
about the location of any zero. ``mp.zetazero`` calls the same Rosser-block
routines, so counting its outputs and comparing to ``nzeros`` is circular for
the same reason.

:func:`count_zeros_argument_principle` is a genuine (A): a contour integral of
``zeta'/zeta`` that never mentions theta, Gram points, ``Z``, or Rosser
blocks.

How far the independence goes
-----------------------------

Far enough to be worth having, not far enough to call the result certified.
(A) and (B) reach their answers by different routes -- a contour integral
versus sign changes along a line -- so a bug in one is unlikely to produce the
same wrong answer as a bug in the other.

But both bottom out in ``mp.zeta`` at finite precision with no enclosure. A
systematic error there, or a loss of precision severe enough to flip a sign,
could corrupt both counts at once and their agreement would not reveal it. The
procedures differ; the arithmetic underneath them does not. Hence the standing
wording in :data:`critical_line.constants.COUNTING_COMPARISON_NOTE`, which is
attached to every result record this module returns.
"""

from __future__ import annotations

from typing import Any, Callable, TypedDict

import mpmath as mp

from .constants import (
    COUNTING_COMPARISON_NOTE,
    DEFAULT_CONTOUR_DPS,
    DEFAULT_DPS,
    DEFAULT_GRID_STEP,
    DEFAULT_MAXDEGREE,
    DEFAULT_RESIDUAL_TOLERANCE,
    Scalar,
)
from .hardy_z import hardy_z
from .zeta_tools import log_derivative, working_precision

__all__ = [
    "ResidualToleranceError",
    "rvm_main_term",
    "rvm_main_term_stirling",
    "S",
    "count_zeros_argument_principle",
    "count_zeros_argument_principle_int",
    "count_sign_changes_hardy_z",
    "refine_sign_change",
    "experimentally_check_up_to",
    "ExperimentalCheck",
    "COUNTING_COMPARISON_NOTE",
]


class ResidualToleranceError(RuntimeError):
    """
    Raised when a contour-integral zero count is too far from an integer to be
    rounded responsibly.

    The argument principle returns an integer in exact arithmetic. A computed
    value that is not near an integer means the quadrature has not converged,
    and rounding it would manufacture a confident answer out of a failed
    computation. This package raises instead.

    Attributes
    ----------
    value:
        The raw complex value returned by the quadrature.
    residual:
        Distance from the nearest integer.
    tolerance:
        The tolerance that was not met.
    """

    def __init__(self, value: "mp.mpc", residual: "mp.mpf", tolerance: "mp.mpf") -> None:
        self.value = value
        self.residual = residual
        self.tolerance = tolerance
        super().__init__(
            f"argument-principle quadrature did not converge: value {mp.nstr(value, 12)} "
            f"is {mp.nstr(residual, 4)} from the nearest integer, exceeding the "
            f"tolerance {mp.nstr(tolerance, 4)}. Increase `maxdegree` or `dps`, or "
            f"move the contour away from a zero. Refusing to round an unconverged "
            f"result, since a small residual is only a numerical diagnostic and a "
            f"large one is evidence the count is wrong."
        )


# --------------------------------------------------------------------------
# (C) the smooth main term
# --------------------------------------------------------------------------

def rvm_main_term(T: Scalar, dps: int = DEFAULT_DPS) -> "mp.mpf":
    """
    The Riemann-von Mangoldt main term ``theta(T)/pi + 1``.

    This is a smooth explicit function. It is *not* a zero count and must not
    be compared against anything derived from theta -- see the module
    docstring.
    """
    with working_precision(dps):
        T = mp.mpmathify(T)
        return mp.siegeltheta(T) / mp.pi + 1


def rvm_main_term_stirling(T: Scalar, dps: int = DEFAULT_DPS) -> "mp.mpf":
    """
    The familiar textbook form of the main term,

    .. math:: \\frac{T}{2\\pi}\\log\\frac{T}{2\\pi} - \\frac{T}{2\\pi} + \\frac78,

    which is :func:`rvm_main_term` with Stirling's expansion applied to theta.
    The two agree to ``O(1/T)``, so this form should not be used as the main
    term at small ``T``: at ``T = 30`` they already differ by more than 0.01.
    Assumes ``T > 2 pi``; below that the logarithm's argument drops under 1 and
    the expansion is meaningless.
    """
    with working_precision(dps):
        T = mp.mpmathify(T)
        return (T / (2 * mp.pi)) * mp.log(T / (2 * mp.pi)) - T / (2 * mp.pi) + mp.mpf(7) / 8


def S(T: Scalar, dps: int = DEFAULT_DPS) -> "mp.mpf":
    """
    ``S(T) = N(T) - theta(T)/pi - 1``, the fluctuating part of the counting
    function, equal to ``(1/pi) arg zeta(1/2 + iT)``.

    Because S is *defined* as this difference, reporting it as the discrepancy
    between two independent computations is a category error. It is included
    here for reference and for the methodology demonstration in
    ``scripts/verify_core_claims.py``.
    """
    with working_precision(dps):
        return mp.backlunds(mp.mpmathify(T))


# --------------------------------------------------------------------------
# (A) an independent count: the argument principle
# --------------------------------------------------------------------------

def count_zeros_argument_principle(
    T: Scalar,
    *,
    sigma: Scalar = 2,
    y0: Scalar = "0.5",
    maxdegree: int = DEFAULT_MAXDEGREE,
    dps: int = DEFAULT_CONTOUR_DPS,
) -> "mp.mpc":
    """
    Count zeros of zeta in a rectangle by the argument principle.

    Evaluates

    .. math::

        \\frac{1}{2\\pi i} \\oint \\frac{\\zeta'(s)}{\\zeta(s)}\\, ds

    around the rectangle ``Re(s) in [1-sigma, sigma]``, ``Im(s) in [y0, T]``.

    With ``sigma = 2`` the rectangle contains the entire critical strip, so
    for ``0 < y0 < 14.13`` this counts exactly the nontrivial zeros with
    ``Im(s) in (0, T]``, with multiplicity, and **without any assumption about
    where in the strip they lie**. The pole at ``s = 1`` and every trivial
    zero sit on the real axis and are excluded by ``y0 > 0``.

    The integrand mentions only zeta and its derivative -- no theta, no Gram
    points, no Hardy Z, no Rosser blocks -- which is what makes this an
    independent count.

    .. warning::

       Choose ``T`` away from a zero ordinate. The upper edge of the contour
       runs through the strip at height ``T``; if a zero sits close to it,
       ``zeta`` is near zero on the path, the integrand is large, and the
       quadrature loses accuracy. At ``T = 150`` (where ``Z(150) = -0.09``,
       nearly a zero) the residual degrades by several orders of magnitude
       relative to nearby well-separated heights.

    .. note::

       Cost grows with ``T`` because the integrand oscillates more, so this is
       practical to heights of order hundreds, not millions. Large-scale
       verification uses Turing's method instead; see ``docs/certification.md``.

    Parameters
    ----------
    T:
        Upper height of the contour.
    sigma:
        Right edge of the rectangle; the left edge is ``1 - sigma``.
    y0:
        Lower edge. Must be positive and below the first zero ordinate.
    maxdegree:
        Quadrature refinement level passed to ``mp.quad``.
    dps:
        Decimal working precision.

    Returns
    -------
    mp.mpc
        The raw value, which should be a near-integer. Use
        :func:`count_zeros_argument_principle_int` for a checked integer.
    """
    with working_precision(dps):
        T = mp.mpmathify(T)
        y0 = mp.mpmathify(y0)
        sigma = mp.mpmathify(sigma)
        a, b = mp.mpf(1) - sigma, sigma

        corners = [
            mp.mpc(a, y0),
            mp.mpc(b, y0),
            mp.mpc(b, T),
            mp.mpc(a, T),
            mp.mpc(a, y0),
        ]
        total = mp.mpc(0)
        for p, q in zip(corners[:-1], corners[1:]):
            total += mp.quad(
                lambda u, p=p, q=q: log_derivative(p + u * (q - p)) * (q - p),
                [0, 1],
                maxdegree=maxdegree,
            )
        return total / (2j * mp.pi)


def count_zeros_argument_principle_int(
    T: Scalar,
    *,
    residual_tolerance: "mp.mpf" = mp.mpf(DEFAULT_RESIDUAL_TOLERANCE),
    **kwargs: Any,
) -> tuple[int, "mp.mpf"]:
    """
    The argument-principle count as a checked integer.

    In exact arithmetic the contour integral is an integer. This function
    rounds the computed value **only if** it is within ``residual_tolerance``
    of an integer, and raises :class:`ResidualToleranceError` otherwise rather
    than returning a confidently-wrong number.

    .. note::

       A small residual is a **numerical diagnostic, not a rigorous proof of
       convergence**. It says the quadrature is self-consistent at this
       refinement level. It does not bound the true error, because ``mp.quad``
       returns no enclosure. A *large* residual is much more informative than
       a small one: it is positive evidence that the count is wrong.

    Parameters
    ----------
    T:
        Upper height of the contour.
    residual_tolerance:
        Maximum accepted distance from the nearest integer. Keyword-only.
    **kwargs:
        Forwarded to :func:`count_zeros_argument_principle`
        (``sigma``, ``y0``, ``maxdegree``, ``dps``).

    Returns
    -------
    tuple
        ``(count, residual)``.

    Raises
    ------
    ResidualToleranceError
        If the residual exceeds ``residual_tolerance``.
    """
    residual_tolerance = mp.mpmathify(residual_tolerance)
    value = count_zeros_argument_principle(T, **kwargs)
    n = int(mp.nint(value.real))
    residual = abs(value - n)
    if residual > residual_tolerance:
        raise ResidualToleranceError(value, residual, residual_tolerance)
    return n, residual


# --------------------------------------------------------------------------
# (B) on-line zeros: sign changes of Hardy's Z
# --------------------------------------------------------------------------

def count_sign_changes_hardy_z(
    T: Scalar,
    *,
    t0: Scalar = "0.1",
    step: Scalar = DEFAULT_GRID_STEP,
    dps: int = DEFAULT_DPS,
    z_function: Callable[..., "mp.mpf"] | None = None,
) -> int:
    """
    Count sign changes of Hardy's ``Z`` on ``(t0, T]`` by uniform grid sampling.

    Each sign change brackets at least one zero of zeta on the critical line,
    so the result is a **lower bound** on the number of on-line zeros.

    .. warning::

       A uniform grid can step over a pair of zeros closer together than
       ``step`` -- the Lehmer phenomenon -- and a zero of even order would
       produce no sign change at all. Grid sampling alone therefore proves
       nothing. It becomes evidence only when paired with an independent total
       count; see :func:`experimentally_check_up_to`, where an undercount
       shows up as a disagreement rather than passing silently.

    Parameters
    ----------
    T:
        Upper limit.
    t0:
        Lower limit; must be positive and below the first zero ordinate.
    step:
        Grid spacing.
    dps:
        Decimal working precision.
    z_function:
        Override for the Z implementation, mainly for testing.
    """
    with working_precision(dps):
        # The precision context is entered once for the whole scan; at one grid
        # point every `step` the per-call context switch would dominate.
        if z_function is None:
            raw = mp.siegelz
        else:
            raw = lambda x: z_function(x, dps=dps)  # noqa: E731

        T = mp.mpmathify(T)
        step = mp.mpmathify(step)
        t = mp.mpmathify(t0)

        previous = raw(t)
        count = 0
        while t < T:
            t = min(t + step, T)
            current = raw(t)
            if previous == 0:
                count += 1
            elif previous * current < 0:
                count += 1
            previous = current
        return count


def refine_sign_change(
    lo: Scalar,
    hi: Scalar,
    *,
    dps: int = DEFAULT_DPS,
    tol: Scalar | None = None,
    coarse_dps: int | None = None,
    z_function: Callable[..., "mp.mpf"] | None = None,
) -> "mp.mpf":
    """
    Locate the zero ordinate inside a bracketed sign change of ``Z``.

    Two phases: cheap bisection at reduced precision to tighten the bracket,
    then a secant polish at full precision. See the inline comments for why
    the split matters.

    Parameters
    ----------
    lo, hi:
        Endpoints bracketing a sign change of ``Z``.
    dps:
        Decimal working precision for the final answer.
    tol:
        Absolute width at which to stop. Defaults to roughly the working
        precision.
    coarse_dps:
        Precision for the initial bisection. Defaults to ``min(dps, 25)``.
    z_function:
        Override for the Z implementation.

    Returns
    -------
    mp.mpf
        The located ordinate. This is where the bracket collapsed or where the
        solver converged -- not a proven enclosure. Use
        :func:`critical_line.zeros.zero_residual` to see what ``|zeta|``
        actually is there.
    """
    coarse_precision = min(dps, 25) if coarse_dps is None else coarse_dps
    tolerance = (
        mp.mpf(10) ** (-(dps - 3)) if tol is None else mp.mpmathify(tol)
    )

    # Phase 1: tighten the bracket at LOW precision. Only the sign of Z matters
    # here, and away from the zero |Z| is comfortably far from 0. This is worth
    # doing cheaply: at high dps mpmath's siegelz falls back to Euler-Maclaurin
    # summation, which costs an order of magnitude more per evaluation than the
    # same call at 25 digits.
    with working_precision(coarse_precision):
        if z_function is None:
            raw = mp.siegelz
        else:
            raw = lambda x: z_function(x, dps=coarse_precision)  # noqa: E731

        lo, hi = mp.mpmathify(lo), mp.mpmathify(hi)
        f_lo = raw(lo)
        # Only enough bisection to put the secant solver safely inside its
        # basin of convergence: from a grid step of 0.05 that is about six
        # evaluations. Bisecting further is wasted work, since the secant then
        # converges superlinearly and reaches full precision in roughly nine
        # more.
        coarse_target = max(tolerance, mp.mpf(10) ** -3)
        while hi - lo > coarse_target:
            mid = (lo + hi) / 2
            f_mid = raw(mid)
            if f_lo * f_mid <= 0:
                hi = mid
            else:
                lo, f_lo = mid, f_mid
        lo, hi = mp.mpf(lo), mp.mpf(hi)

    if hi - lo <= tolerance:
        return (lo + hi) / 2

    # Phase 2: polish at full precision. The secant method converges
    # superlinearly, so this costs a handful of expensive evaluations rather
    # than the ~3.3 * dps that bisecting all the way would need. The result is
    # accepted only if it stays inside the bracket known to contain a sign
    # change; otherwise fall back to bisection, which is slower but cannot
    # wander off.
    zfun = z_function if z_function is not None else hardy_z
    with working_precision(dps):
        try:
            root = mp.findroot(
                lambda t: zfun(t, dps=dps), (lo, hi), solver="secant",
            )
            if lo <= root <= hi:
                return mp.mpf(root)
        except (ValueError, ZeroDivisionError):
            pass

        f_lo = zfun(lo, dps=dps)
        while hi - lo > tolerance:
            mid = (lo + hi) / 2
            f_mid = zfun(mid, dps=dps)
            if f_lo * f_mid <= 0:
                hi = mid
            else:
                lo, f_lo = mid, f_mid
        return (lo + hi) / 2


# --------------------------------------------------------------------------
# putting (A) and (B) together
# --------------------------------------------------------------------------

class ExperimentalCheck(TypedDict):
    """Result record from :func:`experimentally_check_up_to`."""

    T: "mp.mpf"
    strip_count: int
    strip_residual: "mp.mpf"
    on_line_count: int
    agree: bool
    main_term: "mp.mpf"
    S: "mp.mpf"
    note: str


def experimentally_check_up_to(
    T: Scalar,
    *,
    step: Scalar = DEFAULT_GRID_STEP,
    maxdegree: int = DEFAULT_MAXDEGREE,
    residual_tolerance: "mp.mpf" = mp.mpf(DEFAULT_RESIDUAL_TOLERANCE),
    dps: int = DEFAULT_DPS,
    contour_dps: int = DEFAULT_CONTOUR_DPS,
    verbose: bool = False,
) -> ExperimentalCheck:
    """
    Compare an independent count of zeros in the strip against the number
    found on the critical line, both up to height ``T``.

    This function is deliberately **not** named ``certify_up_to``. This
    project computes with ``mpmath``, which does not provide interval-certified
    arithmetic, so "certify" would misdescribe what the output is. What this
    returns is experimental evidence at a stated precision.

    If the two counts agree, then -- subject to that caveat -- every zero of
    zeta with ``0 < Im(s) <= T`` lies on the critical line and is simple: the
    strip count bounds the total, the on-line count exhibits that many zeros
    on the line, and equality leaves no room for anything else.

    %s

    And it is silent about every zero above height ``T``. No finite
    computation can be otherwise; RH is a statement about infinitely many
    zeros.

    Parameters
    ----------
    T:
        Height up to which to check.
    step:
        Grid step for the Z scan.
    maxdegree:
        Quadrature refinement for the contour integral.
    residual_tolerance:
        Maximum accepted distance from an integer for the contour count.
    dps:
        Working precision for the Z scan.
    contour_dps:
        Working precision for the contour integral.
    verbose:
        Print a human-readable report.

    Returns
    -------
    ExperimentalCheck
        Record with both counts, the quadrature residual, whether they agree,
        and the main term and ``S(T)`` for reference.

    Raises
    ------
    ResidualToleranceError
        If the contour quadrature did not converge to within tolerance.
    """
    strip, residual = count_zeros_argument_principle_int(
        T,
        residual_tolerance=residual_tolerance,
        maxdegree=maxdegree,
        dps=contour_dps,
    )
    on_line = count_sign_changes_hardy_z(T, step=step, dps=dps)

    record: ExperimentalCheck = {
        "T": mp.mpmathify(T),
        "strip_count": strip,
        "strip_residual": residual,
        "on_line_count": on_line,
        "agree": strip == on_line,
        "main_term": rvm_main_term(T, dps=dps),
        "S": S(T, dps=dps),
        "note": COUNTING_COMPARISON_NOTE,
    }

    if verbose:
        print(f"  T = {float(record['T']):g}")
        print(f"    zeros in strip (argument principle) : {strip}")
        print(f"      quadrature residual              : {mp.nstr(residual, 3)}"
              f"   (numerical diagnostic, not a proof of convergence)")
        print(f"    zeros on the line (Z sign changes)  : {on_line}")
        print(f"    counts agree                        : {record['agree']}")
        print(f"    reference: theta(T)/pi + 1 = {float(record['main_term']):.4f}, "
              f"S(T) = {float(record['S']):+.4f}")

    return record


experimentally_check_up_to.__doc__ = (
    experimentally_check_up_to.__doc__ or ""
) % COUNTING_COMPARISON_NOTE

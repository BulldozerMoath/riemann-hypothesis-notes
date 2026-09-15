"""
Shared constants, default settings, and the project's standing disclaimer.

Keeping the disclaimer text in one place means every figure, script, and
report carries exactly the same wording, and that wording cannot drift
independently across the project.
"""

from __future__ import annotations

from typing import Final, Union

import mpmath as mp

# Anything that can be handed to mp.mpmathify.
Scalar = Union[int, float, str, "mp.mpf"]

# --------------------------------------------------------------------------
# defaults
# --------------------------------------------------------------------------

#: Default working precision, in decimal digits.
DEFAULT_DPS: Final[int] = 50

#: Default precision for contour integration, where the cost of extra digits
#: is paid on every quadrature node.
DEFAULT_CONTOUR_DPS: Final[int] = 30

#: Default quadrature refinement level for the argument-principle contour.
#: Empirically, ``maxdegree=6`` is not enough past T ~ 40; ``8`` holds to
#: T ~ 120 at :data:`DEFAULT_CONTOUR_DPS`. See docs/numerical_methods.md.
DEFAULT_MAXDEGREE: Final[int] = 8

#: Default residual tolerance for accepting a rounded contour-integral count.
#: A residual below this is a *numerical diagnostic*, not a proof of
#: convergence.
DEFAULT_RESIDUAL_TOLERANCE: Final[str] = "1e-8"

#: Default grid step for scanning sign changes of Hardy's Z.
DEFAULT_GRID_STEP: Final[str] = "0.05"

#: Default threshold on the normalized gap below which a consecutive pair of
#: zeros is reported as Lehmer-like.
DEFAULT_LEHMER_THRESHOLD: Final[str] = "0.15"

# --------------------------------------------------------------------------
# disclaimers
# --------------------------------------------------------------------------

#: Short form, for figure captions and one-line script banners.
DISCLAIMER_SHORT: Final[str] = (
    "Numerical exploration only — this figure does not prove the "
    "Riemann Hypothesis."
)

#: Long form, for script headers and report footers.
DISCLAIMER_LONG: Final[str] = (
    "This project performs finite-precision numerical exploration. Its "
    "outputs are experimental evidence, not certified results, and nothing "
    "it computes bears on the truth of the Riemann Hypothesis, which "
    "remains open."
)

#: Standing wording for why two differing counting procedures agreeing is
#: evidence rather than proof. Required verbatim wherever the comparison is
#: described.
COUNTING_COMPARISON_NOTE: Final[str] = (
    "These use different numerical counting procedures, but both rely on "
    "finite-precision numerical evaluation and therefore provide "
    "experimental evidence rather than a certified result."
)

#: Standing wording about Lehmer pairs.
LEHMER_NOTE: Final[str] = (
    "Lehmer pairs are unusually close consecutive zeros on the critical "
    "line. They are not counterexamples to the Riemann Hypothesis, but they "
    "make naive zero searches more difficult."
)

# --------------------------------------------------------------------------
# reference data
# --------------------------------------------------------------------------

#: Imaginary parts of the first ten nontrivial zeros, to 15 decimal places.
#: Source: NIST DLMF §25.10 and standard tables; independently reproduced by
#: this project's grid-and-bisection finder (see tests/test_zeros.py).
FIRST_ZERO_ORDINATES: Final[tuple[str, ...]] = (
    "14.134725141734693",
    "21.022039638771555",
    "25.010857580145688",
    "30.424876125859513",
    "32.935061587739190",
    "37.586178158825671",
    "40.918719012147495",
    "43.327073280914999",
    "48.005150881167159",
    "49.773832477672302",
)

#: The smallest n at which Gram's law is known to fail.
FIRST_GRAM_LAW_FAILURE: Final[int] = 126


def wigner_surmise_gue(s: Scalar) -> "mp.mpf":
    """
    The GUE Wigner surmise density for normalized spacings,

    .. math:: p(s) = \\frac{32}{\\pi^2} s^2 \\exp(-4 s^2 / \\pi).

    Used only as a visual overlay when plotting observed zero spacings.
    Agreement is an empirical numerical comparison; Montgomery's pair
    correlation conjecture remains a conjecture.

    Assumes ``s`` is an already-unfolded (unit-mean) spacing; feeding it raw
    gaps is meaningless. This is the Wigner *surmise*, a 2x2 approximation to
    the true GUE spacing density -- close, but not the exact distribution, so
    small deviations in a histogram overlay are expected and are not evidence
    of anything.
    """
    s = mp.mpmathify(s)
    return (32 / (mp.pi ** 2)) * s ** 2 * mp.exp(-4 * s ** 2 / mp.pi)


def poisson_density(s: Scalar) -> "mp.mpf":
    """
    Density ``exp(-s)`` of gaps in a unit-rate Poisson process, plotted for
    contrast with :func:`wigner_surmise_gue`.

    Assumes ``s`` is an unfolded (unit-mean) spacing. This is what *uncorrelated*
    points would give; the zeros visibly do not follow it. That visible
    difference is an observation about computed zeros in a finite range, not a
    theorem about all of them.
    """
    return mp.exp(-mp.mpmathify(s))

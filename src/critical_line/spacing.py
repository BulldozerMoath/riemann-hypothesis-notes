"""
Spacing statistics of zeros on the critical line.

Consecutive zero ordinates get closer together as height increases, at a rate
given by the density of zeros. To compare spacings across heights they must be
*unfolded* -- rescaled by the local mean spacing

.. math:: \\frac{2\\pi}{\\log(\\gamma / 2\\pi)}

so the rescaled gaps have mean 1 by construction.

The unfolded gaps are then compared against the Gaussian Unitary Ensemble
prediction from random matrix theory. **Every such comparison in this project
is an empirical numerical comparison, not a proof of anything.** Montgomery's
pair correlation conjecture is a conjecture; the striking agreement Odlyzko
found numerically is evidence for it, and evidence is not proof. Nor would
proving it prove RH.
"""

from __future__ import annotations

from typing import Sequence, TypedDict

import mpmath as mp

from .constants import DEFAULT_LEHMER_THRESHOLD, LEHMER_NOTE, Scalar

__all__ = [
    "zero_gaps",
    "normalized_gaps",
    "spacing_summary",
    "lehmer_like_pairs",
    "local_mean_spacing",
    "SpacingSummary",
    "LehmerPair",
    "LEHMER_NOTE",
]


class SpacingSummary(TypedDict):
    """Summary statistics returned by :func:`spacing_summary`."""

    n_zeros: int
    n_gaps: int
    height_min: "mp.mpf"
    height_max: "mp.mpf"
    raw_gap_min: "mp.mpf"
    raw_gap_max: "mp.mpf"
    normalized_mean: "mp.mpf"
    normalized_min: "mp.mpf"
    normalized_max: "mp.mpf"
    normalized_std: "mp.mpf"


class LehmerPair(TypedDict):
    """One unusually close consecutive pair, from :func:`lehmer_like_pairs`."""

    index: int
    gamma_lower: "mp.mpf"
    gamma_upper: "mp.mpf"
    raw_gap: "mp.mpf"
    normalized_gap: "mp.mpf"


def local_mean_spacing(gamma: Scalar) -> "mp.mpf":
    """
    The local mean spacing between zeros at height ``gamma``,

    .. math:: \\frac{2\\pi}{\\log(\\gamma / 2\\pi)}.

    This is the reciprocal of the zero density, obtained by differentiating
    the Riemann-von Mangoldt main term.

    Assumes ``gamma > 2 pi`` (the logarithm is negative below that and the
    "spacing" it returns is meaningless) and is an *asymptotic* density: it
    describes the average local rate, so individual gaps scatter widely around
    it. Using it to normalize is a convention, not a measurement -- see
    :func:`normalized_gaps`.
    """
    gamma = mp.mpmathify(gamma)
    return 2 * mp.pi / mp.log(gamma / (2 * mp.pi))


def zero_gaps(zeros: Sequence[Scalar]) -> list["mp.mpf"]:
    """
    Raw differences between consecutive ordinates.

    Parameters
    ----------
    zeros:
        Ordinates in increasing order.

    Assumes the input is sorted in increasing order; this is not checked, and
    unsorted input silently yields negative gaps. Gaps are only as accurate as
    the ordinates they are formed from, and subtracting two nearby ordinates
    cancels leading digits -- a gap carries noticeably fewer significant
    figures than its endpoints do.

    Returns
    -------
    list of mp.mpf
        Length ``len(zeros) - 1``; empty if fewer than two zeros given.
    """
    vals = [mp.mpmathify(z) for z in zeros]
    return [b - a for a, b in zip(vals[:-1], vals[1:])]


def normalized_gaps(zeros: Sequence[Scalar]) -> list["mp.mpf"]:
    """
    Gaps rescaled by the local mean spacing, so the mean is 1 by construction.

    .. math::

        \\delta_n = (\\gamma_{n+1} - \\gamma_n)
                    \\cdot \\frac{\\log(\\gamma_n / 2\\pi)}{2\\pi}

    A mean near 1 is therefore *not* a check on anything -- it is what the
    normalization was designed to produce. What carries information is the
    *shape* of the distribution, which is what
    :func:`critical_line.plots.plot_normalized_gap_histogram` displays.

    Parameters
    ----------
    zeros:
        Ordinates in increasing order.
    """
    vals = [mp.mpmathify(z) for z in zeros]
    return [
        (b - a) / local_mean_spacing(a)
        for a, b in zip(vals[:-1], vals[1:])
    ]


def spacing_summary(zeros: Sequence[Scalar]) -> SpacingSummary:
    """
    Summary statistics of raw and normalized gaps.

    Parameters
    ----------
    zeros:
        Ordinates in increasing order; at least two required.

    Returns
    -------
    SpacingSummary

    Raises
    ------
    ValueError
        If fewer than two ordinates are supplied, since no gap exists.
    """
    vals = [mp.mpmathify(z) for z in zeros]
    if len(vals) < 2:
        raise ValueError(
            f"need at least 2 zeros to form a gap, got {len(vals)}"
        )

    raw = zero_gaps(vals)
    norm = normalized_gaps(vals)
    mean = sum(norm) / len(norm)
    variance = sum((g - mean) ** 2 for g in norm) / len(norm)

    return {
        "n_zeros": len(vals),
        "n_gaps": len(raw),
        "height_min": vals[0],
        "height_max": vals[-1],
        "raw_gap_min": min(raw),
        "raw_gap_max": max(raw),
        "normalized_mean": mean,
        "normalized_min": min(norm),
        "normalized_max": max(norm),
        "normalized_std": mp.sqrt(variance),
    }


def lehmer_like_pairs(
    zeros: Sequence[Scalar],
    threshold: "mp.mpf" = mp.mpf(DEFAULT_LEHMER_THRESHOLD),
) -> list[LehmerPair]:
    """
    Consecutive pairs whose normalized gap falls below ``threshold``.

    Lehmer pairs are unusually close consecutive zeros on the critical line.
    They are not counterexamples to the Riemann Hypothesis, but they make
    naive zero searches more difficult: a grid coarser than the gap will step
    straight over both zeros and report neither, and the resulting undercount
    is silent unless an independent total count is available to contradict it.

    They are also where the de Bruijn-Newman constant result bites. Rodgers
    and Tao proved ``Lambda >= 0`` in 2020; RH is equivalent to
    ``Lambda <= 0``. Close pairs are the configurations that would have to
    exist in profusion for ``Lambda`` to be strictly positive, which is why
    they get attention out of proportion to their number.

    Parameters
    ----------
    zeros:
        Ordinates in increasing order.
    threshold:
        Normalized-gap cutoff. The default ``0.15`` is a convention, not a
        definition; Lehmer's original 1956 example near ``t = 7005`` has a
        normalized gap well under ``0.1``.

    Returns
    -------
    list of LehmerPair
        Possibly empty. Low-lying zeros are well separated, so an empty result
        at modest heights is expected.
    """
    threshold = mp.mpmathify(threshold)
    vals = [mp.mpmathify(z) for z in zeros]
    raw = zero_gaps(vals)
    norm = normalized_gaps(vals)

    return [
        {
            "index": i,
            "gamma_lower": vals[i],
            "gamma_upper": vals[i + 1],
            "raw_gap": raw[i],
            "normalized_gap": norm[i],
        }
        for i in range(len(norm))
        if norm[i] < threshold
    ]

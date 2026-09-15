"""
Figure generation.

Every figure produced here carries a visible disclaimer annotation. That is
deliberate: figures get screenshotted, pasted into slides, and separated from
their captions, and a plot of zeros on the critical line is exactly the kind
of image that gets mistaken for evidence of something it does not show.

Palette
-------
Three categorical slots, validated for colour-vision deficiency separation and
contrast against the light chart surface:

==========  =========  ==============================================
slot        hex        role
==========  =========  ==============================================
series 1    ``#2a78d6``  observed / computed data
series 2    ``#eb6834``  theoretical reference being compared against
series 3    ``#1baf7a``  secondary reference (dashed; see note)
==========  =========  ==============================================

Slot 3 sits below 3:1 contrast on the light surface, so it is used only for
curves that also carry a legend label and a distinct dash pattern -- identity
never rests on hue alone.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Sequence

import matplotlib

matplotlib.use("Agg")  # headless; must precede pyplot import

import matplotlib.pyplot as plt  # noqa: E402
import mpmath as mp  # noqa: E402
import numpy as np  # noqa: E402

from .constants import (  # noqa: E402
    DEFAULT_DPS,
    DISCLAIMER_SHORT,
    Scalar,
    poisson_density,
    wigner_surmise_gue,
)
from .counting import rvm_main_term  # noqa: E402
from .explicit_formula import explicit_formula_partial_sums, psi  # noqa: E402
from .hardy_z import hardy_z  # noqa: E402
from .spacing import normalized_gaps  # noqa: E402

__all__ = [
    "plot_hardy_z",
    "plot_zero_ordinates",
    "plot_normalized_gap_histogram",
    "plot_counting_comparison",
    "plot_explicit_formula_error",
]

# --------------------------------------------------------------------------
# palette (validated; see module docstring)
# --------------------------------------------------------------------------

SURFACE = "#fcfcfb"
INK_PRIMARY = "#0b0b0b"
INK_SECONDARY = "#52514e"
INK_MUTED = "#898781"
GRIDLINE = "#e1e0d9"
AXIS = "#c3c2b7"

SERIES_1 = "#2a78d6"
SERIES_2 = "#eb6834"
SERIES_3 = "#1baf7a"

_LINEWIDTH = 1.8
_MARKERSIZE = 5.0


def _new_figure(figsize: tuple[float, float] = (10.0, 5.6)):
    """Create a figure and axes wearing the project's chart chrome."""
    fig, ax = plt.subplots(figsize=figsize, dpi=150)
    fig.patch.set_facecolor(SURFACE)
    ax.set_facecolor(SURFACE)

    ax.grid(True, color=GRIDLINE, linewidth=0.8, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(AXIS)
        ax.spines[side].set_linewidth(1.0)
    ax.tick_params(colors=INK_MUTED, labelsize=9, length=0)
    return fig, ax


def _finish(fig, ax, title: str, xlabel: str, ylabel: str, output_path) -> Path:
    """Apply titles, the standing disclaimer, and write the PNG."""
    ax.set_title(title, color=INK_PRIMARY, fontsize=13, pad=14, loc="left")
    ax.set_xlabel(xlabel, color=INK_SECONDARY, fontsize=10)
    ax.set_ylabel(ylabel, color=INK_SECONDARY, fontsize=10)

    fig.text(
        0.5, 0.012, DISCLAIMER_SHORT,
        ha="center", va="bottom", color=INK_MUTED, fontsize=8.5, style="italic",
    )
    fig.tight_layout(rect=(0, 0.045, 1, 1))

    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, facecolor=fig.get_facecolor(), bbox_inches="tight")
    plt.close(fig)
    return path


def _style_legend(ax) -> None:
    legend = ax.legend(
        frameon=True, fontsize=9, loc="best",
        facecolor=SURFACE, edgecolor=GRIDLINE, framealpha=0.95,
    )
    for text in legend.get_texts():
        text.set_color(INK_SECONDARY)


# --------------------------------------------------------------------------
# figures
# --------------------------------------------------------------------------

def plot_hardy_z(
    start: Scalar,
    end: Scalar,
    output_path: str | os.PathLike[str],
    dps: int = DEFAULT_DPS,
    *,
    samples: int = 900,
) -> Path:
    """
    Plot Hardy's ``Z(t)`` over ``[start, end]``, marking its sign changes.

    Every crossing of the horizontal axis is a zero of zeta *on the critical
    line*. Z cannot see zeros off the line, so a picture full of crossings is
    not evidence for RH -- it is a picture of the thing RH is about, which is
    not the same claim.

    Parameters
    ----------
    start, end:
        Range of ``t``.
    output_path:
        Destination PNG path; parent directories are created.
    dps:
        Decimal working precision.
    samples:
        Number of evaluation points.

    Returns
    -------
    pathlib.Path
        The path written.
    """
    a, b = float(mp.mpmathify(start)), float(mp.mpmathify(end))
    ts = np.linspace(a, b, samples)
    zs = np.array([float(hardy_z(t, dps=dps)) for t in ts])

    crossings = [
        float(ts[i] - zs[i] * (ts[i + 1] - ts[i]) / (zs[i + 1] - zs[i]))
        for i in range(len(ts) - 1)
        if zs[i] * zs[i + 1] < 0
    ]

    fig, ax = _new_figure()
    ax.axhline(0, color=AXIS, linewidth=1.0, zorder=1)
    ax.plot(ts, zs, color=SERIES_1, linewidth=_LINEWIDTH, zorder=3,
            label="$Z(t)$")
    if crossings:
        ax.plot(
            crossings, [0.0] * len(crossings), "o",
            color=SERIES_2, markersize=_MARKERSIZE, zorder=4,
            markeredgecolor=SURFACE, markeredgewidth=1.2,
            label=f"sign changes ({len(crossings)})",
        )
    _style_legend(ax)
    return _finish(
        fig, ax,
        f"Hardy's $Z(t)$ on $[{a:g}, {b:g}]$ — every crossing is a zero on the critical line",
        "$t$", "$Z(t)$", output_path,
    )


def plot_zero_ordinates(
    zeros: Sequence[Scalar],
    output_path: str | os.PathLike[str],
) -> Path:
    """
    Plot the located zero ordinates as a spectrum of vertical ticks.

    This is the "bar code" view: it shows how the zeros thin out slowly with
    height, and how irregular the local spacing is.

    Parameters
    ----------
    zeros:
        Ordinates in increasing order.
    output_path:
        Destination PNG path.

    Raises
    ------
    ValueError
        If no zeros are supplied.
    """
    vals = [float(mp.mpmathify(z)) for z in zeros]
    if not vals:
        raise ValueError("no zeros supplied; nothing to plot")

    fig, ax = _new_figure(figsize=(10.0, 3.4))
    ax.eventplot(
        vals, orientation="horizontal", colors=SERIES_1,
        lineoffsets=0.5, linelengths=0.8, linewidths=1.1, zorder=3,
    )
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.grid(axis="y", visible=False)
    ax.annotate(
        f"{len(vals)} zeros,  $\\gamma \\in [{vals[0]:.2f},\\ {vals[-1]:.2f}]$",
        xy=(0.995, 0.93), xycoords="axes fraction", ha="right",
        color=INK_SECONDARY, fontsize=9,
    )
    return _finish(
        fig, ax,
        "Zero ordinates on the critical line",
        r"$\gamma$  (imaginary part of $\rho = \frac{1}{2} + i\gamma$)", "",
        output_path,
    )


def plot_normalized_gap_histogram(
    zeros: Sequence[Scalar],
    output_path: str | os.PathLike[str],
    *,
    bins: int = 24,
    xmax: float = 3.0,
) -> Path:
    """
    Histogram of normalized consecutive gaps, against GUE and Poisson curves.

    The GUE (Gaussian Unitary Ensemble) curve is the Wigner surmise from
    random matrix theory; the Poisson curve is what uncorrelated points would
    give. Observed gaps follow GUE closely and Poisson not at all -- the zeros
    repel each other.

    This is an **empirical numerical comparison, not a proof**. Montgomery's
    pair correlation conjecture remains a conjecture, and proving it would not
    prove RH.

    Parameters
    ----------
    zeros:
        Ordinates in increasing order; at least three required.
    output_path:
        Destination PNG path.
    bins:
        Histogram bin count.
    xmax:
        Right edge of the plotted range.

    Raises
    ------
    ValueError
        If fewer than three ordinates are supplied.
    """
    vals = list(zeros)
    if len(vals) < 3:
        raise ValueError(f"need at least 3 zeros to form a gap distribution, got {len(vals)}")

    gaps = [float(g) for g in normalized_gaps(vals)]
    grid = np.linspace(1e-4, xmax, 400)
    gue = np.array([float(wigner_surmise_gue(s)) for s in grid])
    poisson = np.array([float(poisson_density(s)) for s in grid])

    fig, ax = _new_figure()
    ax.hist(
        gaps, bins=bins, range=(0.0, xmax), density=True,
        color=SERIES_1, alpha=0.85, zorder=2,
        edgecolor=SURFACE, linewidth=1.2,
        label=f"observed gaps (n = {len(gaps)})",
    )
    ax.plot(grid, gue, color=SERIES_2, linewidth=_LINEWIDTH, zorder=4,
            label="GUE (Wigner surmise)")
    ax.plot(grid, poisson, color=SERIES_3, linewidth=_LINEWIDTH,
            linestyle=(0, (5, 4)), zorder=3, label="Poisson (uncorrelated)")
    ax.set_xlim(0, xmax)
    _style_legend(ax)
    return _finish(
        fig, ax,
        "Normalized gaps between consecutive zeros — level repulsion",
        r"normalized gap $\delta$", "density", output_path,
    )


def plot_counting_comparison(
    t_values: Sequence[Scalar],
    output_path: str | os.PathLike[str],
    dps: int = DEFAULT_DPS,
    *,
    zeros: Sequence[Scalar] | None = None,
) -> Path:
    """
    Plot the observed on-line zero count against the Riemann-von Mangoldt
    main term ``theta(T)/pi + 1``.

    .. warning::

       This figure is **illustrative, not a verification**. The step function
       counts zeros found on the critical line; the smooth curve is the main
       term of the counting formula. Their difference is essentially ``S(T)``,
       which is *defined* as that difference. Close agreement here is built
       in, not discovered. The comparison that carries information is
       :func:`critical_line.counting.experimentally_check_up_to`, which puts
       an independent strip count beside the on-line count.

    Parameters
    ----------
    t_values:
        Heights at which to evaluate the main term.
    output_path:
        Destination PNG path.
    dps:
        Decimal working precision.
    zeros:
        Ordinates for the step function. If ``None``, they are located by
        scanning up to ``max(t_values)``.
    """
    ts = [float(mp.mpmathify(t)) for t in t_values]
    if not ts:
        raise ValueError("t_values is empty")
    top = max(ts)

    if zeros is None:
        from .zeros import find_on_line_zeros
        zeros = find_on_line_zeros(top, dps=dps)
    gammas = [float(mp.mpmathify(g)) for g in zeros]

    main = [float(rvm_main_term(t, dps=dps)) for t in ts]

    fig, ax = _new_figure()
    if gammas:
        step_x = [0.0]
        step_y = [0.0]
        for i, g in enumerate(gammas, start=1):
            step_x += [g, g]
            step_y += [i - 1, i]
        step_x.append(top)
        step_y.append(len(gammas))
        ax.plot(step_x, step_y, color=SERIES_1, linewidth=_LINEWIDTH, zorder=3,
                label="zeros found on the critical line")
    ax.plot(ts, main, color=SERIES_2, linewidth=_LINEWIDTH,
            linestyle=(0, (6, 4)), zorder=4,
            label=r"$\theta(T)/\pi + 1$  (main term)")
    ax.set_xlim(0, top)
    _style_legend(ax)
    return _finish(
        fig, ax,
        "Counting zeros: observed on-line count vs the smooth main term",
        "$T$", "number of zeros with $0 < \\gamma \\leq T$", output_path,
    )


def plot_explicit_formula_error(
    x: Scalar,
    zeros: Sequence[Scalar],
    output_path: str | os.PathLike[str],
) -> Path:
    """
    Plot partial sums of the explicit formula converging toward ``psi(x)``.

    Each step adds one conjugate pair of zeros. The oscillation is the point:
    the zeros do not approximate psi smoothly, they encode its jumps, and the
    partial sums swing around the true value while closing in.

    Parameters
    ----------
    x:
        Evaluation point. Avoid prime powers, where psi jumps and the formula
        converges to the midpoint.
    zeros:
        Ordinates in increasing order.
    output_path:
        Destination PNG path.

    Raises
    ------
    ValueError
        If no zeros are supplied.
    """
    vals = list(zeros)
    if not vals:
        raise ValueError("no zeros supplied; nothing to plot")

    sums = [float(s) for s in explicit_formula_partial_sums(x, vals)]
    truth = float(psi(x))
    ks = list(range(len(sums)))

    fig, ax = _new_figure()
    ax.plot(ks, sums, color=SERIES_1, linewidth=_LINEWIDTH, zorder=3,
            label="partial sum over zero pairs")
    ax.axhline(truth, color=SERIES_2, linewidth=_LINEWIDTH,
               linestyle=(0, (6, 4)), zorder=4,
               label=rf"$\psi({float(mp.mpmathify(x)):g}) = {truth:.6f}$")
    ax.set_xlim(0, len(sums) - 1)
    _style_legend(ax)
    return _finish(
        fig, ax,
        rf"The explicit formula: rebuilding $\psi({float(mp.mpmathify(x)):g})$ from zeros of $\zeta$",
        "number of conjugate zero pairs summed", r"$\psi(x)$ approximation",
        output_path,
    )

"""
The Critical Line — educational and numerical tooling for the Riemann zeta
function and the Riemann Hypothesis.

This package performs finite-precision numerical exploration. Its outputs are
experimental evidence, not certified results, and nothing it computes bears on
the truth of the Riemann Hypothesis, which remains open.

No function in this package is named ``certify``. ``mpmath`` does not provide
interval-certified arithmetic, so that word would misdescribe what is
produced. The comparison routine is
:func:`critical_line.counting.experimentally_check_up_to`.

Modules
-------
``constants``
    Defaults, disclaimer text, reference data.
``zeta_tools``
    zeta, its logarithmic derivative, the functional equation, xi.
``hardy_z``
    Riemann-Siegel theta and Hardy's Z function.
``zeros``
    Locating zeros on the critical line by grid scan and bisection.
``counting``
    Three ways to count zeros, and which comparisons carry information.
``gram``
    Gram points and Gram's law, as a diagnostic — not a proof technique.
``spacing``
    Gap statistics and GUE comparison.
``explicit_formula``
    Von Mangoldt's formula: primes rebuilt from zeros.
``plots``
    Figure generation, every figure carrying the standing disclaimer.
"""

from __future__ import annotations

from . import (
    constants,
    counting,
    explicit_formula,
    gram,
    hardy_z,
    spacing,
    zeros,
    zeta_tools,
)
from .constants import DISCLAIMER_LONG, DISCLAIMER_SHORT

__version__ = "0.5.0"

__all__ = [
    "constants",
    "counting",
    "explicit_formula",
    "gram",
    "hardy_z",
    "spacing",
    "zeros",
    "zeta_tools",
    "plots",
    "DISCLAIMER_SHORT",
    "DISCLAIMER_LONG",
    "__version__",
]


def __getattr__(name: str):
    """
    Import :mod:`critical_line.plots` lazily.

    It pulls in matplotlib, which is slow to import and unnecessary for the
    purely numerical parts of the package.

    Uses :func:`importlib.import_module` rather than ``from . import plots``:
    the latter re-enters this function through the import system's attribute
    lookup and recurses until the stack is exhausted.
    """
    if name == "plots":
        import importlib

        module = importlib.import_module(".plots", __name__)
        globals()["plots"] = module  # cache, so this runs only once
        return module
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

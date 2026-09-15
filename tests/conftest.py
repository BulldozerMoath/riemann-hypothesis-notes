"""Shared fixtures. Precision is set per-test so no test leaks a setting."""

from __future__ import annotations

import mpmath as mp
import pytest

from critical_line.zeta_tools import set_precision

#: Ordinates of the first five zeros, to 15 decimals (NIST DLMF 25.10).
KNOWN_ORDINATES = [
    "14.134725141734693",
    "21.022039638771555",
    "25.010857580145688",
    "30.424876125859513",
    "32.935061587739190",
]


@pytest.fixture(autouse=True)
def _restore_precision():
    """Restore mpmath's global precision after every test."""
    old = mp.mp.dps
    yield
    set_precision(old)


@pytest.fixture(scope="session")
def zeros_to_60():
    """Zeros on the critical line below T = 60, computed once per session."""
    from critical_line.zeros import find_on_line_zeros
    return find_on_line_zeros(60, dps=25)


@pytest.fixture(scope="session")
def zeros_to_150():
    """Zeros below T = 150, computed once per session."""
    from critical_line.zeros import find_on_line_zeros
    return find_on_line_zeros(150, dps=25)

"""Tests for the grid-and-bisection zero finder."""

from __future__ import annotations

import mpmath as mp
import pytest

from conftest import KNOWN_ORDINATES
from critical_line import zeros as zr


def test_find_zeros_matches_known_ordinates():
    found = zr.find_on_line_zeros(35, dps=25)
    assert len(found) == len(KNOWN_ORDINATES)
    for got, want in zip(found, KNOWN_ORDINATES):
        assert abs(got - mp.mpf(want)) < mp.mpf("1e-12")


def test_find_zeros_agrees_with_mpmath_zetazero(zeros_to_60):
    """
    Second opinion: grid+bisection against mpmath's Rosser-block machinery.
    These are genuinely different algorithms, so agreement is meaningful --
    unlike nzeros vs the counting formula, which is circular.
    """
    for k, got in enumerate(zeros_to_60, start=1):
        assert abs(got - mp.zetazero(k).imag) < mp.mpf("1e-12")


def test_found_zeros_are_increasing(zeros_to_60):
    assert all(b > a for a, b in zip(zeros_to_60[:-1], zeros_to_60[1:]))


def test_zero_residual_is_small(zeros_to_60):
    """|zeta(1/2 + i gamma)| at a located ordinate: small, never exactly zero."""
    for gamma in zeros_to_60[:5]:
        residual = zr.zero_residual(gamma, dps=25)
        assert residual >= 0
        assert residual < mp.mpf("1e-18")


def test_first_n_zeros_returns_exact_count():
    for n in [1, 5, 12]:
        assert len(zr.first_n_zeros(n, dps=25)) == n


def test_first_n_zeros_extends_search_range():
    """More zeros than fit below the initial height forces the range to grow."""
    found = zr.first_n_zeros(12, dps=25, initial_height=30)
    assert len(found) == 12
    assert found[-1] > 30


def test_first_n_zeros_empty_for_nonpositive():
    assert zr.first_n_zeros(0, dps=25) == []


def test_zeros_below_first_ordinate_is_empty():
    assert zr.find_on_line_zeros(14, dps=25) == []

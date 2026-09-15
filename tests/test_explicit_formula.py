"""Tests for the von Mangoldt sieve, psi, and the explicit formula."""

from __future__ import annotations

import math

import mpmath as mp
import pytest

from critical_line import explicit_formula as ef


def test_von_mangoldt_values():
    lam = ef.von_mangoldt_sieve(20)
    assert abs(lam[2] - math.log(2)) < 1e-12
    assert abs(lam[4] - math.log(2)) < 1e-12    # 2^2
    assert abs(lam[8] - math.log(2)) < 1e-12    # 2^3
    assert abs(lam[9] - math.log(3)) < 1e-12    # 3^2
    assert lam[6] == 0.0                        # 2 * 3, not a prime power
    assert lam[1] == 0.0
    assert lam[0] == 0.0


def test_von_mangoldt_rejects_negative():
    with pytest.raises(ValueError, match="non-negative"):
        ef.von_mangoldt_sieve(-1)


def test_psi_small_value():
    """psi(10) = 3 log 2 + 2 log 3 + log 5 + log 7."""
    expected = 3 * math.log(2) + 2 * math.log(3) + math.log(5) + math.log(7)
    assert abs(ef.psi(10) - expected) < 1e-12


def test_psi_below_two_is_zero():
    assert ef.psi(1) == 0.0
    assert ef.psi(0) == 0.0


def test_psi_is_nondecreasing():
    values = [ef.psi(x) for x in [10, 20, 50, 100]]
    assert all(b >= a for a, b in zip(values[:-1], values[1:]))


def test_baseline_overshoots_psi():
    """With no zeros summed the formula gives roughly x, above psi(x)."""
    x = 100.0
    baseline = ef.explicit_formula_partial_sums(x, [], dps=25)[0]
    assert baseline > ef.psi(x)
    assert abs(baseline - x) < 5


def test_partial_sums_length(zeros_to_150):
    sums = ef.explicit_formula_partial_sums(100.0, zeros_to_150, dps=25)
    assert len(sums) == len(zeros_to_150) + 1


def test_explicit_formula_converges(zeros_to_150):
    """
    More zeros means a better approximation on average. The partial sums
    oscillate, so early and late errors are compared in aggregate rather than
    demanding monotonicity.
    """
    x = 100.0
    sums = ef.explicit_formula_partial_sums(x, zeros_to_150, dps=25)
    truth = mp.mpf(ef.psi(x))
    errors = [abs(s - truth) for s in sums]
    early = sum(errors[1:11]) / 10
    late = sum(errors[-10:]) / 10
    assert late < early


def test_explicit_formula_error_tuple(zeros_to_150):
    approximation, truth, error = ef.explicit_formula_error(100.0, zeros_to_150, dps=25)
    assert abs(abs(approximation - truth) - error) < mp.mpf("1e-20")
    assert error < mp.mpf("2.0")

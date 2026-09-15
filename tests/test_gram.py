"""Tests for Gram points and Gram's law."""

from __future__ import annotations

import mpmath as mp
import pytest

from critical_line import gram
from critical_line.constants import FIRST_GRAM_LAW_FAILURE
from critical_line.hardy_z import riemann_siegel_theta


@pytest.mark.parametrize("n", [0, 1, 2, 5, 20, 40])
def test_gram_points_satisfy_defining_equation(n):
    """
    theta(g_n) = n * pi, to within the working precision.

    The comparison must happen *inside* a matching precision context: mp.pi
    evaluated at the ambient 15 digits would limit the residual to ~1e-16 and
    say nothing about the accuracy of the root.
    """
    g = gram.gram_point(n, dps=30)
    with mp.workdps(30):
        assert abs(riemann_siegel_theta(g, dps=30) - n * mp.pi) < mp.mpf("1e-20")


def test_first_gram_point_known_value():
    """g_0 = 17.8455... is the standard tabulated value."""
    assert abs(gram.gram_point(0, dps=25) - mp.mpf("17.8455995404")) < mp.mpf("1e-9")


def test_gram_points_are_increasing():
    points = [gram.gram_point(n, dps=25) for n in range(0, 12)]
    assert all(b > a for a, b in zip(points[:-1], points[1:]))


def test_gram_interval_report_shape():
    rows = gram.gram_interval_report(0, 8, dps=25)
    assert len(rows) == 8
    for row in rows:
        assert set(row) == {
            "n", "g_n", "z_at_g_n", "expected_sign",
            "observed_sign", "obeys_gram_law",
        }
        assert row["expected_sign"] in (1, -1)
        assert row["observed_sign"] in (1, -1, 0)


def test_gram_law_holds_for_small_indices():
    """
    No failures below n = 126. An empty result here is expected and is not
    evidence for anything -- Gram's law is known to fail infinitely often.
    """
    assert gram.gram_law_failures(1, 40, dps=25) == []


def test_gram_law_failure_constant_is_documented():
    assert FIRST_GRAM_LAW_FAILURE == 126


def test_report_and_failures_are_consistent():
    rows = gram.gram_interval_report(1, 25, dps=25)
    from_report = [r["n"] for r in rows if not r["obeys_gram_law"]]
    assert from_report == gram.gram_law_failures(1, 25, dps=25)


def test_z_sign_alternates_at_gram_points_in_good_range():
    rows = gram.gram_interval_report(1, 20, dps=25)
    signs = [r["observed_sign"] for r in rows]
    assert all(a * b < 0 for a, b in zip(signs[:-1], signs[1:]))

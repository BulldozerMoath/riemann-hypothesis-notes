"""Tests for zeta, the functional equation, and xi."""

from __future__ import annotations

import mpmath as mp
import pytest

from critical_line import zeta_tools as zt


@pytest.mark.parametrize("s", [
    mp.mpc("0.3", "1.7"),
    mp.mpc("-0.5", "4.0"),
    mp.mpc("0.75", "30.0"),
])
def test_functional_equation(s):
    assert zt.functional_equation_residual(s, dps=25) < mp.mpf("1e-20")


def test_functional_equation_sides_agree():
    with mp.workdps(30):
        lhs, rhs = zt.functional_equation_sides(mp.mpc("0.3", "1.7"), dps=30)
        assert abs(lhs - rhs) < mp.mpf("1e-22")


def test_special_values():
    with zt.working_precision(25):
        assert abs(mp.zeta(0) + mp.mpf("0.5")) < mp.mpf("1e-20")
        assert abs(mp.zeta(2) - mp.pi ** 2 / 6) < mp.mpf("1e-20")
        assert abs(mp.zeta(-1) + mp.mpf(1) / 12) < mp.mpf("1e-20")


@pytest.mark.parametrize("k", [-2, -4, -6, -8])
def test_trivial_zeros(k):
    assert abs(zt.zeta(k)) < mp.mpf("1e-20")


@pytest.mark.parametrize("re_part,im_part", [("0.3", "1.7"), ("0.9", "5.0")])
def test_xi_functional_equation(re_part, im_part):
    """
    xi(s) = xi(1-s).

    Both s and its reflection are constructed inside the precision context:
    forming 1-s at the ambient 15 digits would make the two arguments
    inexact reflections of each other and cap the agreement at ~1e-16.
    """
    with mp.workdps(30):
        s = mp.mpc(re_part, im_part)
        assert abs(zt.xi(s, dps=30) - zt.xi(1 - s, dps=30)) < mp.mpf("1e-22")


@pytest.mark.parametrize("t", [5, "14.134725141734693", "21.022039638771555"])
def test_xi_real_on_critical_line(t):
    with mp.workdps(30):
        value = zt.xi(mp.mpc(mp.mpf("0.5"), mp.mpmathify(t)), dps=30)
        assert abs(value.imag) < mp.mpf("1e-22")


def test_log_derivative_matches_finite_difference():
    """zeta'/zeta against a numerical derivative of log zeta."""
    with mp.workdps(30):
        s = mp.mpc("2.0", "3.0")
        analytic = zt.log_derivative(s)
        numeric = mp.diff(lambda z: mp.log(mp.zeta(z)), s)
        assert abs(analytic - numeric) < mp.mpf("1e-15")


def test_working_precision_restores_on_exit():
    before = mp.mp.dps
    with zt.working_precision(123):
        assert mp.mp.dps == 123
    assert mp.mp.dps == before


def test_working_precision_restores_on_exception():
    before = mp.mp.dps
    with pytest.raises(RuntimeError):
        with zt.working_precision(77):
            raise RuntimeError("boom")
    assert mp.mp.dps == before

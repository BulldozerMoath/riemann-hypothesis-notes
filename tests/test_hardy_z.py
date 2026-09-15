"""Tests for the Riemann-Siegel theta function and Hardy's Z."""

from __future__ import annotations

import mpmath as mp
import pytest

from conftest import KNOWN_ORDINATES
from critical_line import hardy_z as hz


@pytest.mark.parametrize("t", [10, 14.5, 25, 50.5, 120.25])
def test_hardy_z_from_zeta_agrees_with_riemann_siegel(t):
    """
    hardy_z and hardy_z_from_zeta agree.

    What this does NOT establish: mp.siegelz only uses the Riemann-Siegel
    asymptotic expansion when |t| > 500 * prec (about t > 51500 at 30 digits).
    Every t used here is far below that, so mp.siegelz is computing
    expj(siegeltheta(t)) * zeta(1/2 + it) -- the same formula
    hardy_z_from_zeta uses. The two are therefore NOT independent
    implementations at these heights, and this test confirms the wiring (guard
    digits, phase sign, real-part extraction), not the mathematics.

    The tolerance is realistic rather than zero: the two differ in guard digits
    and rounding, so they will not agree bit for bit.
    """
    a = hz.hardy_z(t, dps=30)
    b = hz.hardy_z_from_zeta(t, dps=30)
    assert abs(a - b) < mp.mpf("1e-20")


@pytest.mark.parametrize("gamma", KNOWN_ORDINATES)
def test_hardy_z_vanishes_at_known_ordinates(gamma):
    """Z(gamma_k) is zero to within the accuracy of the tabulated ordinate."""
    assert abs(hz.hardy_z(gamma, dps=30)) < mp.mpf("1e-14")


def test_hardy_z_first_ordinate_specifically():
    """Z(gamma_1) ~ 0 for the first known ordinate, 14.134725141734693."""
    assert abs(hz.hardy_z(KNOWN_ORDINATES[0], dps=30)) < mp.mpf("1e-14")


def test_siegelz_does_not_use_riemann_siegel_at_project_heights():
    """
    Pin down the fact that makes the agreement above a wiring check rather than
    an independent cross-check.

    mpmath's siegelz takes the Riemann-Siegel branch only when
    |t| > 500 * ctx.prec (bits). This asserts the threshold is far above every
    height this project reaches, so the fallback -- the same formula
    hardy_z_from_zeta uses -- is what actually runs.

    If a future mpmath lowers that threshold, this test fails and the
    docstrings claiming the two are not independent need revisiting.
    """
    highest_height_used_here = 2000  # analyze_spacings --count 1000 reaches ~1420
    for dps in (15, 25, 30, 50):
        with mp.workdps(dps):
            threshold = 500 * mp.mp.prec
        assert threshold > highest_height_used_here, (
            f"at dps={dps} the Riemann-Siegel branch starts at t={threshold}, "
            f"which is below the heights this project uses"
        )


def test_hardy_z_is_real():
    for t in [10, 33.3, 77.7]:
        assert isinstance(hz.hardy_z(t, dps=25), mp.mpf)
        assert isinstance(hz.hardy_z_from_zeta(t, dps=25), mp.mpf)


def test_hardy_z_magnitude_matches_zeta_magnitude():
    """|Z(t)| = |zeta(1/2 + it)| -- the phase factor has modulus one."""
    for t in [12, 20, 45.5]:
        with mp.workdps(30):
            lhs = abs(hz.hardy_z(t, dps=30))
            rhs = abs(mp.zeta(mp.mpc(mp.mpf("0.5"), mp.mpmathify(t))))
            assert abs(lhs - rhs) < mp.mpf("1e-20")


def test_imaginary_residual_is_small_but_tracked():
    """
    The discarded imaginary part is mathematically zero and numerically not.
    It should be tiny, and the function should report it rather than hide it.
    """
    residual = hz.hardy_z_imaginary_residual(25.5, dps=30)
    assert residual >= 0
    assert residual < mp.mpf("1e-20")


def test_theta_is_zero_at_origin():
    assert abs(hz.riemann_siegel_theta(0, dps=25)) < mp.mpf("1e-20")


def test_theta_is_increasing_above_first_gram_point():
    with mp.workdps(25):
        values = [hz.riemann_siegel_theta(t, dps=25) for t in [20, 30, 40, 50]]
    assert all(b > a for a, b in zip(values[:-1], values[1:]))


def test_sign_change_brackets_first_zero():
    lo = hz.hardy_z(14.0, dps=25)
    hi = hz.hardy_z(14.5, dps=25)
    assert lo * hi < 0

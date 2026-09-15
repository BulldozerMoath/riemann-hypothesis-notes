"""
Tests for zero counting, including the circularity result that motivated
this module and the residual-tolerance guard.
"""

from __future__ import annotations

import mpmath as mp
import pytest

from critical_line import counting as ct


# --------------------------------------------------------------------------
# the correction that motivated this module
# --------------------------------------------------------------------------

@pytest.mark.parametrize("T", [50, 100, 200])
def test_nzeros_is_not_independent_of_the_formula(T):
    """
    mp.nzeros(T) - (theta(T)/pi + 1) IS S(T), identically, to working
    precision. So comparing nzeros against the Riemann-von Mangoldt formula is
    not a cross-check of two computations -- it is one identity rearranged.

    This is pinned down as a test because an earlier version of this project
    reported that difference as though it were evidence about zero locations.
    """
    with mp.workdps(30):
        difference = mp.mpf(mp.nzeros(T)) - ct.rvm_main_term(T, dps=30)
        assert abs(difference - ct.S(T, dps=30)) < mp.mpf("1e-15")


def test_main_term_variants_agree_asymptotically():
    """theta(T)/pi + 1 and its Stirling expansion agree to O(1/T)."""
    for T in [100, 500, 1000]:
        a = ct.rvm_main_term(T, dps=25)
        b = ct.rvm_main_term_stirling(T, dps=25)
        assert abs(a - b) < mp.mpf("0.01")


# --------------------------------------------------------------------------
# the independent count
# --------------------------------------------------------------------------

@pytest.mark.parametrize("T,expected", [(20, 1), (30, 3)])
def test_argument_principle_count(T, expected):
    n, residual = ct.count_zeros_argument_principle_int(
        T, maxdegree=6, dps=25, residual_tolerance=mp.mpf("1e-6")
    )
    assert n == expected
    assert residual < mp.mpf("1e-6")


def test_argument_principle_below_first_zero():
    n, _ = ct.count_zeros_argument_principle_int(
        14, maxdegree=6, dps=25, residual_tolerance=mp.mpf("1e-6")
    )
    assert n == 0


def test_argument_principle_raises_on_unrealistic_tolerance():
    """
    An impossibly tight tolerance must produce a clear exception rather than a
    silently rounded integer. Refusing to round an unconverged result is the
    whole point of the guard.
    """
    with pytest.raises(ct.ResidualToleranceError) as excinfo:
        ct.count_zeros_argument_principle_int(
            20, maxdegree=4, dps=20, residual_tolerance=mp.mpf("1e-300")
        )
    error = excinfo.value
    assert error.residual > error.tolerance
    assert "did not converge" in str(error)
    assert "Refusing to round" in str(error)


def test_residual_tolerance_error_carries_diagnostics():
    try:
        ct.count_zeros_argument_principle_int(
            20, maxdegree=4, dps=20, residual_tolerance=mp.mpf("1e-300")
        )
    except ct.ResidualToleranceError as error:
        assert error.value is not None
        assert error.residual >= 0
        assert error.tolerance == mp.mpf("1e-300")
    else:
        pytest.fail("expected ResidualToleranceError")


# --------------------------------------------------------------------------
# on-line counting
# --------------------------------------------------------------------------

@pytest.mark.parametrize("T,expected", [(30, 3), (50, 10), (100, 29)])
def test_sign_change_count(T, expected):
    assert ct.count_sign_changes_hardy_z(T, dps=25) == expected


def test_refine_sign_change_finds_first_zero():
    t = ct.refine_sign_change(14.0, 14.5, dps=25)
    assert abs(t - mp.mpf("14.134725141734693")) < mp.mpf("1e-12")


# --------------------------------------------------------------------------
# the two together
# --------------------------------------------------------------------------

@pytest.mark.parametrize("T", [30, 50])
def test_experimental_check_agrees(T):
    record = ct.experimentally_check_up_to(
        T, dps=25, contour_dps=30, maxdegree=8,
        residual_tolerance=mp.mpf("1e-8"),
    )
    assert record["agree"], (
        f"strip count {record['strip_count']} != on-line count "
        f"{record['on_line_count']}; far more likely a bug or a quadrature "
        f"failure than a counterexample"
    )
    assert record["strip_count"] == record["on_line_count"]


def test_experimental_check_record_shape():
    record = ct.experimentally_check_up_to(
        20, dps=20, contour_dps=25, maxdegree=6,
        residual_tolerance=mp.mpf("1e-6"),
    )
    assert set(record) == {
        "T", "strip_count", "strip_residual", "on_line_count",
        "agree", "main_term", "S", "note",
    }


def test_experimental_check_carries_the_standing_note():
    """The comparison must always ship with its own caveat attached."""
    record = ct.experimentally_check_up_to(
        20, dps=20, contour_dps=25, maxdegree=6,
        residual_tolerance=mp.mpf("1e-6"),
    )
    assert "experimental evidence rather than a certified result" in record["note"]


def test_no_function_named_certify():
    """
    mpmath is not interval-certified, so 'certify' must not appear as an
    action name anywhere in the package's public API.
    """
    import critical_line

    for module_name in critical_line.__all__:
        module = getattr(critical_line, module_name, None)
        if module is None or not hasattr(module, "__dict__"):
            continue
        for attr in vars(module):
            if callable(getattr(module, attr, None)):
                assert not attr.startswith("certify"), (
                    f"{module_name}.{attr} uses 'certify' as an action name"
                )

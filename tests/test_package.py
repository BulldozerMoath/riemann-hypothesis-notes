"""
Tests that the installed package namespace works as documented.

These exercise `from critical_line import ...` rather than any path hack, so
they fail if the package is not installed (``pip install -e .``).
"""

from __future__ import annotations

import pytest


def test_import_counting_through_package_namespace():
    from critical_line import counting

    assert hasattr(counting, "experimentally_check_up_to")


@pytest.mark.parametrize("name", [
    "constants", "counting", "explicit_formula", "gram",
    "hardy_z", "spacing", "zeros", "zeta_tools",
])
def test_submodules_importable(name):
    import importlib

    module = importlib.import_module(f"critical_line.{name}")
    assert module.__name__ == f"critical_line.{name}"


def test_plots_is_lazily_available():
    import critical_line

    assert critical_line.plots.__name__ == "critical_line.plots"


def test_unknown_attribute_raises():
    import critical_line

    with pytest.raises(AttributeError):
        critical_line.no_such_module


def test_renamed_function_is_present_and_old_name_is_gone():
    from critical_line import counting

    assert callable(counting.experimentally_check_up_to)
    assert not hasattr(counting, "certify_up_to")


def test_package_exposes_disclaimers():
    import critical_line

    assert "does not" in critical_line.DISCLAIMER_SHORT.lower() or \
           "not prove" in critical_line.DISCLAIMER_SHORT.lower()
    assert "experimental evidence" in critical_line.DISCLAIMER_LONG


def test_version_is_set():
    import critical_line

    assert critical_line.__version__

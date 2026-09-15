"""
Tests for figure generation.

Figures are written to pytest's tmp_path, never to the repository, so the test
suite leaves no artifacts behind.
"""

from __future__ import annotations

import pytest

PNG_MAGIC = b"\x89PNG\r\n\x1a\n"
MIN_BYTES = 2000


@pytest.fixture(scope="module")
def plots():
    return pytest.importorskip("critical_line.plots")


def _assert_real_png(path):
    assert path.exists(), f"{path} was not created"
    data = path.read_bytes()
    assert len(data) > MIN_BYTES, f"{path} is suspiciously small ({len(data)} bytes)"
    assert data.startswith(PNG_MAGIC), f"{path} is not a PNG"


def test_plot_hardy_z(plots, tmp_path):
    out = plots.plot_hardy_z(10, 40, tmp_path / "hardy_z.png", dps=20, samples=120)
    _assert_real_png(out)


def test_plot_zero_ordinates(plots, tmp_path, zeros_to_60):
    out = plots.plot_zero_ordinates(zeros_to_60, tmp_path / "ordinates.png")
    _assert_real_png(out)


def test_plot_zero_ordinates_rejects_empty(plots, tmp_path):
    with pytest.raises(ValueError, match="no zeros"):
        plots.plot_zero_ordinates([], tmp_path / "empty.png")


def test_plot_normalized_gap_histogram(plots, tmp_path, zeros_to_150):
    out = plots.plot_normalized_gap_histogram(zeros_to_150, tmp_path / "gaps.png")
    _assert_real_png(out)


def test_plot_normalized_gap_histogram_needs_three_zeros(plots, tmp_path):
    with pytest.raises(ValueError, match="at least 3"):
        plots.plot_normalized_gap_histogram(["14.1", "21.0"], tmp_path / "gaps.png")


def test_plot_counting_comparison(plots, tmp_path, zeros_to_60):
    out = plots.plot_counting_comparison(
        [10, 20, 30, 40, 50, 60], tmp_path / "counting.png",
        dps=20, zeros=zeros_to_60,
    )
    _assert_real_png(out)


def test_plot_counting_comparison_rejects_empty_heights(plots, tmp_path):
    with pytest.raises(ValueError, match="empty"):
        plots.plot_counting_comparison([], tmp_path / "counting.png", dps=20, zeros=[])


def test_plot_explicit_formula_error(plots, tmp_path, zeros_to_60):
    out = plots.plot_explicit_formula_error(100.0, zeros_to_60, tmp_path / "ef.png")
    _assert_real_png(out)


def test_plot_explicit_formula_rejects_empty(plots, tmp_path):
    with pytest.raises(ValueError, match="no zeros"):
        plots.plot_explicit_formula_error(100.0, [], tmp_path / "ef.png")


def test_plots_create_missing_parent_directories(plots, tmp_path, zeros_to_60):
    nested = tmp_path / "a" / "b" / "c" / "ordinates.png"
    _assert_real_png(plots.plot_zero_ordinates(zeros_to_60, nested))

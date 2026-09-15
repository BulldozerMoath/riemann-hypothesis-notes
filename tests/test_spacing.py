"""Tests for gap statistics."""

from __future__ import annotations

import mpmath as mp
import pytest

from critical_line import spacing


def test_zero_gaps_are_differences():
    gaps = spacing.zero_gaps(["1.0", "3.5", "4.0"])
    assert len(gaps) == 2
    assert abs(gaps[0] - mp.mpf("2.5")) < mp.mpf("1e-25")
    assert abs(gaps[1] - mp.mpf("0.5")) < mp.mpf("1e-25")


def test_gaps_empty_for_single_zero():
    assert spacing.zero_gaps(["14.13"]) == []
    assert spacing.normalized_gaps(["14.13"]) == []


def test_normalized_gaps_are_positive(zeros_to_150):
    """Ordinates increase, and the local mean spacing is positive."""
    gaps = spacing.normalized_gaps(zeros_to_150)
    assert len(gaps) == len(zeros_to_150) - 1
    assert all(g > 0 for g in gaps)


def test_raw_gaps_are_positive(zeros_to_150):
    assert all(g > 0 for g in spacing.zero_gaps(zeros_to_150))


def test_local_mean_spacing_decreases_with_height():
    """Zeros get denser higher up, so the mean spacing shrinks."""
    values = [spacing.local_mean_spacing(t) for t in [50, 500, 5000]]
    assert all(b < a for a, b in zip(values[:-1], values[1:]))


def test_normalized_mean_is_near_one(zeros_to_150):
    """
    Near 1 by construction -- the unfolding is designed to produce it. This
    checks the arithmetic, not the zeros.
    """
    summary = spacing.spacing_summary(zeros_to_150)
    assert mp.mpf("0.8") < summary["normalized_mean"] < mp.mpf("1.2")


def test_spacing_summary_shape(zeros_to_150):
    summary = spacing.spacing_summary(zeros_to_150)
    assert set(summary) == {
        "n_zeros", "n_gaps", "height_min", "height_max",
        "raw_gap_min", "raw_gap_max", "normalized_mean",
        "normalized_min", "normalized_max", "normalized_std",
    }
    assert summary["n_gaps"] == summary["n_zeros"] - 1
    assert summary["normalized_min"] <= summary["normalized_mean"]
    assert summary["normalized_mean"] <= summary["normalized_max"]


def test_spacing_summary_requires_two_zeros():
    with pytest.raises(ValueError, match="at least 2"):
        spacing.spacing_summary(["14.13"])


def test_lehmer_pairs_empty_at_low_heights(zeros_to_150):
    """Low-lying zeros are well separated; close pairs appear much higher."""
    assert spacing.lehmer_like_pairs(zeros_to_150) == []


def test_lehmer_pairs_detects_a_close_pair():
    """A synthetic close pair must be reported, with its record populated."""
    close = ["1000.0", "1000.01", "1002.0"]
    pairs = spacing.lehmer_like_pairs(close, threshold=mp.mpf("0.5"))
    assert len(pairs) == 1
    assert pairs[0]["index"] == 0
    assert pairs[0]["normalized_gap"] < mp.mpf("0.5")
    assert pairs[0]["gamma_lower"] == mp.mpf("1000.0")


def test_lehmer_threshold_is_respected(zeros_to_150):
    """A threshold above every observed gap catches everything."""
    gaps = spacing.normalized_gaps(zeros_to_150)
    generous = max(gaps) + 1
    assert len(spacing.lehmer_like_pairs(zeros_to_150, threshold=generous)) == len(gaps)

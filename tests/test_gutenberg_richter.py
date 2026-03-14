"""Tests for src/gutenberg_richter.py — b-value estimation and changepoints."""

import numpy as np
import pytest

from src.gutenberg_richter import (
    estimate_b_value,
    bootstrap_b_value,
    segment_cost,
    detect_changepoints,
    ks_test_magnitudes,
)


class TestEstimateBValue:
    def test_known_b_value(self):
        # Generate magnitudes from GR distribution with b=1.0, mc=2.0
        rng = np.random.default_rng(123)
        b_true = 1.0
        mc = 2.0
        # Magnitudes above mc follow exponential: M - mc ~ Exp(b * ln(10))
        mags = mc + rng.exponential(scale=1.0 / (b_true * np.log(10)), size=5000)

        b_est, b_std = estimate_b_value(mags, mc)
        assert abs(b_est - b_true) < 0.1, f"Expected ~1.0, got {b_est}"
        assert b_std > 0

    def test_returns_tuple(self):
        b, std = estimate_b_value([3.0, 3.5, 4.0, 4.5, 5.0], mc=2.5)
        assert isinstance(b, float)
        assert isinstance(std, float)

    def test_raises_on_too_few_events(self):
        with pytest.raises(ValueError, match="at least 2"):
            estimate_b_value([1.0], mc=2.0)

    def test_filters_below_mc(self):
        mags = [1.0, 1.5, 3.0, 3.5, 4.0]
        b, _ = estimate_b_value(mags, mc=2.5)
        # Should only use the 3 events >= 2.5
        assert b > 0


class TestBootstrapBValue:
    def test_ci_is_reasonable(self):
        rng = np.random.default_rng(456)
        b_true = 1.0
        mc = 2.0
        mags = mc + rng.exponential(scale=1.0 / (b_true * np.log(10)), size=2000)

        b_mean, b_low, b_high = bootstrap_b_value(mags, mc, n_bootstrap=500)
        # CI should be narrow and centered near the true value
        assert b_low < b_mean < b_high
        assert abs(b_mean - b_true) < 0.2
        assert (b_high - b_low) < 0.5


class TestSegmentCost:
    def test_zero_variance_returns_zero(self):
        assert segment_cost([5.0, 5.0, 5.0], 0, 3) == 0.0

    def test_positive_cost_for_variable_data(self):
        cost = segment_cost([1.0, 2.0, 3.0, 4.0, 5.0], 0, 5)
        assert cost > 0

    def test_short_segment_returns_zero(self):
        assert segment_cost([1.0], 0, 1) == 0.0


class TestDetectChangepoints:
    def test_detects_mean_shift(self):
        # Two segments with different means
        rng = np.random.default_rng(789)
        seg1 = rng.normal(0, 1, 100)
        seg2 = rng.normal(5, 1, 100)
        data = np.concatenate([seg1, seg2])

        cps = detect_changepoints(data, min_size=20)
        # Should find a changepoint near index 100
        assert len(cps) >= 1
        assert any(80 <= cp <= 120 for cp in cps)

    def test_no_changepoint_in_constant_data(self):
        data = np.ones(100)
        cps = detect_changepoints(data, min_size=20)
        assert len(cps) == 0

    def test_too_short_returns_empty(self):
        cps = detect_changepoints([1.0, 2.0, 3.0], min_size=20)
        assert cps == []


class TestKsTestMagnitudes:
    def test_same_distribution_high_pvalue(self):
        rng = np.random.default_rng(0)
        mags = rng.exponential(1.0, 500) + 2.5
        stat, pval = ks_test_magnitudes(mags[:250], mags[250:], mc=2.5)
        assert pval > 0.05

    def test_different_distribution_low_pvalue(self):
        rng = np.random.default_rng(0)
        m1 = rng.exponential(0.5, 500) + 2.5  # b ~ 0.87
        m2 = rng.exponential(2.0, 500) + 2.5  # b ~ 0.22
        stat, pval = ks_test_magnitudes(m1, m2, mc=2.5)
        assert pval < 0.05

    def test_raises_on_empty_sample(self):
        with pytest.raises(ValueError):
            ks_test_magnitudes([1.0], [1.0], mc=5.0)

"""Tests for src.gutenberg_richter module."""

import numpy as np
import pytest

from src.gutenberg_richter import estimate_b_value


def _synthetic_magnitudes(b_true, mc, n, seed=42):
    """Generate magnitudes from the inverse GR CDF: M = mc - log10(U) / b."""
    rng = np.random.default_rng(seed)
    u = rng.uniform(0, 1, n)
    return mc - np.log10(u) / b_true


class TestEstimateBValue:
    def test_recovery_known_b(self):
        """Recover a known b-value from synthetic data within tolerance."""
        b_true = 1.0
        mc = 2.5
        mags = _synthetic_magnitudes(b_true, mc, n=5000)
        b_est, _ = estimate_b_value(mags, mc)
        assert abs(b_est - b_true) < 0.15, f"Expected ~{b_true}, got {b_est}"

    def test_recovery_high_b(self):
        """Recover a higher b-value (1.5)."""
        b_true = 1.5
        mc = 2.0
        mags = _synthetic_magnitudes(b_true, mc, n=5000)
        b_est, _ = estimate_b_value(mags, mc)
        assert abs(b_est - b_true) < 0.25

    def test_std_error_formula(self):
        """Verify that std error is reasonable."""
        mc = 2.5
        mags = _synthetic_magnitudes(1.0, mc, n=1000)
        b, se = estimate_b_value(mags, mc)
        # Std error should be positive and much smaller than b
        assert se > 0
        assert se < b * 0.5

    def test_too_few_events_raises(self):
        """Fewer than 2 events above mc should raise ValueError."""
        with pytest.raises(ValueError):
            estimate_b_value([1.0, 1.5], mc=5.0)

    def test_single_event_raises(self):
        """A single event above mc should raise ValueError."""
        with pytest.raises(ValueError):
            estimate_b_value([3.0], mc=2.0)

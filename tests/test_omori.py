"""Tests for src/omori.py — Omori Law fitting and recovery analysis."""

import numpy as np
import pytest

from src.omori import (
    _omori_rate,
    compute_recovery_fraction,
    fit_recovery_exponential,
)


class TestOmoriRate:
    def test_decay_with_known_params(self):
        # R(t) = K * (t + c)^(-p)
        K, c, p = 100.0, 1.0, 1.0
        assert _omori_rate(0, K, c, p) == pytest.approx(100.0)
        assert _omori_rate(9, K, c, p) == pytest.approx(10.0)

    def test_rate_decreases_with_time(self):
        K, c, p = 50.0, 0.5, 1.2
        rates = [_omori_rate(t, K, c, p) for t in [1, 10, 100, 1000]]
        assert all(rates[i] > rates[i + 1] for i in range(len(rates) - 1))


class TestComputeRecoveryFraction:
    def test_full_recovery(self):
        # b goes from perturbed back to baseline
        b_baseline = 1.0
        b_timeseries = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
        fractions = compute_recovery_fraction(b_timeseries, b_baseline)
        assert fractions[-1] == pytest.approx(1.0)
        assert fractions[0] == pytest.approx(0.0)

    def test_no_perturbation(self):
        fractions = compute_recovery_fraction([1.0, 1.0, 1.0], 1.0)
        assert all(f == pytest.approx(1.0) for f in fractions)

    def test_monotonic_recovery(self):
        b_timeseries = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
        fractions = compute_recovery_fraction(b_timeseries, 1.0)
        assert all(fractions[i] <= fractions[i + 1]
                   for i in range(len(fractions) - 1))


class TestFitRecoveryExponential:
    def test_recovers_known_tau(self):
        tau_true = 100.0
        t = np.linspace(1, 500, 50)
        f = 1.0 - np.exp(-t / tau_true)
        # Add tiny noise
        rng = np.random.default_rng(42)
        f += rng.normal(0, 0.01, len(f))

        result = fit_recovery_exponential(f, t)
        assert result is not None
        assert abs(result["tau_b"] - tau_true) < 20
        assert result["r_squared"] > 0.9

    def test_returns_none_on_insufficient_data(self):
        result = fit_recovery_exponential([0.5, 0.6], [1.0, 2.0])
        assert result is None

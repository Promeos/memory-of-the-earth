"""Tests for src.interevent module."""

import numpy as np
import pytest
from scipy import stats

from src.interevent import fit_distributions, classify_regime


class TestFitDistributions:
    def test_weibull_data_detected(self):
        """Data drawn from Weibull should be best-fit by Weibull (or close)."""
        rng = np.random.default_rng(42)
        data = stats.weibull_min.rvs(0.7, scale=100, size=2000,
                                      random_state=rng)
        result = fit_distributions(data)
        assert result["best_aic"] in ("weibull", "gamma", "lognormal")

    def test_exponential_data(self):
        """Exponential data should be best-fit by exponential or Weibull k~1."""
        rng = np.random.default_rng(42)
        data = stats.expon.rvs(scale=50, size=2000, random_state=rng)
        result = fit_distributions(data)
        best = result["best_aic"]
        if best == "weibull":
            k = result["weibull"]["params"][0]
            assert 0.8 < k < 1.2
        else:
            assert best in ("exponential", "gamma")

    def test_returns_all_distributions(self):
        """Result should contain entries for all four distributions."""
        data = np.random.default_rng(42).exponential(10, 500)
        result = fit_distributions(data)
        for name in ("exponential", "gamma", "lognormal", "weibull"):
            assert name in result
            assert "aic" in result[name]


class TestClassifyRegime:
    def test_clustered_low_k(self):
        """Weibull with k < 0.9 should classify as clustered."""
        rng = np.random.default_rng(42)
        data = stats.weibull_min.rvs(0.5, scale=100, size=5000,
                                      random_state=rng)
        result = fit_distributions(data)
        regime = classify_regime(result)
        assert regime == "clustered"

    def test_quasi_periodic_high_k(self):
        """Weibull with k > 1.1 should classify as quasi_periodic."""
        rng = np.random.default_rng(42)
        data = stats.weibull_min.rvs(2.0, scale=100, size=5000,
                                      random_state=rng)
        result = fit_distributions(data)
        regime = classify_regime(result)
        assert regime == "quasi_periodic"

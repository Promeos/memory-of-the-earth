"""Tests for src/interevent.py — inter-event times and distribution fitting."""

import numpy as np
import pandas as pd
import pytest

from src.interevent import (
    compute_interevent_times,
    fit_distributions,
    classify_regime,
)


class TestComputeIntereventTimes:
    def test_basic_computation(self):
        times = pd.to_datetime([
            "2020-01-01 00:00:00",
            "2020-01-01 00:01:00",
            "2020-01-01 00:03:00",
        ])
        iet = compute_interevent_times(times)
        assert len(iet) == 2
        assert iet[0] == pytest.approx(60.0)
        assert iet[1] == pytest.approx(120.0)

    def test_returns_seconds(self):
        times = pd.to_datetime([
            "2020-01-01 00:00:00",
            "2020-01-01 01:00:00",
        ])
        iet = compute_interevent_times(times)
        assert iet[0] == pytest.approx(3600.0)

    def test_single_event_returns_empty(self):
        times = pd.to_datetime(["2020-01-01"])
        iet = compute_interevent_times(times)
        assert len(iet) == 0


class TestFitDistributions:
    def test_returns_all_distributions(self):
        rng = np.random.default_rng(42)
        iet = rng.exponential(scale=100, size=500)
        result = fit_distributions(iet)

        for name in ["exponential", "gamma", "lognormal", "weibull"]:
            assert name in result
            assert "aic" in result[name]
            assert "bic" in result[name]
            assert "loglik" in result[name]

        assert "best_aic" in result
        assert "best_bic" in result

    def test_exponential_data_prefers_exponential(self):
        rng = np.random.default_rng(42)
        iet = rng.exponential(scale=100, size=2000)
        result = fit_distributions(iet)
        # Exponential or Weibull with k~1 should win
        assert result["best_aic"] in ("exponential", "weibull")


class TestClassifyRegime:
    def test_poisson_classification(self):
        # Mock a fit result where exponential clearly wins
        fit_result = {
            "exponential": {"aic": 100, "params": ()},
            "gamma": {"aic": 200, "params": (1.0, 0, 1)},
            "lognormal": {"aic": 200, "params": (1.0, 0, 1)},
            "weibull": {"aic": 200, "params": (1.0, 0, 1)},
        }
        assert classify_regime(fit_result) == "poisson"

    def test_clustered_classification(self):
        fit_result = {
            "exponential": {"aic": 200, "params": ()},
            "gamma": {"aic": 200, "params": (0.5, 0, 1)},
            "lognormal": {"aic": 200, "params": (1.0, 0, 1)},
            "weibull": {"aic": 100, "params": (0.5, 0, 1)},  # k < 0.9
        }
        assert classify_regime(fit_result) == "clustered"

    def test_quasi_periodic_classification(self):
        fit_result = {
            "exponential": {"aic": 200, "params": ()},
            "gamma": {"aic": 200, "params": (1.0, 0, 1)},
            "lognormal": {"aic": 200, "params": (1.0, 0, 1)},
            "weibull": {"aic": 100, "params": (3.0, 0, 1)},  # k > 1.1
        }
        assert classify_regime(fit_result) == "quasi_periodic"

    def test_ambiguous_when_close_aic(self):
        fit_result = {
            "exponential": {"aic": 100, "params": ()},
            "gamma": {"aic": 105, "params": (1.0, 0, 1)},
            "lognormal": {"aic": 110, "params": (1.0, 0, 1)},
            "weibull": {"aic": 108, "params": (1.0, 0, 1)},
        }
        assert classify_regime(fit_result) == "ambiguous"

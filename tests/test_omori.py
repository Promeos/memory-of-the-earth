"""Tests for src.omori module."""

import numpy as np
import pandas as pd
import pytest

from src.omori import fit_omori, compute_recovery_fraction


class TestFitOmori:
    def test_synthetic_omori_decay(self):
        """Fit synthetic Omori-decay data and recover approximate p."""
        rng = np.random.default_rng(42)
        K_true, c_true, p_true = 50.0, 1.0, 1.2

        mainshock = pd.Timestamp("2020-01-01", tz="UTC")
        events = []
        for day in range(1, 500):
            rate = K_true * (day + c_true) ** (-p_true)
            n_events = rng.poisson(rate)
            for _ in range(n_events):
                offset = rng.uniform(0, 1)
                events.append(mainshock + pd.Timedelta(days=day + offset))

        if len(events) < 20:
            pytest.skip("Too few synthetic events generated")

        result = fit_omori(pd.Series(events), mainshock)
        assert result is not None
        assert 0.5 < result["p"] < 3.0
        assert result["r_squared"] > 0.3

    def test_too_few_events_returns_none(self):
        """Fewer than 10 events should return None."""
        mainshock = pd.Timestamp("2020-01-01", tz="UTC")
        events = pd.Series([mainshock + pd.Timedelta(days=i) for i in range(1, 5)])
        result = fit_omori(events, mainshock)
        assert result is None


class TestComputeRecoveryFraction:
    def test_full_recovery(self):
        """If b returns to baseline, recovery fraction should approach 1."""
        b_series = np.array([0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
        fractions = compute_recovery_fraction(b_series, b_baseline=1.0)
        assert abs(fractions[-1] - 1.0) < 0.01

    def test_no_perturbation(self):
        """If b starts at baseline, fractions should all be 1."""
        b_series = np.array([1.0, 1.0, 1.0])
        fractions = compute_recovery_fraction(b_series, b_baseline=1.0)
        np.testing.assert_allclose(fractions, 1.0)

    def test_no_recovery(self):
        """If b stays at initial perturbed value, fraction stays at 0."""
        b_series = np.array([0.5, 0.5, 0.5])
        fractions = compute_recovery_fraction(b_series, b_baseline=1.0)
        np.testing.assert_allclose(fractions, 0.0)

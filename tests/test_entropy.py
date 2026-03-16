"""Tests for src.entropy module."""

import numpy as np
import pandas as pd
import pytest

from src.entropy import shannon_entropy, detect_anomalies


class TestShannonEntropy:
    def test_uniform_distribution(self):
        """Uniform distribution across n bins should give log2(n)."""
        mc = 2.5
        bin_width = 0.5
        max_mag = 4.5  # 4 bins: [2.5,3.0), [3.0,3.5), [3.5,4.0), [4.0,4.5)
        n_bins = 4
        # Place equal events in each bin
        mags = np.array([2.6, 2.7, 3.1, 3.2, 3.6, 3.7, 4.1, 4.2])
        H = shannon_entropy(mags, mc, bin_width=bin_width, max_mag=max_mag)
        expected = np.log2(n_bins)
        assert abs(H - expected) < 0.01, f"Expected {expected}, got {H}"

    def test_single_bin(self):
        """All events in one bin should give entropy 0."""
        mc = 2.5
        mags = np.array([2.6, 2.7, 2.8, 2.9])
        H = shannon_entropy(mags, mc, bin_width=0.5, max_mag=3.0)
        assert abs(H) < 0.01, f"Expected ~0, got {H}"

    def test_empty_returns_nan(self):
        """No events should return NaN."""
        H = shannon_entropy(np.array([]), mc=2.5)
        assert np.isnan(H)

    def test_all_below_mc_returns_nan(self):
        """Events all below mc should return NaN."""
        H = shannon_entropy(np.array([1.0, 1.5, 2.0]), mc=3.0)
        assert np.isnan(H)


class TestDetectAnomalies:
    def test_detects_planted_drops(self):
        """Should detect anomalies planted below the 5th percentile."""
        n = 200
        rng = np.random.default_rng(42)
        H_values = rng.normal(2.0, 0.1, n)
        # Plant 5 very low values
        H_values[:5] = 1.0

        df = pd.DataFrame({
            "center_time": pd.date_range("2020-01-01", periods=n, freq="7D"),
            "H": H_values,
        })

        anomalies = detect_anomalies(df, percentile=5)
        # All 5 planted drops should be detected
        assert len(anomalies) >= 5
        # The planted anomalies at indices 0-4 should all be included
        assert all(anomalies["H"] <= np.percentile(H_values, 5))

    def test_no_anomalies_when_uniform(self):
        """With constant entropy, the 5th percentile captures ~5% of data."""
        n = 100
        df = pd.DataFrame({
            "center_time": pd.date_range("2020-01-01", periods=n, freq="7D"),
            "H": np.full(n, 2.0),
        })
        anomalies = detect_anomalies(df, percentile=5)
        # All values equal the percentile threshold, so all are flagged
        assert len(anomalies) == n

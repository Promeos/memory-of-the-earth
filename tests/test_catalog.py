"""Tests for src/catalog.py — cleaning, deduplication, and Mc estimation."""

import numpy as np
import pandas as pd
import pytest

from src.catalog import clean_catalog, deduplicate, estimate_mc


class TestCleanCatalog:
    def test_filters_non_earthquake_types(self):
        df = pd.DataFrame({
            "time": ["2020-01-01T00:00:00", "2020-01-02T00:00:00",
                     "2020-01-03T00:00:00"],
            "mag": [3.0, 2.5, 4.0],
            "type": ["earthquake", "quarry blast", "earthquake"],
        })
        result = clean_catalog(df)
        assert len(result) == 2
        assert all(result["type"] == "earthquake")

    def test_drops_nan_magnitudes(self):
        df = pd.DataFrame({
            "time": ["2020-01-01T00:00:00", "2020-01-02T00:00:00"],
            "mag": [3.0, np.nan],
            "type": ["earthquake", "earthquake"],
        })
        result = clean_catalog(df)
        assert len(result) == 1

    def test_parses_time_to_utc(self):
        df = pd.DataFrame({
            "time": ["2020-06-15T12:00:00"],
            "mag": [5.0],
        })
        result = clean_catalog(df)
        assert result["time"].dt.tz is not None

    def test_sorts_by_time(self):
        df = pd.DataFrame({
            "time": ["2020-03-01", "2020-01-01", "2020-02-01"],
            "mag": [3.0, 2.0, 4.0],
        })
        result = clean_catalog(df)
        assert result["time"].is_monotonic_increasing

    def test_works_without_type_column(self):
        df = pd.DataFrame({
            "time": ["2020-01-01T00:00:00"],
            "mag": [3.0],
        })
        result = clean_catalog(df)
        assert len(result) == 1


class TestDeduplicate:
    def test_removes_duplicate_events(self):
        base_time = pd.Timestamp("2020-01-01 00:00:00", tz="UTC")
        df = pd.DataFrame({
            "time": [base_time, base_time + pd.Timedelta(seconds=1)],
            "latitude": [35.0, 35.01],
            "longitude": [-97.0, -97.01],
            "mag": [3.0, 3.1],
            "nst": [10, 20],
        })
        result = deduplicate(df)
        assert len(result) == 1
        assert result.iloc[0]["nst"] == 20  # keeps higher station count

    def test_keeps_distinct_events(self):
        base_time = pd.Timestamp("2020-01-01 00:00:00", tz="UTC")
        df = pd.DataFrame({
            "time": [base_time, base_time + pd.Timedelta(hours=1)],
            "latitude": [35.0, 40.0],
            "longitude": [-97.0, -120.0],
            "mag": [3.0, 4.0],
            "nst": [10, 15],
        })
        result = deduplicate(df)
        assert len(result) == 2


class TestEstimateMc:
    def test_maxc_with_known_distribution(self):
        # Create magnitudes peaked at 2.5
        rng = np.random.default_rng(0)
        mags = np.concatenate([
            rng.uniform(2.0, 2.5, 200),  # below completeness
            rng.exponential(0.5, 800) + 2.5,  # above completeness
        ])
        mc = estimate_mc(mags)
        # Mc should be near 2.5 + 0.2 = 2.7 (mode bin + correction)
        assert 2.3 <= mc <= 3.2

    def test_returns_float(self):
        mc = estimate_mc([2.0, 2.5, 3.0, 3.5, 4.0])
        assert isinstance(mc, float)

    def test_empty_array_returns_zero(self):
        mc = estimate_mc([])
        assert mc == 0.0

    def test_rejects_unknown_method(self):
        with pytest.raises(ValueError, match="Unsupported method"):
            estimate_mc([3.0, 4.0], method="unknown")

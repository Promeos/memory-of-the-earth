"""Tests for src.catalog module."""

import numpy as np
import pandas as pd
import pytest

from src.catalog import estimate_mc, deduplicate, clean_catalog


class TestEstimateMc:
    def test_known_mode(self):
        """Mc should be the mode of the FMD + 0.2."""
        rng = np.random.default_rng(42)
        mags = np.concatenate([
            rng.uniform(1.5, 2.05, 500),
            rng.uniform(2.0, 5.0, 200),
        ])
        mc = estimate_mc(mags)
        assert 1.5 <= mc <= 3.0

    def test_empty_returns_nan(self):
        """Empty magnitude array returns NaN."""
        assert np.isnan(estimate_mc(np.array([])))

    def test_unsupported_method_raises(self):
        with pytest.raises(ValueError, match="Unsupported method"):
            estimate_mc([1.0, 2.0, 3.0], method="gft")


class TestDeduplicate:
    def test_removes_duplicates(self):
        """Events within 2s and 0.05 deg should be deduplicated."""
        df = pd.DataFrame({
            "time": pd.to_datetime([
                "2020-01-01 00:00:00",
                "2020-01-01 00:00:01",
                "2020-01-01 12:00:00",
            ], utc=True),
            "latitude": [35.0, 35.01, 40.0],
            "longitude": [-100.0, -100.01, -80.0],
            "mag": [3.0, 3.5, 4.0],
            "nst": [10, 20, 15],
        })
        result = deduplicate(df)
        assert len(result) == 2
        assert result.iloc[0]["mag"] == 3.5

    def test_no_duplicates_unchanged(self):
        """Catalog with no duplicates should be unchanged."""
        df = pd.DataFrame({
            "time": pd.to_datetime([
                "2020-01-01", "2020-06-01", "2020-12-01",
            ], utc=True),
            "latitude": [10.0, 20.0, 30.0],
            "longitude": [100.0, 110.0, 120.0],
            "mag": [3.0, 4.0, 5.0],
        })
        result = deduplicate(df)
        assert len(result) == 3


class TestCleanCatalog:
    def test_filters_non_earthquakes(self):
        """Non-earthquake types should be removed."""
        df = pd.DataFrame({
            "time": ["2020-01-01", "2020-01-02", "2020-01-03"],
            "mag": [3.0, 2.5, 4.0],
            "type": ["earthquake", "quarry blast", "earthquake"],
        })
        result = clean_catalog(df)
        assert len(result) == 2

    def test_drops_nan_mag(self):
        """Events with NaN magnitude should be dropped."""
        df = pd.DataFrame({
            "time": ["2020-01-01", "2020-01-02"],
            "mag": [3.0, np.nan],
            "type": ["earthquake", "earthquake"],
        })
        result = clean_catalog(df)
        assert len(result) == 1

    def test_sorts_by_time(self):
        """Output should be sorted by time."""
        df = pd.DataFrame({
            "time": ["2020-06-01", "2020-01-01", "2020-03-01"],
            "mag": [3.0, 2.0, 4.0],
        })
        result = clean_catalog(df)
        assert result["time"].is_monotonic_increasing

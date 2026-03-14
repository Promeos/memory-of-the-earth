"""Shared fixtures for the test suite."""

import numpy as np
import pandas as pd
import pytest


@pytest.fixture
def rng():
    """Seeded random generator for reproducible tests."""
    return np.random.default_rng(42)


@pytest.fixture
def synthetic_catalog(rng):
    """Small synthetic earthquake catalog for testing."""
    n = 500
    base = pd.Timestamp("2020-01-01", tz="UTC")
    times = base + pd.to_timedelta(rng.uniform(0, 365 * 3, n), unit="D")
    times = times.sort_values()

    return pd.DataFrame({
        "time": times,
        "latitude": rng.uniform(-90, 90, n),
        "longitude": rng.uniform(-180, 180, n),
        "mag": rng.exponential(scale=1.0, size=n) + 2.5,
        "type": "earthquake",
        "nst": rng.integers(5, 50, n),
    })

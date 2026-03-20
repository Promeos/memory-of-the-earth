"""Earthquake catalog data cleaning and preparation.

Provides utilities for loading raw CSV earthquake catalogs, cleaning and
deduplicating events, estimating magnitude of completeness, and identifying
temporal gaps.
"""

from pathlib import Path

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------


def load_raw_csvs(raw_dir):
    """Load and concatenate all CSV files from *raw_dir*.

    Parameters
    ----------
    raw_dir : str or Path
        Directory containing one or more ``.csv`` files.

    Returns
    -------
    pd.DataFrame
        Concatenated dataframe of all CSV contents.

    Raises
    ------
    FileNotFoundError
        If *raw_dir* does not exist or contains no CSV files.
    """
    raw_dir = Path(raw_dir)
    if not raw_dir.is_dir():
        raise FileNotFoundError(f"Directory not found: {raw_dir}")

    csv_files = sorted(raw_dir.glob("*.csv"))
    if not csv_files:
        raise FileNotFoundError(f"No CSV files found in {raw_dir}")

    frames = [pd.read_csv(f) for f in csv_files]
    return pd.concat(frames, ignore_index=True)


# ---------------------------------------------------------------------------
# Cleaning
# ---------------------------------------------------------------------------


def clean_catalog(df):
    """Clean a raw earthquake catalog dataframe.

    Steps performed:
    1. Filter rows to ``type == "earthquake"`` (if a *type* column exists).
    2. Parse the *time* column to ``datetime64`` (UTC).
    3. Drop rows with missing magnitude (*mag*).
    4. Sort by *time* ascending.

    Parameters
    ----------
    df : pd.DataFrame
        Raw catalog dataframe.  Expected columns include at least
        ``time`` and ``mag``.

    Returns
    -------
    pd.DataFrame
        Cleaned copy of the catalog.
    """
    df = df.copy()

    # Keep only earthquake-type rows when a type column is present.
    if "type" in df.columns:
        df = df[df["type"].str.lower().str.strip() == "earthquake"]

    # Parse time column.
    df["time"] = pd.to_datetime(df["time"], utc=True)

    # Drop rows without a magnitude value.
    df = df.dropna(subset=["mag"])

    # Sort chronologically.
    df = df.sort_values("time").reset_index(drop=True)

    return df


# ---------------------------------------------------------------------------
# Deduplication
# ---------------------------------------------------------------------------


def deduplicate(df, time_tol_sec=2, lat_tol=0.05, lon_tol=0.05):
    """Remove duplicate events reported by multiple networks.

    Events are grouped by rounding *time* (to the nearest ``time_tol_sec``
    seconds), *latitude*, and *longitude* to the specified tolerances.
    Within each group the event with the highest *nst* (station count) is
    kept.

    Parameters
    ----------
    df : pd.DataFrame
        Cleaned catalog with ``time``, ``latitude``, ``longitude``,
        ``mag``, and (optionally) ``nst`` columns.
    time_tol_sec : float
        Rounding tolerance for origin time in seconds.
    lat_tol : float
        Rounding tolerance for latitude in degrees.
    lon_tol : float
        Rounding tolerance for longitude in degrees.

    Returns
    -------
    pd.DataFrame
        Deduplicated catalog.
    """
    df = df.copy()

    # Create grouping keys by rounding coordinates and time.
    time_round = df["time"].dt.round(f"{time_tol_sec}s")
    lat_round = (df["latitude"] / lat_tol).round() * lat_tol
    lon_round = (df["longitude"] / lon_tol).round() * lon_tol

    df["_time_grp"] = time_round
    df["_lat_grp"] = lat_round
    df["_lon_grp"] = lon_round

    # If nst is missing fill with 0 so sorting still works.
    if "nst" not in df.columns:
        df["nst"] = 0
    df["nst"] = df["nst"].fillna(0)

    # Within each group keep the row with the largest station count.
    df = df.sort_values("nst", ascending=False)
    df = df.drop_duplicates(
        subset=["_time_grp", "_lat_grp", "_lon_grp"], keep="first"
    )

    # Clean up helper columns and re-sort.
    df = df.drop(columns=["_time_grp", "_lat_grp", "_lon_grp"])
    df = df.sort_values("time").reset_index(drop=True)

    return df


# ---------------------------------------------------------------------------
# Magnitude of completeness
# ---------------------------------------------------------------------------


def estimate_mc(magnitudes, method="maxc"):
    """Estimate magnitude of completeness.

    Uses the maximum curvature method (Wiemer & Wyss 2000): the bin centre
    of the most populated 0.1-unit magnitude bin, **plus 0.2**.

    Parameters
    ----------
    magnitudes : array-like
        Vector of earthquake magnitudes.
    method : str
        Estimation method.  Currently only ``"maxc"`` is supported.

    Returns
    -------
    float
        Estimated magnitude of completeness.

    Raises
    ------
    ValueError
        If an unsupported *method* is requested.
    """
    if method != "maxc":
        raise ValueError(f"Unsupported method: {method!r}. Use 'maxc'.")

    magnitudes = np.asarray(magnitudes, dtype=float)
    magnitudes = magnitudes[~np.isnan(magnitudes)]

    if len(magnitudes) == 0:
        return np.nan

    # Build frequency-magnitude distribution with 0.1-unit bins.
    bin_min = np.floor(magnitudes.min() * 10) / 10
    bin_max = np.ceil(magnitudes.max() * 10) / 10
    bins = np.arange(bin_min, bin_max + 0.1, 0.1)

    counts, edges = np.histogram(magnitudes, bins=bins)

    if len(counts) == 0 or counts.sum() == 0:
        return float(np.round(magnitudes.min() + 0.2, 1))

    # Bin centres.
    centres = (edges[:-1] + edges[1:]) / 2

    # MAXC: the centre of the bin with the maximum count + 0.2.
    mc = centres[np.argmax(counts)] + 0.2
    return float(np.round(mc, 1))


# ---------------------------------------------------------------------------
# Temporal gap detection
# ---------------------------------------------------------------------------


def identify_temporal_gaps(df, threshold_hours=24):
    """Identify temporal gaps in the catalog that exceed a threshold.

    Parameters
    ----------
    df : pd.DataFrame
        Catalog with a ``time`` column (datetime).
    threshold_hours : float
        Minimum gap duration in hours to report (default 24).

    Returns
    -------
    pd.DataFrame
        Columns: ``gap_start``, ``gap_end``, ``gap_hours``.
        One row per gap exceeding the threshold.
    """
    times = pd.to_datetime(df["time"]).sort_values().values
    deltas = np.diff(times).astype("timedelta64[s]").astype(float) / 3600.0

    mask = deltas > threshold_hours
    if not mask.any():
        return pd.DataFrame(columns=["gap_start", "gap_end", "gap_hours"])

    indices = np.where(mask)[0]
    gaps = []
    for i in indices:
        gaps.append({
            "gap_start": pd.Timestamp(times[i]),
            "gap_end": pd.Timestamp(times[i + 1]),
            "gap_hours": deltas[i],
        })

    return pd.DataFrame(gaps)

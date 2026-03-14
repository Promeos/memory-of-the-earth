"""Earthquake catalog data cleaning and preparation.

Provides utilities for loading raw CSV earthquake catalogs, cleaning and
deduplicating events, estimating magnitude of completeness, and writing
processed catalogs to Parquet format.
"""

from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

REGION_COLORS: dict[str, str] = {
    "oklahoma": "#E63946",
    "permian": "#F4A261",
    "socal": "#457B9D",
    "global": "#2A9D8F",
}

REGION_BOUNDS: dict[str, dict] = {
    "oklahoma": {
        "lat_range": (33.5, 37.5),
        "lon_range": (-100.0, -94.5),
        "min_mag": 1.0,
    },
    "permian": {
        "lat_range": (30.5, 33.5),
        "lon_range": (-105.0, -100.5),
        "min_mag": 1.0,
    },
    "socal": {
        "lat_range": (32.0, 36.5),
        "lon_range": (-121.0, -114.5),
        "min_mag": 1.0,
    },
    "global": {
        "lat_range": None,
        "lon_range": None,
        "min_mag": 2.5,
    },
}

# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------


def load_raw_csvs(raw_dir: str | Path) -> pd.DataFrame:
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


def clean_catalog(df: pd.DataFrame) -> pd.DataFrame:
    """Clean a raw earthquake catalog dataframe.

    Steps performed:
    1. Filter rows to ``type == "earthquake"`` (if a *type* column exists).
    2. Parse the *time* column to ``datetime64``.
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


def deduplicate(
    df: pd.DataFrame,
    time_tol_sec: float = 2,
    lat_tol: float = 0.05,
    lon_tol: float = 0.05,
) -> pd.DataFrame:
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
    df = df.drop_duplicates(subset=["_time_grp", "_lat_grp", "_lon_grp"], keep="first")

    # Clean up helper columns and re-sort.
    df = df.drop(columns=["_time_grp", "_lat_grp", "_lon_grp"])
    df = df.sort_values("time").reset_index(drop=True)

    return df


# ---------------------------------------------------------------------------
# Magnitude of completeness
# ---------------------------------------------------------------------------


def estimate_mc(
    magnitudes: np.ndarray | pd.Series,
    method: str = "maxc",
) -> float:
    """Estimate magnitude of completeness.

    Parameters
    ----------
    magnitudes : array-like
        Vector of earthquake magnitudes.
    method : str
        Estimation method.  Currently only ``"maxc"`` (maximum curvature)
        is supported.

    Returns
    -------
    float
        Estimated magnitude of completeness (bin centre of the most
        populated 0.1-unit magnitude bin).

    Raises
    ------
    ValueError
        If an unsupported *method* is requested.
    """
    if method != "maxc":
        raise ValueError(f"Unsupported method: {method!r}. Use 'maxc'.")

    magnitudes = np.asarray(magnitudes, dtype=float)
    magnitudes = magnitudes[~np.isnan(magnitudes)]

    # Build frequency-magnitude distribution with 0.1-unit bins.
    bin_min = np.floor(magnitudes.min() * 10) / 10
    bin_max = np.ceil(magnitudes.max() * 10) / 10
    bins = np.arange(bin_min, bin_max + 0.1, 0.1)

    counts, edges = np.histogram(magnitudes, bins=bins)

    # Bin centres.
    centres = (edges[:-1] + edges[1:]) / 2

    # MAXC: the centre of the bin with the maximum count.
    mc = centres[np.argmax(counts)]
    return float(np.round(mc, 1))


# ---------------------------------------------------------------------------
# Processing pipelines
# ---------------------------------------------------------------------------


def process_region(raw_dir: str | Path, output_path: str | Path) -> pd.DataFrame:
    """Load raw CSVs for a region, clean, deduplicate, and save as Parquet.

    Parameters
    ----------
    raw_dir : str or Path
        Directory containing raw CSV files for a single region.
    output_path : str or Path
        Destination Parquet file path.

    Returns
    -------
    pd.DataFrame
        The processed catalog.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    df = load_raw_csvs(raw_dir)
    df = clean_catalog(df)
    df = deduplicate(df)
    df.to_parquet(output_path, index=False)

    return df


def process_all_regions(base_dir: str | Path = "data") -> dict[str, pd.DataFrame]:
    """Process all four regions from raw CSVs to processed Parquet files.

    Expected directory layout::

        <base_dir>/
            raw/
                oklahoma/   *.csv
                permian/    *.csv
                socal/      *.csv
                global/     *.csv
            processed/
                oklahoma.parquet
                permian.parquet
                socal.parquet
                global.parquet

    Parameters
    ----------
    base_dir : str or Path
        Root data directory (default ``"data"``).

    Returns
    -------
    dict[str, pd.DataFrame]
        Mapping of region name to its processed catalog dataframe.
    """
    base_dir = Path(base_dir)
    results: dict[str, pd.DataFrame] = {}

    for region in REGION_BOUNDS:
        raw_dir = base_dir / "raw" / region
        output_path = base_dir / "processed" / f"{region}.parquet"
        results[region] = process_region(raw_dir, output_path)

    return results

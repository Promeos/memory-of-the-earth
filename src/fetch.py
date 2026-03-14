"""USGS ComCat earthquake API client with pagination and resume support.

Fetches earthquake catalogs from the USGS ComCat FDSN Event Web Service,
supporting monthly pagination, retry logic with exponential backoff,
and optional on-disk caching for resumable downloads.
"""

import io
import time
from datetime import datetime
from pathlib import Path

import pandas as pd
import requests

BASE_URL = "https://earthquake.usgs.gov/fdsnws/event/1/query"

REGIONS = {
    "global": {
        "min_mag": 2.5,
        "years": range(2000, 2026),
        "lat_range": None,
        "lon_range": None,
    },
    "oklahoma": {
        "min_mag": 1.0,
        "years": range(2000, 2026),
        "lat_range": (33.5, 37.5),
        "lon_range": (-100.0, -94.5),
    },
}


def fetch_month(year, month, min_magnitude=2.5, lat_range=None, lon_range=None):
    """Fetch one month of earthquake data from USGS ComCat.

    Parameters
    ----------
    year : int
        Calendar year.
    month : int
        Calendar month (1-12).
    min_magnitude : float, optional
        Minimum earthquake magnitude to include (default 2.5).
    lat_range : tuple of (float, float) or None, optional
        Latitude bounding box as (min_lat, max_lat).
    lon_range : tuple of (float, float) or None, optional
        Longitude bounding box as (min_lon, max_lon).

    Returns
    -------
    pd.DataFrame
        DataFrame with earthquake records for the requested month.
    """
    start = datetime(year, month, 1)
    if month == 12:
        end = datetime(year + 1, 1, 1)
    else:
        end = datetime(year, month + 1, 1)

    params = {
        "format": "csv",
        "starttime": start.strftime("%Y-%m-%dT%H:%M:%S"),
        "endtime": end.strftime("%Y-%m-%dT%H:%M:%S"),
        "minmagnitude": min_magnitude,
        "orderby": "time-asc",
    }

    if lat_range is not None:
        params["minlatitude"] = lat_range[0]
        params["maxlatitude"] = lat_range[1]
    if lon_range is not None:
        params["minlongitude"] = lon_range[0]
        params["maxlongitude"] = lon_range[1]

    resp = requests.get(BASE_URL, params=params, timeout=120)
    resp.raise_for_status()

    df = pd.read_csv(io.StringIO(resp.text))
    return df


def fetch_region(name, years, min_mag, lat_range=None, lon_range=None,
                 output_dir=None):
    """Fetch a full regional earthquake catalog with retry logic.

    Iterates over every month in the given year range, fetching CSV data
    from USGS ComCat.  Uses 3 retry attempts with exponential backoff on
    failure and a 0.5-second pause between successive requests.

    If ``output_dir`` is provided, each monthly CSV is saved to disk and
    already-existing files are skipped, enabling resumable downloads.

    Parameters
    ----------
    name : str
        Human-readable region name (used for logging and file naming).
    years : iterable of int
        Years to fetch (e.g. ``range(2000, 2026)``).
    min_mag : float
        Minimum magnitude threshold.
    lat_range : tuple of (float, float) or None, optional
        Latitude bounding box as (min_lat, max_lat).
    lon_range : tuple of (float, float) or None, optional
        Longitude bounding box as (min_lon, max_lon).
    output_dir : str or Path or None, optional
        Directory for saving raw monthly CSVs.  If ``None``, data is
        returned in memory only.

    Returns
    -------
    pd.DataFrame
        Concatenated DataFrame of all successfully fetched months.
    """
    if output_dir is not None:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

    frames = []
    max_retries = 3

    for year in years:
        for month in range(1, 13):
            label = f"{name} {year}-{month:02d}"

            # --- resumption: skip if file already on disk ---
            if output_dir is not None:
                csv_path = output_dir / f"{name}_{year}_{month:02d}.csv"
                if csv_path.exists():
                    print(f"[skip] {label} — file exists")
                    df = pd.read_csv(csv_path)
                    frames.append(df)
                    continue

            # --- fetch with retries ---
            df = None
            for attempt in range(1, max_retries + 1):
                try:
                    df = fetch_month(
                        year, month,
                        min_magnitude=min_mag,
                        lat_range=lat_range,
                        lon_range=lon_range,
                    )
                    break
                except Exception as exc:
                    wait = 2 ** attempt
                    print(
                        f"[retry] {label} attempt {attempt}/{max_retries} "
                        f"failed ({exc}), waiting {wait}s"
                    )
                    if attempt < max_retries:
                        time.sleep(wait)
                    else:
                        print(
                            f"[error] {label} — giving up after "
                            f"{max_retries} attempts"
                        )

            if df is not None:
                print(f"[ok]   {label} — {len(df)} events")
                if output_dir is not None:
                    df.to_csv(csv_path, index=False)
                frames.append(df)

            # polite pause between requests
            time.sleep(0.5)

    if not frames:
        return pd.DataFrame()

    return pd.concat(frames, ignore_index=True)


def fetch_all_regions(base_dir="data/raw"):
    """Fetch earthquake catalogs for all predefined regions.

    Convenience wrapper around :func:`fetch_region` that iterates over
    the standard study regions (Global M2.5+ and Oklahoma M1.0+) and
    saves raw CSVs under ``base_dir``.

    Parameters
    ----------
    base_dir : str or Path, optional
        Root directory for raw CSV output.  Each region gets its own
        subdirectory (default ``"data/raw"``).

    Returns
    -------
    dict[str, pd.DataFrame]
        Mapping of region name to its concatenated DataFrame.
    """
    base_dir = Path(base_dir)
    results = {}

    for region_name, cfg in REGIONS.items():
        print(f"\n{'='*60}")
        print(f"Fetching region: {region_name}")
        print(f"{'='*60}")

        region_dir = base_dir / region_name

        df = fetch_region(
            name=region_name,
            years=cfg["years"],
            min_mag=cfg["min_mag"],
            lat_range=cfg["lat_range"],
            lon_range=cfg["lon_range"],
            output_dir=region_dir,
        )

        results[region_name] = df
        print(f"[done] {region_name}: {len(df)} total events")

    return results

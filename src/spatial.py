"""Spatial utilities: gridding, haversine distances, and spatial queries.

Provides tools for dividing the globe into uniform grid cells, computing
great-circle distances, and performing radius-based spatial queries on
earthquake catalogs.
"""

import numpy as np
import pandas as pd


def haversine(lat1, lon1, lat2, lon2):
    """Compute great-circle distance between points using the Haversine formula.

    All inputs may be scalars or numpy arrays (vectorized).

    Parameters
    ----------
    lat1, lon1 : float or array-like
        Latitude and longitude of the first point(s) in decimal degrees.
    lat2, lon2 : float or array-like
        Latitude and longitude of the second point(s) in decimal degrees.

    Returns
    -------
    float or ndarray
        Distance(s) in kilometres.
    """
    R = 6371.0  # Earth mean radius in km

    lat1_r = np.radians(lat1)
    lat2_r = np.radians(lat2)
    dlat = np.radians(lat2 - lat1)
    dlon = np.radians(lon2 - lon1)

    a = (
        np.sin(dlat / 2.0) ** 2
        + np.cos(lat1_r) * np.cos(lat2_r) * np.sin(dlon / 2.0) ** 2
    )
    c = 2.0 * np.arctan2(np.sqrt(a), np.sqrt(1.0 - a))

    return R * c


def build_grid(lat_range=(-90, 90), lon_range=(-180, 180), cell_size=2.0):
    """Generate a uniform latitude/longitude grid.

    Parameters
    ----------
    lat_range : tuple of (float, float)
        (min_lat, max_lat) in degrees.
    lon_range : tuple of (float, float)
        (min_lon, max_lon) in degrees.
    cell_size : float
        Width and height of each cell in degrees (default 2.0).

    Returns
    -------
    list of dict
        Each dict has keys ``lat_min``, ``lat_max``, ``lon_min``, ``lon_max``,
        ``lat_center``, ``lon_center``.
    """
    cells = []
    lat_edges = np.arange(lat_range[0], lat_range[1], cell_size)
    lon_edges = np.arange(lon_range[0], lon_range[1], cell_size)

    for lat_min in lat_edges:
        lat_max = lat_min + cell_size
        for lon_min in lon_edges:
            lon_max = lon_min + cell_size
            cells.append({
                "lat_min": lat_min,
                "lat_max": lat_max,
                "lon_min": lon_min,
                "lon_max": lon_max,
                "lat_center": lat_min + cell_size / 2,
                "lon_center": lon_min + cell_size / 2,
            })

    return cells


def assign_cells(df, cell_size=2.0):
    """Assign each event to a spatial grid cell.

    Adds ``cell_lat`` and ``cell_lon`` columns representing the centre of the
    cell that contains each event.

    Parameters
    ----------
    df : pd.DataFrame
        Must contain ``latitude`` and ``longitude`` columns.
    cell_size : float
        Grid cell size in degrees (default 2.0).

    Returns
    -------
    pd.DataFrame
        Copy of *df* with ``cell_lat`` and ``cell_lon`` columns added.
    """
    df = df.copy()
    df["cell_lat"] = (np.floor(df["latitude"] / cell_size) * cell_size
                      + cell_size / 2)
    df["cell_lon"] = (np.floor(df["longitude"] / cell_size) * cell_size
                      + cell_size / 2)
    return df


def events_in_radius(df, lat, lon, radius_km):
    """Select events within a given radius of a reference point.

    Parameters
    ----------
    df : pd.DataFrame
        Must contain ``latitude`` and ``longitude`` columns.
    lat : float
        Reference latitude in degrees.
    lon : float
        Reference longitude in degrees.
    radius_km : float
        Search radius in kilometres.

    Returns
    -------
    pd.DataFrame
        Subset of *df* within the specified radius.
    """
    dists = haversine(df["latitude"].values, df["longitude"].values, lat, lon)
    return df[dists <= radius_km].copy()

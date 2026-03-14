"""
Aftershock cascade analysis using nearest-neighbor clustering.

Implements the rescaled space-time-magnitude distance metric from
Zaliapin & Ben-Zion (2013) to identify aftershock triggering trees
and compute cascade statistics for population comparison.
"""

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree, ConvexHull
from scipy.stats import mannwhitneyu


def haversine(lat1, lon1, lat2, lon2):
    """Compute great-circle distance between two points using the Haversine formula.

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
        Distance(s) in kilometers.
    """
    R = 6371.0  # Earth mean radius in km

    lat1_r = np.radians(lat1)
    lat2_r = np.radians(lat2)
    dlat = np.radians(lat2 - lat1)
    dlon = np.radians(lon2 - lon1)

    a = np.sin(dlat / 2.0) ** 2 + np.cos(lat1_r) * np.cos(lat2_r) * np.sin(dlon / 2.0) ** 2
    c = 2.0 * np.arctan2(np.sqrt(a), np.sqrt(1.0 - a))

    return R * c


def nearest_neighbor_distance(t_ij, r_ij, b, m_i, d_f=1.6):
    """Compute the rescaled space-time-magnitude nearest-neighbor distance.

    Uses the metric from Zaliapin & Ben-Zion (2013):
        eta = t_ij * r_ij^d_f * 10^(-b * m_i)

    Parameters
    ----------
    t_ij : float or array-like
        Interevent time in seconds between events i and j.
    r_ij : float or array-like
        Spatial distance in kilometers between events i and j.
    b : float
        Gutenberg-Richter b-value for the catalog.
    m_i : float or array-like
        Magnitude of the parent event i.
    d_f : float, optional
        Fractal dimension of the epicenter distribution. Default is 1.6.

    Returns
    -------
    float or ndarray
        Rescaled nearest-neighbor distance eta.
    """
    t_ij = np.asarray(t_ij, dtype=float)
    r_ij = np.asarray(r_ij, dtype=float)
    m_i = np.asarray(m_i, dtype=float)

    eta = t_ij * (r_ij ** d_f) * (10.0 ** (-b * m_i))
    return eta


def build_cascade_forest(catalog_df, b_value, min_mainshock_mag=4.0,
                         time_window_days=30, distance_window_km=200,
                         min_mag=2.5):
    """Build aftershock triggering trees using nearest-neighbor clustering.

    For each event j in the catalog, find its nearest parent i (occurring
    before j) using the rescaled space-time-magnitude distance. Events are
    grouped into cascades rooted at mainshocks.

    Parameters
    ----------
    catalog_df : pandas.DataFrame
        Earthquake catalog with columns: event_id, time (datetime64),
        latitude, longitude, magnitude, depth_km.
    b_value : float
        Gutenberg-Richter b-value for nearest-neighbor distance computation.
    min_mainshock_mag : float, optional
        Minimum magnitude to consider an event as a potential cascade root.
        Default is 4.0.
    time_window_days : float, optional
        Maximum look-back time in days when searching for a parent event.
        Default is 30.
    distance_window_km : float, optional
        Maximum spatial distance in km when searching for a parent event.
        Default is 200.
    min_mag : float, optional
        Minimum magnitude threshold; events below this are excluded.
        Default is 2.5.

    Returns
    -------
    pandas.DataFrame
        Columns: event_id, parent_id, cluster_id, generation.
        generation=0 for mainshocks (cascade roots).
    """
    # Filter to minimum magnitude and sort by time
    cat = catalog_df[catalog_df["magnitude"] >= min_mag].copy()
    cat = cat.sort_values("time").reset_index(drop=True)

    n = len(cat)
    if n == 0:
        return pd.DataFrame(columns=["event_id", "parent_id", "cluster_id", "generation"])

    # Precompute arrays for fast access
    event_ids = cat["event_id"].values
    times = cat["time"].values.astype("datetime64[s]").astype(np.float64)  # epoch seconds
    lats = cat["latitude"].values.astype(np.float64)
    lons = cat["longitude"].values.astype(np.float64)
    mags = cat["magnitude"].values.astype(np.float64)

    time_window_s = time_window_days * 86400.0

    # Build a cKDTree on (lat, lon) in radians for approximate spatial queries.
    # Convert distance_window_km to an angular threshold (radians) for the tree.
    angular_threshold = distance_window_km / 6371.0  # radians
    coords_rad = np.column_stack([np.radians(lats), np.radians(lons)])
    tree = cKDTree(coords_rad)

    parent_id = np.full(n, -1, dtype=np.int64)
    best_parent_idx = np.full(n, -1, dtype=np.int64)

    # For each event j, find its nearest parent i (i < j) within windows
    for j in range(1, n):
        # Spatial candidates from cKDTree
        candidate_indices = tree.query_ball_point(coords_rad[j], angular_threshold)
        # Keep only predecessors
        candidate_indices = [ci for ci in candidate_indices if ci < j]
        if not candidate_indices:
            continue

        ci = np.array(candidate_indices)

        # Time filter
        dt = times[j] - times[ci]
        time_mask = (dt > 0) & (dt <= time_window_s)
        ci = ci[time_mask]
        dt = dt[time_mask]

        if len(ci) == 0:
            continue

        # Compute exact haversine distances
        dists = haversine(lats[ci], lons[ci], lats[j], lons[j])

        # Distance filter
        dist_mask = dists <= distance_window_km
        ci = ci[dist_mask]
        dt = dt[dist_mask]
        dists = dists[dist_mask]

        if len(ci) == 0:
            continue

        # Compute rescaled nearest-neighbor distance for each candidate
        eta = nearest_neighbor_distance(dt, dists, b_value, mags[ci])

        # Select the candidate with the smallest eta
        best = np.argmin(eta)
        best_parent_idx[j] = ci[best]

    # Assign parent event_ids
    has_parent = best_parent_idx >= 0
    parent_id_arr = np.where(has_parent, event_ids[best_parent_idx], -1)

    # Build cluster assignments by walking parent links.
    # A root is any event with no parent OR a mainshock (mag >= min_mainshock_mag
    # that has no parent within the search window).
    cluster_ids = np.full(n, -1, dtype=np.int64)
    generations = np.full(n, -1, dtype=np.int64)

    # Identify roots: events with no parent
    roots = np.where(~has_parent)[0]
    cluster_counter = 0
    for r in roots:
        if mags[r] >= min_mainshock_mag:
            cluster_ids[r] = cluster_counter
            generations[r] = 0
            cluster_counter += 1

    # Build children map for forward traversal
    children = {i: [] for i in range(n)}
    for j in range(n):
        if best_parent_idx[j] >= 0:
            children[best_parent_idx[j]].append(j)

    # BFS from each root to assign cluster_id and generation
    for r in range(n):
        if cluster_ids[r] < 0 or generations[r] != 0:
            continue
        stack = [r]
        while stack:
            node = stack.pop()
            for child in children[node]:
                cluster_ids[child] = cluster_ids[node]
                generations[child] = generations[node] + 1
                stack.append(child)

    # Build result DataFrame (only events assigned to a cascade)
    result = pd.DataFrame({
        "event_id": event_ids,
        "parent_id": parent_id_arr,
        "cluster_id": cluster_ids,
        "generation": generations,
    })
    result = result[result["cluster_id"] >= 0].reset_index(drop=True)

    return result


def cascade_statistics(cascade_df, catalog_df):
    """Compute summary statistics for each aftershock cascade.

    Parameters
    ----------
    cascade_df : pandas.DataFrame
        Output of build_cascade_forest with columns: event_id, parent_id,
        cluster_id, generation.
    catalog_df : pandas.DataFrame
        Original catalog with columns: event_id, time (datetime64),
        latitude, longitude, magnitude, depth_km.

    Returns
    -------
    pandas.DataFrame
        One row per cascade (cluster_id) with columns:
        - cluster_id
        - tree_depth: maximum generation in the cascade
        - branching_ratio: mean number of children per non-leaf node
        - n_events: total number of events in the cascade
        - spatial_footprint_km2: approximate convex hull area
        - temporal_span_days: time from mainshock to last aftershock
        - magnitude_deficit: mainshock magnitude minus largest aftershock
          magnitude (Bath's law)
        - mean_depth_km: mean hypocentral depth of cascade events
    """
    merged = cascade_df.merge(catalog_df, on="event_id", how="left")

    records = []
    for cid, grp in merged.groupby("cluster_id"):
        n_events = len(grp)
        tree_depth = grp["generation"].max()

        # Branching ratio: count children per parent, then average over non-leaf nodes
        parent_counts = grp.loc[grp["generation"] > 0, "parent_id"].value_counts()
        if len(parent_counts) > 0:
            branching_ratio = parent_counts.mean()
        else:
            branching_ratio = 0.0

        # Spatial footprint via convex hull
        coords = grp[["latitude", "longitude"]].values
        if len(coords) >= 3:
            try:
                hull = ConvexHull(coords)
                # Approximate area: convert hull area from deg^2 to km^2
                # At the centroid latitude, 1 deg lat ~ 111 km,
                # 1 deg lon ~ 111 * cos(lat) km
                mean_lat = np.radians(coords[:, 0].mean())
                area_km2 = hull.volume * (111.0 ** 2) * np.cos(mean_lat)
            except Exception:
                area_km2 = 0.0
        elif len(coords) == 2:
            dist = haversine(coords[0, 0], coords[0, 1], coords[1, 0], coords[1, 1])
            area_km2 = 0.0  # degenerate; line, not area
        else:
            area_km2 = 0.0

        # Temporal span
        times = pd.to_datetime(grp["time"])
        temporal_span_days = (times.max() - times.min()).total_seconds() / 86400.0

        # Magnitude deficit (Bath's law)
        mainshock_mag = grp.loc[grp["generation"] == 0, "magnitude"].max()
        aftershock_mags = grp.loc[grp["generation"] > 0, "magnitude"]
        if len(aftershock_mags) > 0:
            magnitude_deficit = mainshock_mag - aftershock_mags.max()
        else:
            magnitude_deficit = np.nan

        # Mean depth
        mean_depth_km = grp["depth_km"].mean() if "depth_km" in grp.columns else np.nan

        records.append({
            "cluster_id": cid,
            "tree_depth": tree_depth,
            "branching_ratio": branching_ratio,
            "n_events": n_events,
            "spatial_footprint_km2": area_km2,
            "temporal_span_days": temporal_span_days,
            "magnitude_deficit": magnitude_deficit,
            "mean_depth_km": mean_depth_km,
        })

    return pd.DataFrame(records)


def compare_populations(stats_induced, stats_tectonic):
    """Compare induced and tectonic cascade populations using Mann-Whitney U tests.

    Tests whether the distributions of branching_ratio and tree_depth differ
    significantly between two populations of cascade statistics.

    Parameters
    ----------
    stats_induced : pandas.DataFrame
        Cascade statistics for induced seismicity (output of cascade_statistics).
    stats_tectonic : pandas.DataFrame
        Cascade statistics for tectonic seismicity (output of cascade_statistics).

    Returns
    -------
    dict
        Keys:
        - branching_ratio_U: Mann-Whitney U statistic for branching_ratio
        - branching_ratio_p: two-sided p-value for branching_ratio
        - tree_depth_U: Mann-Whitney U statistic for tree_depth
        - tree_depth_p: two-sided p-value for tree_depth
    """
    results = {}

    for metric in ["branching_ratio", "tree_depth"]:
        x = stats_induced[metric].dropna().values
        y = stats_tectonic[metric].dropna().values

        if len(x) == 0 or len(y) == 0:
            results[f"{metric}_U"] = np.nan
            results[f"{metric}_p"] = np.nan
        else:
            u_stat, p_val = mannwhitneyu(x, y, alternative="two-sided")
            results[f"{metric}_U"] = u_stat
            results[f"{metric}_p"] = p_val

    return results

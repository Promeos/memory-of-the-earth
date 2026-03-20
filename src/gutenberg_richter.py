"""Gutenberg-Richter b-value analysis and changepoint detection.

Provides maximum-likelihood b-value estimation (Aki, 1965), bootstrap
confidence intervals, rolling temporal analysis, and BIC-penalized
changepoint detection for seismicity rate and b-value time series.
"""

import numpy as np
import pandas as pd
from scipy import stats


# ---------------------------------------------------------------------------
# b-value estimation
# ---------------------------------------------------------------------------

def estimate_b_value(magnitudes, mc, delta_m=0.1):
    """Maximum-likelihood b-value estimation (Aki, 1965).

    Parameters
    ----------
    magnitudes : array_like
        Earthquake magnitudes.
    mc : float
        Completeness magnitude.  Only events with M >= mc are used.
    delta_m : float, optional
        Magnitude binning width (default 0.1).

    Returns
    -------
    b : float
        Estimated b-value.
    std_error : float
        Standard error via the Shi & Bolt (1982) formula.

    Raises
    ------
    ValueError
        If fewer than 2 events remain after filtering.
    """
    mags = np.asarray(magnitudes, dtype=float)
    mags = mags[mags >= mc]
    n = len(mags)
    if n < 2:
        raise ValueError(
            f"Need at least 2 events above mc={mc}, got {n}."
        )

    m_mean = mags.mean()
    denom = m_mean - mc + delta_m / 2.0
    if denom <= 0:
        raise ValueError(
            f"Mean magnitude {m_mean:.3f} too close to mc={mc} "
            f"(denom={denom:.6f}). Cannot estimate b-value."
        )
    b = np.log10(np.e) / denom

    # Shi & Bolt (1982) standard error
    variance_sum = np.sum((mags - m_mean) ** 2)
    std_error = 2.3 * b ** 2 * np.sqrt(variance_sum / (n * (n - 1)))

    return b, std_error


# ---------------------------------------------------------------------------
# Bootstrap confidence intervals
# ---------------------------------------------------------------------------

def bootstrap_b_value(magnitudes, mc, delta_m=0.1, n_bootstrap=1000):
    """Bootstrap confidence intervals on the b-value.

    Parameters
    ----------
    magnitudes : array_like
        Earthquake magnitudes.
    mc : float
        Completeness magnitude.
    delta_m : float, optional
        Magnitude binning width (default 0.1).
    n_bootstrap : int, optional
        Number of bootstrap resamples (default 1000).

    Returns
    -------
    b_mean : float
        Mean b-value across bootstrap samples.
    b_ci_low : float
        Lower bound of 95 % confidence interval.
    b_ci_high : float
        Upper bound of 95 % confidence interval.
    """
    mags = np.asarray(magnitudes, dtype=float)
    mags = mags[mags >= mc]
    n = len(mags)
    if n < 2:
        raise ValueError(
            f"Need at least 2 events above mc={mc}, got {n}."
        )

    rng = np.random.default_rng()
    b_samples = np.empty(n_bootstrap)

    for i in range(n_bootstrap):
        sample = rng.choice(mags, size=n, replace=True)
        m_mean = sample.mean()
        denom = m_mean - mc + delta_m / 2.0
        if denom <= 0:
            b_samples[i] = np.nan
        else:
            b_samples[i] = np.log10(np.e) / denom

    b_samples = b_samples[~np.isnan(b_samples)]
    b_mean = float(np.mean(b_samples))
    b_ci_low = float(np.percentile(b_samples, 2.5))
    b_ci_high = float(np.percentile(b_samples, 97.5))

    return b_mean, b_ci_low, b_ci_high


# ---------------------------------------------------------------------------
# Rolling b-value
# ---------------------------------------------------------------------------

def rolling_b_value(times, magnitudes, mc, window_size=200, step=50):
    """Compute b-values in sliding windows of fixed event count.

    Parameters
    ----------
    times : array_like
        Origin times (datetime-like or numeric).
    magnitudes : array_like
        Earthquake magnitudes corresponding to *times*.
    mc : float
        Completeness magnitude.
    window_size : int, optional
        Number of events per window (default 200).
    step : int, optional
        Step size in number of events (default 50).

    Returns
    -------
    pandas.DataFrame
        Columns: center_time, b_value, b_std, rate (events/day), n_events.
    """
    times = np.asarray(times)
    mags = np.asarray(magnitudes, dtype=float)

    # Sort by time
    order = np.argsort(times)
    times = times[order]
    mags = mags[order]

    results = []
    n_total = len(mags)

    for start in range(0, n_total - window_size + 1, step):
        end = start + window_size
        win_mags = mags[start:end]
        win_times = times[start:end]

        above_mc = win_mags[win_mags >= mc]
        n_above = len(above_mc)
        if n_above < 2:
            continue

        b, b_std = estimate_b_value(win_mags, mc)

        # Centre time
        center_time = win_times[window_size // 2]

        # Rate: events per day in the window span
        t_start = win_times[0]
        t_end = win_times[-1]
        try:
            span_days = (
                pd.Timestamp(t_end) - pd.Timestamp(t_start)
            ).total_seconds() / 86400.0
        except Exception:
            span_days = float(t_end - t_start)

        rate = window_size / span_days if span_days > 0 else np.nan

        results.append(
            {
                "center_time": center_time,
                "b_value": b,
                "b_std": b_std,
                "rate": rate,
                "n_events": n_above,
            }
        )

    return pd.DataFrame(results)


# ---------------------------------------------------------------------------
# Changepoint detection helpers
# ---------------------------------------------------------------------------

def segment_cost(data, start, end):
    """Gaussian cost of a segment for changepoint detection.

    Cost is defined as n * log(variance) for the segment data[start:end].

    Parameters
    ----------
    data : array_like
        1-D data array.
    start : int
        Start index (inclusive).
    end : int
        End index (exclusive).

    Returns
    -------
    float
        Segment cost.  Returns 0.0 if the segment has zero variance.
    """
    seg = np.asarray(data[start:end], dtype=float)
    n = len(seg)
    if n < 2:
        return 0.0
    var = np.var(seg, ddof=0)
    if var <= 0:
        return 0.0
    return n * np.log(var)


def detect_changepoints(data, min_size=20, penalty_factor=3.0):
    """BIC-penalized changepoint detection.

    If the ``ruptures`` package is available, PELT with an RBF kernel is
    used.  Otherwise an optimal-partitioning dynamic-programming algorithm
    is run with penalty = penalty_factor * log(n).

    Parameters
    ----------
    data : array_like
        1-D time series (e.g. rolling b-values or event rates).
    min_size : int, optional
        Minimum segment length (default 20).
    penalty_factor : float, optional
        Multiplier for log(n) penalty (default 3.0).

    Returns
    -------
    list of int
        Changepoint indices (positions where a new segment begins).
    """
    data = np.asarray(data, dtype=float)
    n = len(data)
    if n < 2 * min_size:
        return []

    penalty = penalty_factor * np.log(n)

    # Try ruptures first
    try:
        import ruptures as rpt

        model = rpt.Pelt(model="rbf", min_size=min_size).fit(data)
        bkps = model.predict(pen=penalty)
        # ruptures returns breakpoints with the last element == n; drop it
        changepoints = [bp for bp in bkps if bp < n]
        return sorted(changepoints)
    except ImportError:
        pass

    # -----------------------------------------------------------------
    # Fallback: optimal partitioning via dynamic programming
    # -----------------------------------------------------------------
    # cost_cache[i, j] = segment_cost(data, i, j)
    # F[t] = optimal cost for data[0:t]
    # F[0] = -penalty  (so the first segment has net cost = cost + penalty)
    F = np.full(n + 1, np.inf)
    F[0] = -penalty
    last_change = np.zeros(n + 1, dtype=int)

    for t in range(min_size, n + 1):
        for s in range(max(0, t - n), t - min_size + 1):
            seg_len = t - s
            if seg_len < min_size:
                continue
            cost = segment_cost(data, s, t)
            candidate = F[s] + cost + penalty
            if candidate < F[t]:
                F[t] = candidate
                last_change[t] = s

    # Back-track to recover changepoints
    changepoints = []
    idx = n
    while idx > 0:
        cp = last_change[idx]
        if cp > 0:
            changepoints.append(cp)
        idx = cp

    return sorted(changepoints)


# ---------------------------------------------------------------------------
# Kolmogorov-Smirnov test
# ---------------------------------------------------------------------------

def ks_test_magnitudes(mags1, mags2, mc):
    """Two-sample Kolmogorov-Smirnov test on magnitude distributions.

    Both samples are filtered to events with M >= mc before comparison.

    Parameters
    ----------
    mags1 : array_like
        First set of magnitudes.
    mags2 : array_like
        Second set of magnitudes.
    mc : float
        Completeness magnitude.

    Returns
    -------
    statistic : float
        KS test statistic.
    p_value : float
        Two-sided p-value.

    Raises
    ------
    ValueError
        If either filtered sample is empty.
    """
    m1 = np.asarray(mags1, dtype=float)
    m2 = np.asarray(mags2, dtype=float)
    m1 = m1[m1 >= mc]
    m2 = m2[m2 >= mc]

    if len(m1) == 0 or len(m2) == 0:
        raise ValueError(
            "Both magnitude samples must contain at least one event "
            f"above mc={mc}."
        )

    result = stats.ks_2samp(m1, m2)
    return result.statistic, result.pvalue


# ---------------------------------------------------------------------------
# Spatial gridding helpers
# ---------------------------------------------------------------------------

def compute_bvalue_grid(df, cell_size=2.0, min_events=50, delta_m=0.1):
    """Compute b-values on a spatial grid.

    Parameters
    ----------
    df : pd.DataFrame
        Catalog with ``latitude``, ``longitude``, and ``mag`` columns.
        Should already be cleaned and filtered.
    cell_size : float
        Grid cell size in degrees (default 2.0).
    min_events : int
        Minimum number of events above Mc for a valid b-estimate
        (default 50).
    delta_m : float
        Magnitude bin width (default 0.1).

    Returns
    -------
    pd.DataFrame
        One row per cell with columns: ``cell_lat``, ``cell_lon``,
        ``b_value``, ``b_std``, ``mc``, ``n_events``.
    """
    from .spatial import assign_cells
    from .catalog import estimate_mc as _estimate_mc

    df = assign_cells(df, cell_size=cell_size)

    rows = []
    for (clat, clon), grp in df.groupby(["cell_lat", "cell_lon"]):
        mags = grp["mag"].values
        mc = _estimate_mc(mags)
        above = mags[mags >= mc]
        n = len(above)

        if n < min_events:
            continue

        try:
            b, b_std = estimate_b_value(mags, mc, delta_m)
        except ValueError:
            continue

        rows.append({
            "cell_lat": clat,
            "cell_lon": clon,
            "b_value": b,
            "b_std": b_std,
            "mc": mc,
            "n_events": n,
        })

    return pd.DataFrame(rows)


def compute_cv_b(df, cell_size=2.0, window_years=3, stride_years=1,
                 min_events_per_window=50, min_events_total=200,
                 delta_m=0.1):
    """Compute the temporal coefficient of variation of b (CV_b) per cell.

    For each cell with enough data, computes b in rolling time windows
    and returns CV_b = std(b_windows) / mean(b_windows).

    Parameters
    ----------
    df : pd.DataFrame
        Catalog with ``latitude``, ``longitude``, ``mag``, and ``time``
        columns.
    cell_size : float
        Grid cell size in degrees (default 2.0).
    window_years : int
        Width of each rolling window in years (default 3).
    stride_years : int
        Step between windows in years (default 1).
    min_events_per_window : int
        Minimum events above Mc per window (default 50).
    min_events_total : int
        Minimum total events in a cell to attempt CV_b (default 200).
    delta_m : float
        Magnitude bin width (default 0.1).

    Returns
    -------
    pd.DataFrame
        One row per cell with columns: ``cell_lat``, ``cell_lon``,
        ``cv_b``, ``mean_b``, ``std_b``, ``n_windows``.
    """
    from .spatial import assign_cells
    from .catalog import estimate_mc as _estimate_mc

    df = assign_cells(df, cell_size=cell_size)
    df = df.copy()
    df["time"] = pd.to_datetime(df["time"], utc=True)

    rows = []
    for (clat, clon), grp in df.groupby(["cell_lat", "cell_lon"]):
        if len(grp) < min_events_total:
            continue

        mags = grp["mag"].values
        mc = _estimate_mc(mags)

        # Build rolling year windows
        year_min = grp["time"].dt.year.min()
        year_max = grp["time"].dt.year.max()

        b_values = []
        for y_start in range(year_min, year_max - window_years + 2,
                             stride_years):
            y_end = y_start + window_years
            t_start = pd.Timestamp(f"{y_start}-01-01", tz="UTC")
            t_end = pd.Timestamp(f"{y_end}-01-01", tz="UTC")

            mask = (grp["time"] >= t_start) & (grp["time"] < t_end)
            win_mags = grp.loc[mask, "mag"].values
            above = win_mags[win_mags >= mc]

            if len(above) < min_events_per_window:
                continue

            try:
                b, _ = estimate_b_value(win_mags, mc, delta_m)
                b_values.append(b)
            except ValueError:
                continue

        if len(b_values) < 2:
            continue

        b_arr = np.array(b_values)
        mean_b = np.mean(b_arr)
        std_b = np.std(b_arr, ddof=1)
        cv_b = std_b / mean_b if mean_b > 0 else np.nan

        rows.append({
            "cell_lat": clat,
            "cell_lon": clon,
            "cv_b": cv_b,
            "mean_b": mean_b,
            "std_b": std_b,
            "n_windows": len(b_values),
        })

    return pd.DataFrame(rows)

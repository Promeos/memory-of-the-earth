"""Shannon entropy of earthquake magnitude distributions.

Provides tools for computing the Shannon entropy of magnitude histograms
in rolling time windows, detecting entropy anomalies, and testing for
temporal associations between entropy drops and large earthquakes via
null models.
"""

import numpy as np
import pandas as pd


def shannon_entropy(magnitudes, mc, bin_width=0.5, max_mag=10.0):
    """Compute Shannon entropy of a magnitude distribution.

    H = -Σ p_i · log₂(p_i)

    Parameters
    ----------
    magnitudes : array-like
        Earthquake magnitudes (already filtered to events above *mc*).
    mc : float
        Magnitude of completeness (lower bound of first bin).
    bin_width : float
        Width of each magnitude bin (default 0.5).
    max_mag : float
        Upper edge of the last bin (default 7.5).

    Returns
    -------
    float
        Shannon entropy in bits.  Returns ``NaN`` if no events fall in
        any bin.
    """
    magnitudes = np.asarray(magnitudes, dtype=float)
    magnitudes = magnitudes[magnitudes >= mc]

    if len(magnitudes) == 0:
        return np.nan

    bins = np.arange(mc, max_mag + bin_width, bin_width)
    counts, _ = np.histogram(magnitudes, bins=bins)

    # Convert to probabilities
    total = counts.sum()
    if total == 0:
        return np.nan

    probs = counts / total
    # Remove zero-probability bins to avoid log(0)
    probs = probs[probs > 0]

    H = -np.sum(probs * np.log2(probs))
    return float(H)


def rolling_entropy(times, magnitudes, mc, window_days=90, stride_days=7,
                    min_events=100, bin_width=0.5):
    """Compute Shannon entropy in rolling time windows.

    Parameters
    ----------
    times : array-like of datetime
        Event origin times (sorted).
    magnitudes : array-like
        Event magnitudes corresponding to *times*.
    mc : float
        Magnitude of completeness.
    window_days : int
        Window width in days (default 90).
    stride_days : int
        Step between windows in days (default 7).
    min_events : int
        Minimum number of events per window for a valid entropy
        estimate (default 100).
    bin_width : float
        Magnitude bin width (default 0.5).

    Returns
    -------
    pd.DataFrame
        Columns: ``center_time``, ``H`` (entropy), ``n_events``.
    """
    times = pd.to_datetime(times)
    magnitudes = np.asarray(magnitudes, dtype=float)

    # Sort by time
    order = np.argsort(times)
    times = times[order]
    magnitudes = magnitudes[order]

    # Convert to int64 nanoseconds for fast searchsorted
    times_ns = times.astype(np.int64)

    t_min = times.min()
    t_max = times.max()
    window_ns = int(pd.Timedelta(days=window_days).value)
    stride_ns = int(pd.Timedelta(days=stride_days).value)
    half_window_ns = window_ns // 2

    # Pre-compute magnitude bins for entropy calculation
    bins = np.arange(mc, 10.0 + bin_width, bin_width)

    t_start_ns = int(t_min.value) if hasattr(t_min, 'value') else int(
        pd.Timestamp(t_min).value)
    t_max_ns = int(t_max.value) if hasattr(t_max, 'value') else int(
        pd.Timestamp(t_max).value)

    results = []
    t_ns = t_start_ns
    while t_ns + window_ns <= t_max_ns:
        t_end_ns = t_ns + window_ns

        # Use searchsorted for O(log n) window lookup
        i_start = np.searchsorted(times_ns, t_ns, side='left')
        i_end = np.searchsorted(times_ns, t_end_ns, side='left')

        mags_window = magnitudes[i_start:i_end]
        mags_above_mc = mags_window[mags_window >= mc]
        n = len(mags_above_mc)

        if n >= min_events:
            # Inline entropy calculation to avoid function call overhead
            counts, _ = np.histogram(mags_above_mc, bins=bins)
            total = counts.sum()
            if total == 0:
                H = np.nan
            else:
                probs = counts / total
                probs = probs[probs > 0]
                H = float(-np.sum(probs * np.log2(probs)))
        else:
            H = np.nan

        results.append({
            "center_time": pd.Timestamp(t_ns + half_window_ns),
            "H": H,
            "n_events": n,
        })

        t_ns += stride_ns

    return pd.DataFrame(results)


def detect_anomalies(entropy_df, percentile=5):
    """Identify entropy anomaly windows (drops below a percentile threshold).

    Parameters
    ----------
    entropy_df : pd.DataFrame
        Must contain ``center_time`` and ``H`` columns (output of
        :func:`rolling_entropy`).
    percentile : float
        Percentile threshold (default 5).  Windows with H below this
        percentile are flagged as anomalies.

    Returns
    -------
    pd.DataFrame
        Subset of *entropy_df* where H is below the threshold.
    """
    valid = entropy_df.dropna(subset=["H"])
    threshold = np.percentile(valid["H"], percentile)
    return valid[valid["H"] <= threshold].copy()


def null_model_test(times, magnitudes, mc, large_event_times,
                    window_days=90, stride_days=7, min_events=100,
                    association_days=60, percentile=5, n_shuffles=200,
                    seed=42):
    """Test whether entropy anomalies are associated with large earthquakes.

    Shuffles event magnitudes (preserving times and locations) to build a
    null distribution of hit rates, then compares the observed hit rate.

    Parameters
    ----------
    times : array-like of datetime
        Event origin times.
    magnitudes : array-like
        Event magnitudes.
    mc : float
        Magnitude of completeness.
    large_event_times : array-like of datetime
        Times of large earthquakes to test association against.
    window_days, stride_days, min_events : int
        Parameters for :func:`rolling_entropy`.
    association_days : int
        Time window (days) after an anomaly within which a large event
        counts as a "hit" (default 180).
    percentile : float
        Anomaly threshold percentile (default 5).
    n_shuffles : int
        Number of null-model iterations (default 1000).
    seed : int
        Random seed for reproducibility.

    Returns
    -------
    dict
        Keys: ``observed_hit_rate``, ``null_hit_rates`` (array),
        ``p_value`` (fraction of null rates >= observed).
    """
    large_event_times = pd.to_datetime(large_event_times)
    # Pre-convert to int64 nanoseconds for fast comparison
    large_ns = large_event_times.astype(np.int64)
    large_ns_sorted = np.sort(large_ns)
    assoc_ns = int(pd.Timedelta(days=association_days).value)

    def _hit_rate(mags):
        ent = rolling_entropy(times, mags, mc, window_days, stride_days,
                              min_events)
        anomalies = detect_anomalies(ent, percentile)
        if len(anomalies) == 0:
            return 0.0

        # Vectorized hit calculation using searchsorted
        anom_times_ns = pd.to_datetime(
            anomalies["center_time"]).astype(np.int64).values
        hits = 0
        for t_ns in anom_times_ns:
            # Check if any large event falls in [t_ns, t_ns + assoc_ns]
            idx = np.searchsorted(large_ns_sorted, t_ns, side='left')
            if idx < len(large_ns_sorted) and large_ns_sorted[idx] <= t_ns + assoc_ns:
                hits += 1

        return hits / len(anomalies)

    # Observed hit rate
    magnitudes = np.asarray(magnitudes, dtype=float)
    observed = _hit_rate(magnitudes)

    # Null model: shuffle magnitudes
    rng = np.random.default_rng(seed)
    null_rates = np.empty(n_shuffles)
    for i in range(n_shuffles):
        shuffled = rng.permutation(magnitudes)
        null_rates[i] = _hit_rate(shuffled)

    p_value = np.mean(null_rates >= observed)

    return {
        "observed_hit_rate": observed,
        "null_hit_rates": null_rates,
        "p_value": p_value,
    }

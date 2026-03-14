"""Interevent time analysis: distribution fitting, DFA, and autocorrelation.

Provides tools for characterising the temporal structure of earthquake
catalogues, including coefficient-of-variation tracking, detrended
fluctuation analysis (Peng et al., 1994), and classical autocorrelation.
"""

import numpy as np
import pandas as pd
from scipy import stats


# ---------------------------------------------------------------------------
# Interevent times
# ---------------------------------------------------------------------------

def compute_interevent_times(times):
    """Compute time differences between consecutive events in seconds.

    Parameters
    ----------
    times : array-like of datetime or pandas.Timestamp
        Sorted event origin times.

    Returns
    -------
    numpy.ndarray
        Interevent times in seconds (length = len(times) - 1).
    """
    times = pd.to_datetime(times)
    deltas = np.diff(times)
    # Convert timedelta64 to float seconds
    return deltas.astype("timedelta64[ns]").astype(np.float64) / 1e9


# ---------------------------------------------------------------------------
# Distribution fitting
# ---------------------------------------------------------------------------

def fit_interevent_distributions(iet):
    """Fit four parametric models to interevent times via MLE.

    Models
    ------
    * Exponential  (``scipy.stats.expon``)
    * Gamma        (``scipy.stats.gamma``)
    * Log-normal   (``scipy.stats.lognorm``)
    * Weibull      (``scipy.stats.weibull_min``)

    Parameters
    ----------
    iet : array-like
        Positive interevent times (seconds).

    Returns
    -------
    dict
        Keys are distribution names mapping to dicts with ``params``,
        ``loglik``, ``aic``, ``bic``.  Additional keys ``best_aic`` and
        ``best_bic`` give the name of the winning model under each criterion.
    """
    iet = np.asarray(iet, dtype=np.float64)
    n = len(iet)

    distributions = {
        "exponential": stats.expon,
        "gamma": stats.gamma,
        "lognormal": stats.lognorm,
        "weibull": stats.weibull_min,
    }

    results = {}

    for name, dist in distributions.items():
        params = dist.fit(iet)
        loglik = np.sum(dist.logpdf(iet, *params))
        k = len(params)
        aic = 2 * k - 2 * loglik
        bic = k * np.log(n) - 2 * loglik
        results[name] = {
            "params": params,
            "loglik": loglik,
            "aic": aic,
            "bic": bic,
        }

    results["best_aic"] = min(
        (name for name in distributions), key=lambda name: results[name]["aic"]
    )
    results["best_bic"] = min(
        (name for name in distributions), key=lambda name: results[name]["bic"]
    )

    return results


# ---------------------------------------------------------------------------
# Rolling coefficient of variation
# ---------------------------------------------------------------------------

def rolling_cv(times, window_size=200, step=50):
    """Rolling coefficient of variation (std / mean) of interevent times.

    Parameters
    ----------
    times : array-like of datetime or pandas.Timestamp
        Sorted event origin times.
    window_size : int
        Number of interevent times per window.
    step : int
        Step size (in number of interevent times) between windows.

    Returns
    -------
    pandas.DataFrame
        Columns: ``center_time``, ``cv``.
    """
    times = pd.to_datetime(times)
    iet = compute_interevent_times(times)

    centres = []
    cvs = []
    for start in range(0, len(iet) - window_size + 1, step):
        window = iet[start : start + window_size]
        mean = np.mean(window)
        if mean == 0:
            cv = np.nan
        else:
            cv = np.std(window, ddof=1) / mean
        # Centre time: midpoint of the window in the original times array
        center_idx = start + window_size // 2
        centres.append(times[center_idx])
        cvs.append(cv)

    return pd.DataFrame({"center_time": centres, "cv": cvs})


# ---------------------------------------------------------------------------
# Detrended Fluctuation Analysis
# ---------------------------------------------------------------------------

def dfa(x, scales=None):
    """Detrended Fluctuation Analysis (Peng et al., 1994).

    Parameters
    ----------
    x : array-like
        Time series (e.g. interevent times).
    scales : array-like of int, optional
        Box sizes.  If *None*, ~30 log-spaced values from 10 to
        ``len(x) // 4`` are used.

    Returns
    -------
    scales : numpy.ndarray
        Box sizes actually used.
    fluctuations : numpy.ndarray
        RMS fluctuation *F(n)* for each scale.
    H : float
        Estimated Hurst (scaling) exponent — slope of
        log *F(n)* vs log *n*.
    """
    x = np.asarray(x, dtype=np.float64)
    N = len(x)

    # Cumulative sum of demeaned series (the "profile")
    profile = np.cumsum(x - np.mean(x))

    if scales is None:
        max_scale = max(N // 4, 11)
        scales = np.unique(
            np.logspace(np.log10(10), np.log10(max_scale), num=30).astype(int)
        )
    else:
        scales = np.asarray(scales, dtype=int)

    fluctuations = np.zeros(len(scales))

    for i, n in enumerate(scales):
        num_segments = N // n
        if num_segments == 0:
            fluctuations[i] = np.nan
            continue

        rms_segments = np.zeros(num_segments)
        for j in range(num_segments):
            segment = profile[j * n : (j + 1) * n]
            # Linear detrend: fit and subtract a line
            t = np.arange(n)
            coeffs = np.polyfit(t, segment, 1)
            trend = np.polyval(coeffs, t)
            rms_segments[j] = np.sqrt(np.mean((segment - trend) ** 2))

        fluctuations[i] = np.mean(rms_segments)

    # Remove any NaN scales
    valid = ~np.isnan(fluctuations)
    scales = scales[valid]
    fluctuations = fluctuations[valid]

    # Log-log fit for the Hurst exponent
    if len(scales) >= 2:
        coeffs = np.polyfit(np.log(scales), np.log(fluctuations), 1)
        H = coeffs[0]
    else:
        H = np.nan

    return scales, fluctuations, H


# ---------------------------------------------------------------------------
# Autocorrelation
# ---------------------------------------------------------------------------

def compute_autocorrelation(x, max_lag=500):
    """Normalised autocorrelation function of a time series.

    Parameters
    ----------
    x : array-like
        Input time series.
    max_lag : int
        Maximum lag to compute (inclusive).

    Returns
    -------
    numpy.ndarray
        ACF values from lag 0 to *max_lag*.
    """
    x = np.asarray(x, dtype=np.float64)
    x = x - np.mean(x)
    n = len(x)
    max_lag = min(max_lag, n - 1)

    # Full correlation via numpy.correlate
    full_corr = np.correlate(x, x, mode="full")
    # Normalise by the zero-lag value
    mid = n - 1  # index of zero lag in full_corr
    acf = full_corr[mid : mid + max_lag + 1] / full_corr[mid]

    return acf


# ---------------------------------------------------------------------------
# Convenience wrapper
# ---------------------------------------------------------------------------

def hurst_exponent(interevent_times):
    """Estimate the Hurst exponent of interevent times via DFA.

    Parameters
    ----------
    interevent_times : array-like
        Positive interevent times (seconds).

    Returns
    -------
    float
        Hurst exponent *H*.
    """
    _, _, H = dfa(interevent_times)
    return H

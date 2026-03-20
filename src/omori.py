"""Modified Omori Law fitting and stress-recovery analysis.

Provides tools for fitting the Modified Omori Law to aftershock rate decay,
computing b-value recovery fractions, and estimating recovery time constants.
"""

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit


# ---------------------------------------------------------------------------
# Modified Omori Law
# ---------------------------------------------------------------------------

def _omori_rate(t, K, c, p):
    """Modified Omori Law: R(t) = K * (t + c)^(-p)."""
    return K * (t + c) ** (-p)


def fit_omori(event_times, mainshock_time, window_days=30, stride_days=10):
    """Fit the Modified Omori Law to aftershock rate decay.

    Computes event rates in sliding windows after the mainshock, then fits
    R(t) = K * (t + c)^(-p) to the rate-vs-time curve.

    Parameters
    ----------
    event_times : array-like of datetime
        Origin times of aftershock events (already filtered to the
        aftershock zone).
    mainshock_time : datetime-like
        Origin time of the mainshock.
    window_days : float
        Width of each rate-estimation window in days.
    stride_days : float
        Step between successive windows in days.

    Returns
    -------
    dict
        Keys: ``K``, ``c``, ``p`` (fitted parameters), ``t_centers``
        (window centre times in days after mainshock), ``rates``
        (observed rates in events/day), ``r_squared`` (goodness of fit).
        Returns ``None`` if fitting fails.
    """
    event_times = pd.to_datetime(event_times)
    mainshock_time = pd.Timestamp(mainshock_time)

    # Days after mainshock for each event
    dt = (event_times - mainshock_time).dt.total_seconds() / 86400.0
    dt = dt[dt > 0].values

    if len(dt) < 10:
        return None

    max_day = dt.max()
    t_starts = np.arange(0, max_day - window_days, stride_days)

    t_centers = []
    rates = []
    for t0 in t_starts:
        t1 = t0 + window_days
        count = np.sum((dt >= t0) & (dt < t1))
        rates.append(count / window_days)
        t_centers.append(t0 + window_days / 2)

    t_centers = np.array(t_centers)
    rates = np.array(rates)

    # Remove zero-rate windows for fitting
    mask = rates > 0
    if mask.sum() < 3:
        return None

    try:
        popt, pcov = curve_fit(
            _omori_rate,
            t_centers[mask],
            rates[mask],
            p0=[rates[mask][0], 1.0, 1.0],
            bounds=([0, 0.001, 0.01], [np.inf, 1000, 5]),
            maxfev=5000,
        )
    except (RuntimeError, ValueError):
        return None

    K, c, p = popt

    # R² goodness of fit
    predicted = _omori_rate(t_centers[mask], K, c, p)
    ss_res = np.sum((rates[mask] - predicted) ** 2)
    ss_tot = np.sum((rates[mask] - np.mean(rates[mask])) ** 2)
    r_squared = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0

    return {
        "K": K,
        "c": c,
        "p": p,
        "K_std": np.sqrt(pcov[0, 0]) if pcov is not None else np.nan,
        "c_std": np.sqrt(pcov[1, 1]) if pcov is not None else np.nan,
        "p_std": np.sqrt(pcov[2, 2]) if pcov is not None else np.nan,
        "t_centers": t_centers,
        "rates": rates,
        "r_squared": r_squared,
    }


# ---------------------------------------------------------------------------
# b-value recovery
# ---------------------------------------------------------------------------

def compute_recovery_fraction(b_timeseries, b_baseline):
    """Compute the b-value recovery fraction over time.

    f(t) = 1 - |b(t) - b_baseline| / |b(t=0) - b_baseline|

    Parameters
    ----------
    b_timeseries : array-like
        b-values in successive time windows after the mainshock.
    b_baseline : float
        Pre-mainshock baseline b-value.

    Returns
    -------
    numpy.ndarray
        Recovery fraction at each time step.  Values near 1 mean recovery
        is nearly complete; values near 0 mean the perturbation persists.
    """
    b_timeseries = np.asarray(b_timeseries, dtype=float)
    initial_perturbation = abs(b_timeseries[0] - b_baseline)

    if initial_perturbation < 1e-6:
        # No measurable perturbation
        return np.ones(len(b_timeseries))

    fractions = 1.0 - np.abs(b_timeseries - b_baseline) / initial_perturbation
    return fractions


def fit_recovery_exponential(recovery_fractions, t_days):
    """Fit an exponential recovery model to the recovery fraction.

    f(t) = 1 - exp(-t / τ_b)

    Parameters
    ----------
    recovery_fractions : array-like
        Recovery fraction values (from :func:`compute_recovery_fraction`).
    t_days : array-like
        Time in days corresponding to each recovery fraction.

    Returns
    -------
    dict or None
        Keys: ``tau_b`` (recovery time constant in days), ``r_squared``.
        Returns ``None`` if fitting fails.
    """
    f = np.asarray(recovery_fractions, dtype=float)
    t = np.asarray(t_days, dtype=float)

    # Remove NaN and non-finite values
    mask = np.isfinite(f) & np.isfinite(t) & (t > 0)
    f = f[mask]
    t = t[mask]

    if len(f) < 3:
        return None

    def _model(t, tau):
        return 1.0 - np.exp(-t / tau)

    try:
        popt, _ = curve_fit(
            _model, t, f,
            p0=[100.0],
            bounds=([1.0], [5000.0]),
            maxfev=5000,
        )
    except (RuntimeError, ValueError):
        return None

    tau_b = popt[0]

    predicted = _model(t, tau_b)
    ss_res = np.sum((f - predicted) ** 2)
    ss_tot = np.sum((f - np.mean(f)) ** 2)
    r_squared = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0

    return {"tau_b": tau_b, "r_squared": r_squared}

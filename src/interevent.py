"""Inter-event time computation, distribution fitting, and regime classification.

Provides tools for computing inter-event times from earthquake catalogs,
fitting parametric distributions via MLE, classifying temporal regimes
(Poisson-like, clustered, quasi-periodic), and tracking rolling statistics.
"""

import numpy as np
import pandas as pd
from scipy import stats


# ---------------------------------------------------------------------------
# Inter-event times
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
        Inter-event times in seconds (length = len(times) - 1).
    """
    times = pd.to_datetime(times)
    times = times.sort_values() if hasattr(times, 'sort_values') else np.sort(times)
    deltas = np.diff(times)
    return deltas.astype("timedelta64[ns]").astype(np.float64) / 1e9


# ---------------------------------------------------------------------------
# Distribution fitting
# ---------------------------------------------------------------------------

def fit_distributions(iet):
    """Fit four parametric models to inter-event times via MLE.

    Models
    ------
    * Exponential  (``scipy.stats.expon``)
    * Gamma        (``scipy.stats.gamma``)
    * Log-normal   (``scipy.stats.lognorm``)
    * Weibull      (``scipy.stats.weibull_min``)

    Parameters
    ----------
    iet : array-like
        Positive inter-event times (seconds or hours).

    Returns
    -------
    dict
        Keys are distribution names mapping to dicts with ``params``,
        ``loglik``, ``aic``, ``bic``.  Additional keys ``best_aic`` and
        ``best_bic`` give the name of the winning model under each criterion.
    """
    iet = np.asarray(iet, dtype=np.float64)
    iet = iet[iet > 0]
    n = len(iet)

    distributions = {
        "exponential": stats.expon,
        "gamma": stats.gamma,
        "lognormal": stats.lognorm,
        "weibull": stats.weibull_min,
    }

    results = {}

    for name, dist in distributions.items():
        params = dist.fit(iet, floc=0)
        loglik = np.sum(dist.logpdf(iet, *params))
        k = len(params) - 1  # exclude fixed loc parameter
        aic = 2 * k - 2 * loglik
        bic = k * np.log(n) - 2 * loglik
        results[name] = {
            "params": params,
            "loglik": loglik,
            "aic": aic,
            "bic": bic,
        }

    results["best_aic"] = min(
        distributions, key=lambda name: results[name]["aic"]
    )
    results["best_bic"] = min(
        distributions, key=lambda name: results[name]["bic"]
    )

    return results


# ---------------------------------------------------------------------------
# Regime classification
# ---------------------------------------------------------------------------

def classify_regime(fit_result):
    """Assign a distributional regime class based on AIC model selection.

    Classification rules (from README §3.4):

    * **Poisson-like:** Exponential wins, or Weibull wins with k in [0.9, 1.1].
    * **Clustered-like:** Weibull with k < 0.9, or gamma with alpha < 0.9.
    * **Quasi-periodic-like:** Weibull with k > 1.1.
    * **Ambiguous:** No model achieves ΔAIC > 10 over its nearest competitor.

    Parameters
    ----------
    fit_result : dict
        Output of :func:`fit_distributions`.

    Returns
    -------
    str
        One of ``"poisson"``, ``"clustered"``, ``"quasi_periodic"``, or
        ``"ambiguous"``.
    """
    dist_names = ["exponential", "gamma", "lognormal", "weibull"]
    aics = {name: fit_result[name]["aic"] for name in dist_names}
    sorted_aics = sorted(aics.items(), key=lambda x: x[1])
    best_name, best_aic = sorted_aics[0]
    second_aic = sorted_aics[1][1]
    delta_aic = second_aic - best_aic

    # Ambiguous if ΔAIC ≤ 10
    if delta_aic <= 10:
        return "ambiguous"

    if best_name == "exponential":
        return "poisson"

    if best_name == "weibull":
        # Weibull shape parameter k is the first fitted parameter
        k = fit_result["weibull"]["params"][0]
        if 0.9 <= k <= 1.1:
            return "poisson"
        elif k < 0.9:
            return "clustered"
        else:
            return "quasi_periodic"

    if best_name == "gamma":
        # Gamma shape parameter alpha is the first fitted parameter
        alpha = fit_result["gamma"]["params"][0]
        if alpha < 0.9:
            return "clustered"
        return "poisson"

    if best_name == "lognormal":
        return "clustered"

    return "ambiguous"


# ---------------------------------------------------------------------------
# Rolling coefficient of variation
# ---------------------------------------------------------------------------

def rolling_cv(times, window_size=200, step=50):
    """Rolling coefficient of variation (std / mean) of inter-event times.

    Parameters
    ----------
    times : array-like of datetime or pandas.Timestamp
        Sorted event origin times.
    window_size : int
        Number of inter-event times per window.
    step : int
        Step size (in number of inter-event times) between windows.

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
        center_idx = start + window_size // 2
        centres.append(times[center_idx])
        cvs.append(cv)

    return pd.DataFrame({"center_time": centres, "cv": cvs})

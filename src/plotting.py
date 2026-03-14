"""Shared visualization utilities for the Memory of the Earth project.

Provides consistent styling, figure saving, and specialized plot functions
for global maps, recovery curves, regime classification maps, Oklahoma
timelines, and entropy time series.
"""

import warnings
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import numpy as np


# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

DEFAULT_FIGSIZE = (12, 7)
DEFAULT_DPI = 300


# ---------------------------------------------------------------------------
# Style helpers
# ---------------------------------------------------------------------------

def setup_style():
    """Apply a consistent matplotlib style for the project.

    Tries ``seaborn-v0_8-whitegrid`` first (matplotlib >= 3.6), then falls
    back to ``seaborn-whitegrid`` (older matplotlib), and finally uses the
    default style.
    """
    style_applied = False
    for style_name in ("seaborn-v0_8-whitegrid", "seaborn-whitegrid"):
        try:
            plt.style.use(style_name)
            style_applied = True
            break
        except OSError:
            continue

    if not style_applied:
        plt.style.use("default")

    plt.rcParams.update({
        "figure.figsize": DEFAULT_FIGSIZE,
        "figure.dpi": DEFAULT_DPI,
        "axes.titlesize": 14,
        "axes.labelsize": 12,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
    })


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------

def save_figure(fig, name, figures_dir="figures"):
    """Save *fig* as a PNG file inside *figures_dir*.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        The figure to save.
    name : str
        Base filename (with or without ``.png`` extension).
    figures_dir : str or Path
        Directory in which to save.  Created if it does not exist.

    Returns
    -------
    pathlib.Path
        Absolute path to the saved file.
    """
    figures_dir = Path(figures_dir)
    figures_dir.mkdir(parents=True, exist_ok=True)

    if not name.endswith(".png"):
        name = f"{name}.png"

    path = figures_dir / name
    fig.savefig(path, dpi=DEFAULT_DPI, bbox_inches="tight")
    return path.resolve()


# ---------------------------------------------------------------------------
# Global map plots
# ---------------------------------------------------------------------------

def plot_global_map(lats, lons, values, cmap="RdYlBu_r", vmin=None,
                    vmax=None, label="", title="", cell_size=2.0, ax=None):
    """Plot values on a global map using a scatter-style grid.

    Parameters
    ----------
    lats : array-like
        Cell centre latitudes.
    lons : array-like
        Cell centre longitudes.
    values : array-like
        Values to plot (one per cell).
    cmap : str
        Matplotlib colormap name.
    vmin, vmax : float or None
        Color scale limits.
    label : str
        Colorbar label.
    title : str
        Figure title.
    cell_size : float
        Grid cell size in degrees (for marker sizing).
    ax : matplotlib.axes.Axes or None
        Axes to draw on.  A new figure is created when *None*.

    Returns
    -------
    matplotlib.axes.Axes
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(14, 7))

    lats = np.asarray(lats)
    lons = np.asarray(lons)
    values = np.asarray(values)

    sc = ax.scatter(
        lons, lats, c=values, cmap=cmap, vmin=vmin, vmax=vmax,
        s=8, marker="s", edgecolors="none", alpha=0.85,
    )

    cbar = plt.colorbar(sc, ax=ax, pad=0.02, shrink=0.7)
    cbar.set_label(label)

    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    ax.set_xlabel("Longitude (°)")
    ax.set_ylabel("Latitude (°)")
    ax.set_title(title)
    ax.set_aspect("equal")

    return ax


def plot_bvalue_volatility_map(lats, lons, cv_values, ax=None):
    """Plot global b-value volatility (CV_b) map.

    Green (stable) → Yellow → Red (volatile).

    Parameters
    ----------
    lats, lons : array-like
        Cell centre coordinates.
    cv_values : array-like
        CV_b values per cell.
    ax : matplotlib.axes.Axes or None

    Returns
    -------
    matplotlib.axes.Axes
    """
    return plot_global_map(
        lats, lons, cv_values,
        cmap="RdYlGn_r", vmin=0, vmax=0.2,
        label="CV_b (b-value volatility)",
        title="Global b-value Temporal Volatility",
        ax=ax,
    )


# ---------------------------------------------------------------------------
# Recovery curves
# ---------------------------------------------------------------------------

def plot_recovery_gallery(recovery_data, n_rows=3, n_cols=3):
    """Plot a gallery of b-value and rate recovery curves.

    Parameters
    ----------
    recovery_data : list of dict
        Each dict should have keys ``name``, ``t_days``, ``b_values``,
        ``b_baseline``, ``rates`` (optional).
    n_rows, n_cols : int
        Grid layout.

    Returns
    -------
    matplotlib.figure.Figure
    """
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 4 * n_rows))
    axes_flat = axes.flatten()

    for i, ax in enumerate(axes_flat):
        if i >= len(recovery_data):
            ax.set_visible(False)
            continue

        rd = recovery_data[i]
        ax.plot(rd["t_days"], rd["b_values"], "b-", linewidth=1.5,
                label="b(t)")
        ax.axhline(rd["b_baseline"], color="gray", linestyle="--",
                   linewidth=1, label="baseline")
        ax.set_title(rd["name"], fontsize=11)
        ax.set_xlabel("Days after mainshock")
        ax.set_ylabel("b-value")
        ax.legend(fontsize=8)

    fig.suptitle("Stress Recovery Curves", fontsize=14, fontweight="bold")
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# Regime classification map
# ---------------------------------------------------------------------------

def plot_regime_map(lats, lons, classes, ax=None):
    """Plot global distributional class map.

    Colors: Poisson-like = blue, Clustered-like = red,
    Quasi-periodic-like = green, Ambiguous = gray.

    Parameters
    ----------
    lats, lons : array-like
        Cell centre coordinates.
    classes : array-like of str
        Regime class per cell.
    ax : matplotlib.axes.Axes or None

    Returns
    -------
    matplotlib.axes.Axes
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(14, 7))

    color_map = {
        "poisson": "#457B9D",
        "clustered": "#E63946",
        "quasi_periodic": "#2A9D8F",
        "ambiguous": "#AAAAAA",
    }

    colors = [color_map.get(c, "#AAAAAA") for c in classes]

    ax.scatter(lons, lats, c=colors, s=8, marker="s",
               edgecolors="none", alpha=0.85)

    # Legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor="#457B9D", label="Poisson-like"),
        Patch(facecolor="#E63946", label="Clustered-like"),
        Patch(facecolor="#2A9D8F", label="Quasi-periodic-like"),
        Patch(facecolor="#AAAAAA", label="Ambiguous"),
    ]
    ax.legend(handles=legend_elements, loc="lower left", fontsize=9)

    ax.set_xlim(-180, 180)
    ax.set_ylim(-90, 90)
    ax.set_xlabel("Longitude (°)")
    ax.set_ylabel("Latitude (°)")
    ax.set_title("Inter-Event Time Regime Classification")
    ax.set_aspect("equal")

    return ax


# ---------------------------------------------------------------------------
# Oklahoma timeline
# ---------------------------------------------------------------------------

def plot_oklahoma_timeline(monthly_counts, rolling_b, phase_boundaries,
                           ax=None):
    """Dual-axis Oklahoma seismicity timeline.

    Parameters
    ----------
    monthly_counts : pd.DataFrame
        Columns: ``month`` (datetime), ``count``.
    rolling_b : pd.DataFrame
        Columns: ``center_time`` (datetime), ``b_value``.
    phase_boundaries : list of (datetime, str)
        List of (date, label) pairs for phase boundary lines.
    ax : matplotlib.axes.Axes or None

    Returns
    -------
    matplotlib.axes.Axes
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=DEFAULT_FIGSIZE)

    # Bars for event count
    ax.bar(monthly_counts["month"], monthly_counts["count"],
           width=25, color="#E63946", alpha=0.6, label="Monthly events")
    ax.set_xlabel("Date")
    ax.set_ylabel("Monthly event count", color="#E63946")
    ax.tick_params(axis="y", labelcolor="#E63946")

    # Secondary axis for b-value
    ax2 = ax.twinx()
    ax2.plot(rolling_b["center_time"], rolling_b["b_value"],
             color="#457B9D", linewidth=1.5, label="b-value")
    ax2.set_ylabel("b-value", color="#457B9D")
    ax2.tick_params(axis="y", labelcolor="#457B9D")

    # Phase boundaries
    for date, label in phase_boundaries:
        ax.axvline(date, color="black", linestyle="--", linewidth=1,
                   alpha=0.5)
        ax.text(date, ax.get_ylim()[1] * 0.95, f" {label}",
                fontsize=8, ha="left", va="top", style="italic")

    ax.set_title("Oklahoma Seismicity Timeline")

    return ax


# ---------------------------------------------------------------------------
# Entropy time series
# ---------------------------------------------------------------------------

def plot_entropy_timeseries(entropy_df, large_events=None, anomalies=None,
                            title="Shannon Entropy of Magnitude Distribution",
                            ax=None):
    """Plot entropy time series with optional event markers and anomalies.

    Parameters
    ----------
    entropy_df : pd.DataFrame
        Columns: ``center_time``, ``H``.
    large_events : array-like of datetime or None
        Times of large earthquakes to mark as vertical lines.
    anomalies : pd.DataFrame or None
        Entropy anomaly windows to shade.
    title : str
        Figure title.
    ax : matplotlib.axes.Axes or None

    Returns
    -------
    matplotlib.axes.Axes
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=DEFAULT_FIGSIZE)

    ax.plot(entropy_df["center_time"], entropy_df["H"],
            color="#2A9D8F", linewidth=1, alpha=0.8)

    if large_events is not None:
        for t in large_events:
            ax.axvline(t, color="#E63946", linewidth=0.8, alpha=0.5)

    if anomalies is not None and len(anomalies) > 0:
        ymin, ymax = ax.get_ylim()
        for _, row in anomalies.iterrows():
            ax.axvspan(row["center_time"],
                       row["center_time"] + np.timedelta64(7, "D"),
                       alpha=0.15, color="#F4A261")

    ax.set_xlabel("Date")
    ax.set_ylabel("Shannon Entropy H (bits)")
    ax.set_title(title)

    return ax


def plot_superposed_epoch(aligned_curves, null_mean=None, null_ci=None,
                          ax=None):
    """Superposed epoch analysis: average entropy around large events.

    Parameters
    ----------
    aligned_curves : np.ndarray
        Shape (n_events, n_timesteps).  Each row is an entropy curve
        aligned at t=0 (the event time).
    null_mean : np.ndarray or None
        Mean of null-model curves (same length as n_timesteps).
    null_ci : tuple of (np.ndarray, np.ndarray) or None
        (lower, upper) 95% confidence band from null model.
    ax : matplotlib.axes.Axes or None

    Returns
    -------
    matplotlib.axes.Axes
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=DEFAULT_FIGSIZE)

    n_steps = aligned_curves.shape[1]
    t_axis = np.linspace(-365, 365, n_steps)

    mean_curve = np.nanmean(aligned_curves, axis=0)
    ax.plot(t_axis, mean_curve, color="#2A9D8F", linewidth=2,
            label="Observed mean")

    if null_mean is not None:
        ax.plot(t_axis, null_mean, color="gray", linewidth=1.5,
                linestyle="--", label="Null model mean")

    if null_ci is not None:
        ax.fill_between(t_axis, null_ci[0], null_ci[1],
                        alpha=0.15, color="gray", label="Null 95% CI")

    ax.axvline(0, color="#E63946", linewidth=1.5, linestyle="-",
               label="Event (t=0)")
    ax.set_xlabel("Days relative to M7+ event")
    ax.set_ylabel("Shannon Entropy H (bits)")
    ax.set_title("Superposed Epoch Analysis")
    ax.legend()

    return ax

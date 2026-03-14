"""
Consistent visualization utilities for the seismic fingerprints project.

Provides region-aware color schemes, styled figures, and specialized plots
for b-value analysis, phase portraits, magnitude-frequency distributions,
inter-event statistics, cascade trees, and early-warning timelines.
"""

import warnings
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import networkx as nx

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

REGION_COLORS = {
    "oklahoma": "#E63946",
    "permian": "#F4A261",
    "socal": "#457B9D",
    "global": "#2A9D8F",
}

REGION_LABELS = {
    "oklahoma": "Oklahoma",
    "permian": "Permian Basin",
    "socal": "Southern California",
    "global": "Global",
}

DEFAULT_FIGSIZE = (12, 7)
DEFAULT_DPI = 300


# ---------------------------------------------------------------------------
# Style helpers
# ---------------------------------------------------------------------------

def setup_style():
    """Apply a consistent matplotlib style for the project.

    Tries ``seaborn-v0_8-whitegrid`` first (matplotlib >= 3.6), then falls
    back to ``seaborn-whitegrid`` (older matplotlib), and finally uses the
    default style.  Additionally sets default figure size, DPI, and font
    sizes for titles, axes labels, and tick labels.
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
# Plot functions
# ---------------------------------------------------------------------------

def plot_bvalue_timeseries(rolling_df, region, changepoints=None, ax=None):
    """Plot b-value over time with optional changepoint markers.

    Parameters
    ----------
    rolling_df : pandas.DataFrame
        Must contain columns ``time`` (or a datetime index) and ``b_value``.
    region : str
        Key into ``REGION_COLORS`` / ``REGION_LABELS``.
    changepoints : list of datetime-like, optional
        Times at which to draw vertical dashed lines.
    ax : matplotlib.axes.Axes, optional
        Axes to draw on.  A new figure is created when *None*.

    Returns
    -------
    matplotlib.axes.Axes
    """
    if ax is None:
        _, ax = plt.subplots(figsize=DEFAULT_FIGSIZE)

    color = REGION_COLORS.get(region, "#333333")
    label = REGION_LABELS.get(region, region)

    time = rolling_df.index if "time" not in rolling_df.columns else rolling_df["time"]
    ax.plot(time, rolling_df["b_value"], color=color, linewidth=1.5, label=label)

    if changepoints is not None:
        for cp in changepoints:
            ax.axvline(cp, color="grey", linestyle="--", alpha=0.7, linewidth=1)

    ax.set_xlabel("Time")
    ax.set_ylabel("b-value")
    ax.set_title(f"b-value Time Series — {label}")
    ax.legend()
    return ax


def plot_bvalue_phase_portrait(rolling_dfs, regions, ax=None):
    """Plot b-value vs log10(rate) phase portrait for one or more regions.

    This is the signature visualization of the project: each region traces
    a parametric curve in (log10 rate, b-value) space, color-coded by time
    using the *viridis* colormap, with arrows indicating the direction of
    temporal evolution.

    Parameters
    ----------
    rolling_dfs : list of pandas.DataFrame
        Each DataFrame must contain ``b_value`` and ``rate`` columns (and a
        datetime index or ``time`` column for ordering).
    regions : list of str
        Region keys corresponding to each DataFrame.
    ax : matplotlib.axes.Axes, optional

    Returns
    -------
    matplotlib.axes.Axes
    """
    if ax is None:
        _, ax = plt.subplots(figsize=DEFAULT_FIGSIZE)

    cmap = cm.viridis

    for df, region in zip(rolling_dfs, regions):
        b = df["b_value"].values
        rate = df["rate"].values
        log_rate = np.log10(rate + 1e-10)  # avoid log(0)
        n = len(b)
        if n < 2:
            continue

        # Normalise time index to [0, 1] for colormap
        t_norm = np.linspace(0, 1, n)

        # Draw coloured line segments
        for i in range(n - 1):
            ax.plot(
                log_rate[i : i + 2],
                b[i : i + 2],
                color=cmap(t_norm[i]),
                linewidth=1.5,
                alpha=0.8,
            )

        # Add arrows at regular intervals to show direction
        arrow_step = max(1, n // 8)
        for i in range(0, n - 1, arrow_step):
            dx = log_rate[i + 1] - log_rate[i]
            dy = b[i + 1] - b[i]
            ax.annotate(
                "",
                xy=(log_rate[i + 1], b[i + 1]),
                xytext=(log_rate[i], b[i]),
                arrowprops=dict(
                    arrowstyle="->",
                    color=cmap(t_norm[i]),
                    lw=1.5,
                ),
            )

        # Region label at the start of the trajectory
        label = REGION_LABELS.get(region, region)
        ax.annotate(
            label,
            xy=(log_rate[0], b[0]),
            fontsize=9,
            fontweight="bold",
            color=REGION_COLORS.get(region, "#333333"),
        )

    # Colorbar for time
    sm = cm.ScalarMappable(cmap=cmap, norm=matplotlib.colors.Normalize(0, 1))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, pad=0.02)
    cbar.set_label("Relative Time (early → late)")

    ax.set_xlabel("log₁₀(Seismicity Rate)")
    ax.set_ylabel("b-value")
    ax.set_title("b-value Phase Portrait")
    return ax


def plot_magnitude_frequency(magnitudes, mc, region, ax=None):
    """Gutenberg–Richter magnitude–frequency plot.

    Plots the cumulative number of events with magnitude >= M on a
    log-scaled y-axis, marks the completeness magnitude *mc*, and overlays
    the GR fit line above *mc*.

    Parameters
    ----------
    magnitudes : array-like
        Event magnitudes.
    mc : float
        Completeness magnitude.
    region : str
        Region key for colours and labels.
    ax : matplotlib.axes.Axes, optional

    Returns
    -------
    matplotlib.axes.Axes
    """
    if ax is None:
        _, ax = plt.subplots(figsize=DEFAULT_FIGSIZE)

    magnitudes = np.asarray(magnitudes)
    color = REGION_COLORS.get(region, "#333333")
    label = REGION_LABELS.get(region, region)

    # Empirical cumulative distribution (exceedance)
    m_sorted = np.sort(magnitudes)
    n_exceed = np.arange(len(m_sorted), 0, -1)
    ax.scatter(m_sorted, n_exceed, s=8, color=color, alpha=0.6, label=label)
    ax.set_yscale("log")

    # GR fit above Mc
    above = magnitudes[magnitudes >= mc]
    if len(above) > 1:
        b_val = np.log10(np.e) / (np.mean(above) - mc)
        a_val = np.log10(len(above)) + b_val * mc
        m_fit = np.linspace(mc, magnitudes.max(), 50)
        log_n_fit = a_val - b_val * m_fit
        ax.plot(m_fit, 10 ** log_n_fit, "--", color=color, linewidth=2,
                label=f"GR fit (b={b_val:.2f})")

    # Mark Mc
    ax.axvline(mc, color="k", linestyle=":", linewidth=1, label=f"Mc = {mc:.1f}")

    ax.set_xlabel("Magnitude")
    ax.set_ylabel("N (≥ M)")
    ax.set_title(f"Magnitude–Frequency Distribution — {label}")
    ax.legend()
    return ax


def plot_interevent_panel(results_dict, ax=None):
    """Panel comparing inter-event time statistics across regions.

    Displays coefficient of variation (CV), Hurst exponent (H), and
    best-fit distribution name for each region as a grouped bar / text
    summary.

    Parameters
    ----------
    results_dict : dict
        Maps region name (str) to a dict with keys ``cv`` (float),
        ``H`` (float), and ``best_dist`` (str).
    ax : matplotlib.axes.Axes, optional

    Returns
    -------
    matplotlib.axes.Axes
    """
    if ax is None:
        _, ax = plt.subplots(figsize=DEFAULT_FIGSIZE)

    regions = list(results_dict.keys())
    n = len(regions)
    x = np.arange(n)
    width = 0.35

    cvs = [results_dict[r]["cv"] for r in regions]
    hs = [results_dict[r]["H"] for r in regions]
    dists = [results_dict[r]["best_dist"] for r in regions]
    colors = [REGION_COLORS.get(r, "#333333") for r in regions]
    labels = [REGION_LABELS.get(r, r) for r in regions]

    bars_cv = ax.bar(x - width / 2, cvs, width, label="CV", color=colors, alpha=0.8)
    bars_h = ax.bar(x + width / 2, hs, width, label="Hurst H", color=colors, alpha=0.5)

    # Annotate best-fit distribution names above bars
    for i, dist_name in enumerate(dists):
        ymax = max(cvs[i], hs[i])
        ax.text(x[i], ymax + 0.05, dist_name, ha="center", fontsize=9,
                fontstyle="italic")

    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Value")
    ax.set_title("Inter-event Time Statistics by Region")
    ax.legend()
    return ax


def plot_cascade_tree(cascade_df, catalog_df, ax=None):
    """Visualise a single aftershock cascade as a directed graph.

    Node size is proportional to event magnitude; edges are coloured by
    the log10 of the time delay between parent and child events.

    Parameters
    ----------
    cascade_df : pandas.DataFrame
        Must contain ``parent_id`` and ``child_id`` columns identifying
        event pairs, plus a ``delay`` column (time delay in seconds or
        days).
    catalog_df : pandas.DataFrame
        Full catalogue with an ``id`` column and a ``magnitude`` column.
    ax : matplotlib.axes.Axes, optional

    Returns
    -------
    matplotlib.axes.Axes
    """
    if ax is None:
        _, ax = plt.subplots(figsize=DEFAULT_FIGSIZE)

    G = nx.DiGraph()

    # Build magnitude lookup
    mag_lookup = {}
    if "id" in catalog_df.columns and "magnitude" in catalog_df.columns:
        mag_lookup = dict(zip(catalog_df["id"], catalog_df["magnitude"]))

    # Add edges
    delays = []
    for _, row in cascade_df.iterrows():
        parent = row["parent_id"]
        child = row["child_id"]
        delay = row.get("delay", 1.0)
        G.add_edge(parent, child, delay=delay)
        delays.append(delay)

    if len(G) == 0:
        ax.text(0.5, 0.5, "Empty cascade", ha="center", va="center",
                transform=ax.transAxes)
        return ax

    # Node sizes proportional to magnitude
    node_sizes = []
    for node in G.nodes():
        mag = mag_lookup.get(node, 2.0)
        node_sizes.append(30 * (10 ** (mag - 1)))  # scale for visibility

    # Layout
    try:
        pos = nx.nx_agraph.graphviz_layout(G, prog="dot")
    except Exception:
        pos = nx.spring_layout(G, seed=42)

    # Edge colours by log10(delay)
    edge_delays = [G[u][v]["delay"] for u, v in G.edges()]
    log_delays = np.log10(np.array(edge_delays, dtype=float) + 1e-10)
    norm = matplotlib.colors.Normalize(vmin=log_delays.min(), vmax=log_delays.max())
    edge_colors = [cm.plasma(norm(ld)) for ld in log_delays]

    nx.draw_networkx_nodes(G, pos, ax=ax, node_size=node_sizes,
                           node_color="#457B9D", alpha=0.8)
    nx.draw_networkx_edges(G, pos, ax=ax, edge_color=edge_colors,
                           width=1.5, arrows=True, arrowsize=12)

    # Colorbar for delay
    sm = cm.ScalarMappable(cmap=cm.plasma, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, pad=0.02)
    cbar.set_label("log₁₀(Time Delay)")

    ax.set_title("Aftershock Cascade Tree")
    ax.axis("off")
    return ax


def plot_early_warning_timeline(permian_probs, oklahoma_probs,
                                years_permian, years_oklahoma, ax=None):
    """Timeline of induced-seismicity probability for two regions.

    Parameters
    ----------
    permian_probs : array-like
        Probability of induced origin for the Permian Basin over time.
    oklahoma_probs : array-like
        Probability of induced origin for Oklahoma over time.
    years_permian : array-like
        Year (or datetime) values corresponding to *permian_probs*.
    years_oklahoma : array-like
        Year (or datetime) values corresponding to *oklahoma_probs*.
    ax : matplotlib.axes.Axes, optional

    Returns
    -------
    matplotlib.axes.Axes
    """
    if ax is None:
        _, ax = plt.subplots(figsize=DEFAULT_FIGSIZE)

    ax.plot(years_permian, permian_probs,
            color=REGION_COLORS["permian"], linewidth=2,
            label=REGION_LABELS["permian"], marker="o", markersize=4)
    ax.plot(years_oklahoma, oklahoma_probs,
            color=REGION_COLORS["oklahoma"], linewidth=2,
            label=REGION_LABELS["oklahoma"], marker="s", markersize=4)

    # Reference threshold line
    ax.axhline(0.5, color="grey", linestyle="--", linewidth=1, alpha=0.6,
               label="p = 0.5 threshold")

    ax.fill_between(years_permian, permian_probs, alpha=0.15,
                     color=REGION_COLORS["permian"])
    ax.fill_between(years_oklahoma, oklahoma_probs, alpha=0.15,
                     color=REGION_COLORS["oklahoma"])

    ax.set_ylim(-0.05, 1.05)
    ax.set_xlabel("Year")
    ax.set_ylabel("P(induced)")
    ax.set_title("Early-Warning Timeline: Induced Seismicity Probability")
    ax.legend(loc="upper left")
    return ax

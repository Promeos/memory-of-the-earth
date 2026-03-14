---
name: fix-notebook
description: Debug and fix a failing Jupyter notebook by reading error output, diagnosing issues, and editing cells.
---

# Fix Notebook Agent

You are a debugging agent for Jupyter notebooks in the "Memory of the Earth" seismology project.

## Instructions

1. Read the error message or traceback provided
2. Read the relevant notebook cell(s) and source module(s)
3. Diagnose the root cause (API mismatch, data issue, missing guard, etc.)
4. Fix the notebook cell code using the Edit or NotebookEdit tool
5. Verify the fix makes sense given the src module APIs

## Available src modules

- `src.catalog` — load_raw_csvs, clean_catalog, deduplicate, estimate_mc, identify_temporal_gaps
- `src.fetch` — fetch_month, fetch_region, fetch_all_regions
- `src.gutenberg_richter` — estimate_b_value, bootstrap_b_value, rolling_b_value, segment_cost, detect_changepoints, ks_test_magnitudes, compute_bvalue_grid, compute_cv_b
- `src.interevent` — compute_interevent_times, fit_distributions, classify_regime, rolling_cv
- `src.entropy` — shannon_entropy, rolling_entropy, detect_anomalies, null_model_test
- `src.omori` — fit_omori, compute_recovery_fraction, fit_recovery_exponential
- `src.spatial` — haversine, build_grid, assign_cells, events_in_radius
- `src.plotting` — setup_style, save_figure, plot_global_map, plot_bvalue_volatility_map, plot_recovery_gallery, plot_regime_map, plot_oklahoma_timeline, plot_entropy_timeseries, plot_superposed_epoch

## Common fixes

- Add `matplotlib.use('Agg')` before importing pyplot for headless execution
- Add length/NaN checks before operations that assume data availability
- Fix column name mismatches between notebook code and actual CSV columns
- Wrap long-running null model tests with reduced iteration counts

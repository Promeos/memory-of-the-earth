# Seismic Fingerprints - Project Guide

## Overview
Data science project detecting statistical signatures of human-caused earthquakes in 25 years of USGS data. Compares Oklahoma (induced seismicity ground truth), Southern California (tectonic control), and the Permian Basin (prospective test region).

## Tech Stack
- Python 3.10+
- Core: pandas, numpy, scipy, scikit-learn, matplotlib, seaborn, networkx
- Data: USGS ComCat API (CSV format, paginated monthly)
- Storage: Raw CSVs in data/raw/, processed Parquet in data/processed/

## Project Structure
- `src/` — Shared library modules (fetch, catalog, analysis, plotting)
- `notebooks/` — 5 Jupyter notebooks (00-04), each self-contained
- `data/raw/` — Raw monthly CSV pulls from USGS API
- `data/processed/` — Cleaned Parquet files per region
- `figures/` — Publication-quality PNG outputs (300 DPI)

## Conventions

### Region Colors (use consistently everywhere)
- Oklahoma: `#E63946` (red)
- Permian Basin: `#F4A261` (amber/orange)
- Southern California: `#457B9D` (steel blue)
- Global/other: `#2A9D8F` (teal)

### Visualization Standards
- 300 DPI, default `figsize=(12, 7)`
- Use `plt.style.use('seaborn-v0_8-whitegrid')`
- Always label axes with units
- Save figures to both `figures/` and display inline

### Data Quality Rules
- Always estimate Mc (magnitude of completeness) per region-year using MAXC method
- Only analyze events above Mc
- Filter to `type == "earthquake"` (exclude quarry blasts, explosions)
- Deduplicate: group by (time ± 2s, lat ± 0.05°, lon ± 0.05°)
- Use conservative common Mc for cross-region comparisons (M2.5 or M3.0)

### Regional Bounding Boxes
- Oklahoma: lat 33.5–37.5°N, lon 100.0–94.5°W, min mag 1.0
- Permian Basin: lat 30.5–33.5°N, lon 105.0–100.5°W, min mag 1.0
- Southern California: lat 32.0–36.5°N, lon 121.0–114.5°W, min mag 1.0
- Global: no bounds, min mag 2.5

## Running
```bash
pip install -r requirements.txt
# Notebooks should be run in order: 00 → 01 → 02 → 03 → 04
# Notebook 00 fetches data from USGS API (takes ~30-60 min)
```

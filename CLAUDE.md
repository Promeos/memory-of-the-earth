# The Memory of the Earth - Project Guide

## Overview
Data science project quantifying temporal structure and stress dynamics in 25 years of global seismicity using the USGS ANSS Comprehensive Earthquake Catalog (ComCat), 2000–2025.

## Tech Stack
- Python 3.10+
- Core: pandas, numpy, scipy, scikit-learn, matplotlib, seaborn
- Data: USGS ComCat API (CSV format, paginated monthly)
- Storage: Raw CSVs in data/raw/, analysis-ready CSVs in data/

## Project Structure
- `src/` — 9 shared library modules (fetch, catalog, gutenberg_richter, interevent, entropy, omori, spatial, plotting, external_data)
- `00_data_acquisition.py` — Standalone data fetch/clean script (run first)
- `01–05_*.ipynb` — 5 analysis notebooks in project root, each self-contained
- `data/raw/` — Raw monthly CSV pulls from USGS API
- `data/` — Cleaned catalog CSVs (earthquake_catalog_global.csv, earthquake_catalog_oklahoma.csv)
- `figures/` — Publication-quality PNG outputs (300 DPI)

## Conventions

### Color Palette
- Oklahoma / alerts: `#E63946` (red)
- Teal / global: `#2A9D8F`
- Steel blue: `#457B9D`
- Amber: `#F4A261`
- Gray: `#AAAAAA`

### Visualization Standards
- 300 DPI, default `figsize=(12, 7)`
- Use `seaborn-v0_8-whitegrid` style (with fallback)
- Always label axes with units
- Save figures to `figures/` and display inline

### Data Quality Rules
- Always estimate Mc per region using MAXC method (+0.2 correction, Wiemer & Wyss 2000)
- Only analyze events above Mc
- Filter to `type == "earthquake"` (exclude quarry blasts, explosions)
- Deduplicate: group by (time ± 2s, lat ± 0.05°, lon ± 0.05°)

### Regional Bounding Boxes
- Oklahoma: lat 33.5–37.5°N, lon 100.0–94.5°W, min mag 1.0
- Global: no bounds, min mag 2.5

## Running
```bash
pip install -r requirements.txt
python 00_data_acquisition.py  # Fetches/cleans data (~30-60 min if no cache)
jupyter notebook 01_bvalue_stability_atlas.ipynb  # Then run notebooks 01-05
```

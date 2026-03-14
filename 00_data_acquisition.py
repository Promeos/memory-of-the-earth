#!/usr/bin/env python3
"""Data Acquisition and Cleaning — The Memory of the Earth.

Fetches earthquake catalogs from USGS ComCat, cleans, deduplicates, and
validates them, producing analysis-ready CSVs for the global (M2.5+) and
Oklahoma (M1.0+) catalogs.

If raw monthly CSVs already exist in data/raw/{region}/, they are loaded
directly instead of re-fetching from the API (~30-60 min saved).

Output
------
data/earthquake_catalog_global.csv
data/earthquake_catalog_oklahoma.csv
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Add project root to path so src is importable
sys.path.insert(0, str(Path(__file__).resolve().parent))

from src.fetch import fetch_all_regions, REGIONS
from src.catalog import load_raw_csvs, clean_catalog, deduplicate, estimate_mc, identify_temporal_gaps

DATA_DIR = Path("data")
RAW_DIR = DATA_DIR / "raw"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _print_summary(df, name):
    """Print summary statistics for a catalog."""
    print(f"\n{'='*60}")
    print(f"Summary: {name}")
    print(f"{'='*60}")
    print(f"  Total events:     {len(df):,}")
    print(f"  Date range:       {df['time'].min()} → {df['time'].max()}")
    print(f"  Magnitude range:  {df['mag'].min():.1f} → {df['mag'].max():.1f}")
    print(f"  Mean magnitude:   {df['mag'].mean():.2f}")

    # Events per year
    df_year = df.copy()
    df_year["year"] = df_year["time"].dt.year
    yearly = df_year.groupby("year").size()
    print(f"  Events/year:      {yearly.mean():.0f} avg "
          f"(min {yearly.min()}, max {yearly.max()})")

    # Mc estimate
    mc = estimate_mc(df["mag"].values)
    print(f"  Estimated Mc:     {mc:.1f}")

    # Magnitude types
    if "magType" in df.columns:
        top_types = df["magType"].value_counts().head(5)
        print(f"  Top mag types:    {dict(top_types)}")

    print()


def _validate_catalog(df, name):
    """Run data validation checks and print results."""
    print(f"\nValidation: {name}")
    print("-" * 40)

    # Temporal gaps
    gaps = identify_temporal_gaps(df, threshold_hours=24)
    n_gaps = len(gaps)
    if n_gaps > 0:
        print(f"  Temporal gaps (>24h): {n_gaps}")
        if n_gaps <= 5:
            for _, row in gaps.iterrows():
                print(f"    {row['gap_start']} → {row['gap_end']} "
                      f"({row['gap_hours']:.1f} h)")
        else:
            max_gap = gaps.loc[gaps["gap_hours"].idxmax()]
            print(f"    Longest: {max_gap['gap_start']} → "
                  f"{max_gap['gap_end']} ({max_gap['gap_hours']:.1f} h)")
    else:
        print("  Temporal gaps (>24h): none")

    # Magnitude type mixing
    if "magType" in df.columns:
        n_types = df["magType"].nunique()
        print(f"  Distinct magnitude types: {n_types}")
        if n_types > 1:
            mw_frac = (df["magType"].str.lower().str.startswith("mw")).mean()
            print(f"  Fraction moment magnitude (mw*): {mw_frac:.1%}")

    # Duplicates check (should be near-zero after dedup)
    print(f"  NaN magnitudes: {df['mag'].isna().sum()}")
    print()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    # ------------------------------------------------------------------
    # Step 1: Load or fetch global catalog
    # ------------------------------------------------------------------
    global_raw_dirs = [RAW_DIR / "global", RAW_DIR / "global_M2.5"]
    global_df = None

    for d in global_raw_dirs:
        if d.is_dir() and list(d.glob("*.csv")):
            print(f"[load] Found existing raw CSVs in {d}")
            global_df = load_raw_csvs(d)
            break

    if global_df is None:
        print("[fetch] No existing global data found — fetching from USGS API")
        cfg = REGIONS["global"]
        from src.fetch import fetch_region
        global_df = fetch_region(
            name="global",
            years=cfg["years"],
            min_mag=cfg["min_mag"],
            lat_range=cfg["lat_range"],
            lon_range=cfg["lon_range"],
            output_dir=RAW_DIR / "global",
        )

    print(f"[raw]  Global raw events: {len(global_df):,}")

    # ------------------------------------------------------------------
    # Step 2: Clean and deduplicate global catalog
    # ------------------------------------------------------------------
    global_df = clean_catalog(global_df)
    print(f"[clean] After cleaning: {len(global_df):,}")

    global_df = deduplicate(global_df)
    print(f"[dedup] After deduplication: {len(global_df):,}")

    # Save
    global_out = DATA_DIR / "earthquake_catalog_global.csv"
    global_df.to_csv(global_out, index=False)
    print(f"[save] Global catalog → {global_out}")

    _print_summary(global_df, "Global M2.5+")
    _validate_catalog(global_df, "Global M2.5+")

    # ------------------------------------------------------------------
    # Step 3: Load or fetch Oklahoma catalog
    # ------------------------------------------------------------------
    ok_raw_dirs = [RAW_DIR / "oklahoma", RAW_DIR / "oklahoma_M1.0"]
    ok_df = None

    for d in ok_raw_dirs:
        if d.is_dir() and list(d.glob("*.csv")):
            print(f"[load] Found existing raw CSVs in {d}")
            ok_df = load_raw_csvs(d)
            break

    if ok_df is None:
        print("[fetch] No existing Oklahoma data found — fetching from USGS API")
        cfg = REGIONS["oklahoma"]
        from src.fetch import fetch_region
        ok_df = fetch_region(
            name="oklahoma",
            years=cfg["years"],
            min_mag=cfg["min_mag"],
            lat_range=cfg["lat_range"],
            lon_range=cfg["lon_range"],
            output_dir=RAW_DIR / "oklahoma",
        )

    print(f"[raw]  Oklahoma raw events: {len(ok_df):,}")

    # ------------------------------------------------------------------
    # Step 4: Clean and deduplicate Oklahoma catalog
    # ------------------------------------------------------------------
    ok_df = clean_catalog(ok_df)
    print(f"[clean] After cleaning: {len(ok_df):,}")

    ok_df = deduplicate(ok_df)
    print(f"[dedup] After deduplication: {len(ok_df):,}")

    # Save
    ok_out = DATA_DIR / "earthquake_catalog_oklahoma.csv"
    ok_df.to_csv(ok_out, index=False)
    print(f"[save] Oklahoma catalog → {ok_out}")

    _print_summary(ok_df, "Oklahoma M1.0+")
    _validate_catalog(ok_df, "Oklahoma M1.0+")

    # ------------------------------------------------------------------
    # Done
    # ------------------------------------------------------------------
    print("=" * 60)
    print("Data acquisition complete.")
    print(f"  Global:   {global_out}  ({len(global_df):,} events)")
    print(f"  Oklahoma: {ok_out}  ({len(ok_df):,} events)")
    print("=" * 60)


if __name__ == "__main__":
    main()

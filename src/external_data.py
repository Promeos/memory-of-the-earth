"""
External dataset parsers for supplementary data sources.

Handles: Global CMT focal mechanisms, GSRM strain rates,
OCC injection well volumes, SCEDC regional catalog,
IHFC Global Heat Flow Database.
"""

import os
import re
import numpy as np
import pandas as pd
from pathlib import Path

EXTERNAL_DIR = Path(__file__).parent.parent / "data" / "external"


# ---------------------------------------------------------------------------
# 1.  Global CMT NDK parser
# ---------------------------------------------------------------------------

def parse_ndk(filepath: str | Path) -> pd.DataFrame:
    """Parse Global CMT NDK format into a DataFrame with focal mechanism info.

    Each event is a 5-line block. We extract origin time, centroid location,
    Mw, moment tensor components, and both nodal planes (strike/dip/rake).
    """
    filepath = Path(filepath)
    lines = filepath.read_text().strip().split("\n")
    n_events = len(lines) // 5
    records = []

    for i in range(n_events):
        block = lines[i * 5 : i * 5 + 5]
        if len(block) < 5:
            continue

        # Line 1: Hypocenter reference
        l1 = block[0]
        # Line 3: Centroid parameters
        l3 = block[2]
        # Line 4: Moment tensor components
        l4 = block[3]
        # Line 5: Nodal planes & principal axes
        l5 = block[4]

        try:
            # Parse origin time from line 1
            parts1 = l1.split()
            date_str = parts1[1]  # YYYY/MM/DD
            time_str = parts1[2]  # HH:MM:SS.s
            lat_hypo = float(parts1[3])
            lon_hypo = float(parts1[4])
            depth_hypo = float(parts1[5])

            # Parse centroid from line 3
            parts3 = l3.split()
            lat_cmt = float(parts3[3])
            lon_cmt = float(parts3[5])
            depth_cmt = float(parts3[7])

            # Parse moment tensor exponent and scalar moment from line 4
            parts4 = l4.split()
            exponent = int(parts4[0])
            mrr = float(parts4[1]) * 10**exponent
            mtt = float(parts4[3]) * 10**exponent
            mpp = float(parts4[5]) * 10**exponent
            mrt = float(parts4[7]) * 10**exponent
            mrp = float(parts4[9]) * 10**exponent
            mtp = float(parts4[11]) * 10**exponent

            # Scalar moment M0 = max eigenvalue of moment tensor
            M0 = max(abs(mrr), abs(mtt), abs(mpp)) * 1.0
            # Mw from scalar moment
            Mw = (2 / 3) * (np.log10(abs(mrr + mtt + mpp) * 1e7 + 1e10) - 16.1)

            # Parse nodal planes from line 5
            parts5 = l5.split()
            # Format: V10 val1 plunge1 azimuth1 val2 p2 a2 val3 p3 a3 M0 str1 dip1 rake1 str2 dip2 rake2
            scalar_moment = float(parts5[10]) * 10**exponent
            strike1 = float(parts5[11])
            dip1 = float(parts5[12])
            rake1 = float(parts5[13])
            strike2 = float(parts5[14])
            dip2 = float(parts5[15])
            rake2 = float(parts5[16])
            if scalar_moment > 0:
                Mw_calc = (2 / 3) * (np.log10(scalar_moment * 1e-7) - 9.1)
            else:
                continue  # skip events with zero moment

            # Classify faulting style from rake1
            faulting = classify_faulting(rake1)

            # Parse datetime
            dt_str = f"{date_str} {time_str}"
            try:
                origin_time = pd.to_datetime(dt_str, format="%Y/%m/%d %H:%M:%S.%f")
            except ValueError:
                origin_time = pd.to_datetime(dt_str, format="%Y/%m/%d %H:%M:%S")
            origin_time = origin_time.tz_localize("UTC")

            records.append({
                "time": origin_time,
                "latitude": lat_cmt,
                "longitude": lon_cmt,
                "depth": depth_cmt,
                "Mw": round(Mw_calc, 2),
                "scalar_moment": scalar_moment,
                "strike1": strike1, "dip1": dip1, "rake1": rake1,
                "strike2": strike2, "dip2": dip2, "rake2": rake2,
                "faulting": faulting,
            })
        except (ValueError, IndexError):
            continue

    return pd.DataFrame(records)


def classify_faulting(rake: float) -> str:
    """Classify faulting style from rake angle."""
    rake = rake % 360
    if rake > 180:
        rake -= 360
    if -120 <= rake <= -60:
        return "normal"
    elif 60 <= rake <= 120:
        return "thrust"
    elif abs(rake) <= 30 or abs(rake) >= 150:
        return "strike-slip"
    else:
        return "oblique"


def load_gcmt_catalog(min_year: int = 2000) -> pd.DataFrame:
    """Load parsed Global CMT catalog, filtering to min_year+."""
    ndk_path = EXTERNAL_DIR / "gcmt_combined.ndk"
    if not ndk_path.exists():
        # Try individual files
        ndk_path = EXTERNAL_DIR / "gcmt_jan76_dec20.ndk"
    if not ndk_path.exists():
        raise FileNotFoundError(f"No CMT NDK file found in {EXTERNAL_DIR}")

    df = parse_ndk(ndk_path)
    df = df[df["time"].dt.year >= min_year].reset_index(drop=True)
    return df


# ---------------------------------------------------------------------------
# 2.  GSRM strain rate grid parser
# ---------------------------------------------------------------------------

def load_gsrm_strain(filepath: str | Path = None) -> pd.DataFrame:
    """Load GSRM v2.1 strain rate grid.

    Returns DataFrame with columns: lon, lat, e1, e2, azimuth, dilatation, second_invariant
    """
    if filepath is None:
        filepath = EXTERNAL_DIR / "GSRM_average_strain_v2.1.txt"
    filepath = Path(filepath)

    # Skip comment lines starting with #
    rows = []
    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 5:
                try:
                    lat = float(parts[0])
                    lon = float(parts[1])
                    e1 = float(parts[2])   # max principal strain rate (nanostr/yr)
                    e2 = float(parts[3])   # min principal strain rate
                    azimuth = float(parts[4])  # azimuth of e1
                    dilatation = e1 + e2
                    second_invariant = np.sqrt(e1**2 + e2**2)
                    rows.append({
                        "lon": lon, "lat": lat,
                        "e1_nanostr_yr": e1, "e2_nanostr_yr": e2,
                        "azimuth_e1": azimuth,
                        "dilatation": dilatation,
                        "second_invariant": second_invariant,
                    })
                except ValueError:
                    continue
    return pd.DataFrame(rows)


def resample_strain_to_grid(strain_df: pd.DataFrame, grid_size: float = 2.0) -> pd.DataFrame:
    """Resample 0.1° strain rate grid to coarser grid (e.g. 2°×2°).

    Returns mean strain rate per grid cell.
    """
    strain_df = strain_df.copy()
    strain_df["lat_bin"] = np.floor(strain_df["lat"] / grid_size) * grid_size + grid_size / 2
    strain_df["lon_bin"] = np.floor(strain_df["lon"] / grid_size) * grid_size + grid_size / 2

    grouped = strain_df.groupby(["lat_bin", "lon_bin"]).agg(
        mean_second_invariant=("second_invariant", "mean"),
        mean_dilatation=("dilatation", "mean"),
        mean_e1=("e1_nanostr_yr", "mean"),
        mean_e2=("e2_nanostr_yr", "mean"),
    ).reset_index()
    return grouped


# ---------------------------------------------------------------------------
# 3.  OCC injection well volume parser
# ---------------------------------------------------------------------------

def load_occ_uic_volumes(years: list[int] = None) -> pd.DataFrame:
    """Load OCC UIC injection volumes from Excel files.

    Returns DataFrame with columns: API, lat, lon, year, month, volume_bbl, pressure_psi
    """
    if years is None:
        years = list(range(2006, 2025))

    all_records = []

    for year in years:
        if 2006 <= year <= 2010:
            path = EXTERNAL_DIR / "occ_uic_2006_2010.xls"
            if not path.exists():
                continue
            try:
                df = pd.read_excel(path, engine="xlrd")
                df.columns = [c.strip() for c in df.columns]
                df = df[df["YEAR"] == year]
                for _, row in df.iterrows():
                    total_vol = row.get("TOTAL VOLUME", 0)
                    if pd.isna(total_vol):
                        total_vol = 0
                    all_records.append({
                        "api": str(row.get("API_NUMBER", "")),
                        "latitude": row.get("DECIMAL_LAT", np.nan),
                        "longitude": row.get("DECIMAL_LONG", np.nan),
                        "year": int(year),
                        "month": 0,  # annual total
                        "volume_bbl": float(total_vol),
                        "pressure_psi": np.nan,
                    })
            except Exception:
                continue
        else:
            path = EXTERNAL_DIR / f"occ_uic_{year}.xlsx"
            if not path.exists():
                continue
            try:
                df = pd.read_excel(path, engine="openpyxl")
                months = ["Jan", "Feb", "Mar", "Apr", "May", "Jun",
                          "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
                for _, row in df.iterrows():
                    api = str(row.get("API", ""))
                    lat = row.get("LAT", np.nan)
                    lon = row.get("LON", np.nan)
                    for m_idx, month_name in enumerate(months, 1):
                        vol_col = f"{month_name} Vol"
                        psi_col = f"{month_name} PSI"
                        vol = row.get(vol_col, 0)
                        psi = row.get(psi_col, np.nan)
                        if pd.isna(vol):
                            vol = 0
                        all_records.append({
                            "api": api,
                            "latitude": lat,
                            "longitude": lon,
                            "year": int(year),
                            "month": m_idx,
                            "volume_bbl": float(vol),
                            "pressure_psi": float(psi) if not pd.isna(psi) else np.nan,
                        })
            except Exception:
                continue

    return pd.DataFrame(all_records)


def aggregate_oklahoma_injection(uic_df: pd.DataFrame) -> pd.DataFrame:
    """Aggregate injection volumes to monthly totals for Oklahoma.

    Returns DataFrame with columns: year, month, total_volume_bbl, n_wells
    """
    # Filter to wells with valid coordinates in Oklahoma bounding box
    mask = (
        uic_df["latitude"].between(33.5, 37.5) &
        uic_df["longitude"].between(-100.0, -94.5) &
        (uic_df["volume_bbl"] > 0)
    )
    ok_df = uic_df[mask].copy()

    # Split annual (month=0) and monthly data
    annual = ok_df[ok_df["month"] == 0]
    monthly = ok_df[ok_df["month"] > 0]

    parts = []
    if len(annual) > 0:
        ann_grouped = annual.groupby("year").agg(
            total_volume_bbl=("volume_bbl", "sum"),
            n_wells=("api", "nunique"),
        ).reset_index()
        ann_grouped["month"] = 6  # mid-year proxy
        parts.append(ann_grouped)

    if len(monthly) > 0:
        mon_grouped = monthly.groupby(["year", "month"]).agg(
            total_volume_bbl=("volume_bbl", "sum"),
            n_wells=("api", "nunique"),
        ).reset_index()
        parts.append(mon_grouped)

    if not parts:
        return pd.DataFrame(columns=["year", "month", "total_volume_bbl", "n_wells", "date"])

    grouped = pd.concat(parts, ignore_index=True)
    grouped["date"] = pd.to_datetime(
        grouped["year"].astype(str) + "-" + grouped["month"].astype(str).str.zfill(2) + "-15"
    )
    return grouped.sort_values("date").reset_index(drop=True)


# ---------------------------------------------------------------------------
# 4.  SCEDC regional catalog parser
# ---------------------------------------------------------------------------

def load_scedc_catalog(filepath: str | Path = None) -> pd.DataFrame:
    """Load SCEDC Southern California catalog from FDSN text format."""
    if filepath is None:
        filepath = EXTERNAL_DIR / "scedc_combined.txt"
    filepath = Path(filepath)

    rows = []
    with open(filepath) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = [p.strip() for p in line.split("|")]
            if len(parts) < 11:
                continue
            try:
                event_id = parts[0]
                time_str = parts[1].strip()
                lat = float(parts[2])
                lon = float(parts[3])
                depth = float(parts[4])
                mag = float(parts[10])
                mag_type = parts[9].strip()

                time = pd.to_datetime(time_str, utc=True)

                rows.append({
                    "event_id": event_id,
                    "time": time,
                    "latitude": lat,
                    "longitude": lon,
                    "depth": depth,
                    "mag": mag,
                    "magType": mag_type,
                })
            except (ValueError, IndexError):
                continue
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# 5.  USGS Oklahoma (all magnitudes) loader
# ---------------------------------------------------------------------------

def load_usgs_oklahoma_full(filepath: str | Path = None) -> pd.DataFrame:
    """Load USGS Oklahoma catalog with all magnitudes."""
    if filepath is None:
        filepath = EXTERNAL_DIR / "ogs_usgs_oklahoma_full.csv"
    filepath = Path(filepath)
    df = pd.read_csv(filepath)
    df["time"] = pd.to_datetime(df["time"], format="ISO8601", utc=True)
    df = df[df["type"] == "earthquake"].reset_index(drop=True)
    return df


# ---------------------------------------------------------------------------
# 6.  IHFC Global Heat Flow Database parser
# ---------------------------------------------------------------------------

def load_heat_flow(filepath: str | Path = None, min_quality: str = None) -> pd.DataFrame:
    """Load IHFC Global Heat Flow Database.

    Supports both 2024 format (tab-separated, coded headers with real names on row 13)
    and 2021/2023 format (semicolon-separated CSV with readable headers).

    Parameters
    ----------
    filepath : path to the data file. If None, tries 2024 then 2023 then 2021.
    min_quality : if set, filter to entries with quality score <= this value
                  (lower = better in the GHFDB scoring system).

    Returns
    -------
    DataFrame with columns: lat, lon, q_mw_m2 (heat flow), q_uncertainty,
    elevation, environment, q_method, quality_score
    """
    if filepath is None:
        # Try releases in order of preference (newest first)
        candidates = [
            EXTERNAL_DIR / "IHFC_2024_GHFDB_v.2026.03.txt",
            EXTERNAL_DIR / "IHFC_2023_GHFDB.csv",
            EXTERNAL_DIR / "IHFC_2021_GHFDB.csv",
        ]
        for c in candidates:
            if c.exists():
                filepath = c
                break
        if filepath is None:
            raise FileNotFoundError(f"No IHFC heat flow file found in {EXTERNAL_DIR}")
    filepath = Path(filepath)

    if filepath.suffix == ".txt" and "2024" in filepath.name:
        return _load_hf_2024(filepath, min_quality)
    else:
        return _load_hf_2021_2023(filepath, min_quality)


def _load_hf_2024(filepath: Path, min_quality: str = None) -> pd.DataFrame:
    """Parse the 2024 GHFDB release (tab-separated, coded P/C headers)."""
    # Row 13 (0-indexed line 12) has real column names; data starts at row 14
    with open(filepath, encoding="latin-1") as f:
        all_lines = f.readlines()

    # Find the header row with real column names (contains "q\t" at start)
    header_idx = None
    for i, line in enumerate(all_lines):
        stripped = line.strip()
        if stripped.startswith("q\tq_uncertainty\tname"):
            header_idx = i
            break

    if header_idx is None:
        raise ValueError("Could not find real column name row in 2024 heat flow file")

    col_names = all_lines[header_idx].strip().split("\t")
    data_lines = all_lines[header_idx + 1:]

    records = []
    for line in data_lines:
        parts = line.strip().split("\t")
        if len(parts) < len(col_names):
            parts.extend([""] * (len(col_names) - len(parts)))
        row = dict(zip(col_names, parts))

        try:
            q = float(row.get("q", "")) if row.get("q", "").strip() else np.nan
            lat = float(row.get("lat_NS", "")) if row.get("lat_NS", "").strip() else np.nan
            lon = float(row.get("long_EW", "")) if row.get("long_EW", "").strip() else np.nan
        except ValueError:
            continue

        if np.isnan(q) or np.isnan(lat) or np.isnan(lon):
            continue

        q_unc = _safe_float(row.get("q_uncertainty", ""))
        elev = _safe_float(row.get("elevation", ""))
        env = row.get("environment", "").strip()
        q_method = row.get("q_method", "").strip()
        quality = _safe_float(row.get("Quality_Score_Parent", ""))

        records.append({
            "lat": lat, "lon": lon,
            "q_mw_m2": q, "q_uncertainty": q_unc,
            "elevation": elev, "environment": env,
            "q_method": q_method, "quality_score": quality,
        })

    df = pd.DataFrame(records)
    if min_quality is not None and "quality_score" in df.columns:
        df = df[df["quality_score"] <= float(min_quality)].reset_index(drop=True)
    return df


def _load_hf_2021_2023(filepath: Path, min_quality: str = None) -> pd.DataFrame:
    """Parse the 2021 or 2023 GHFDB release (semicolon-separated CSV)."""
    df = pd.read_csv(filepath, sep=";", low_memory=False)

    # Map column names to our standard schema
    col_map = {
        "lat": "lat", "lng": "lon",
        "q": "q_mw_m2", "q_unc": "q_uncertainty",
        "elevation": "elevation", "env": "environment",
        "q_method": "q_method",
    }

    result = pd.DataFrame()
    for src_col, dst_col in col_map.items():
        if src_col in df.columns:
            result[dst_col] = df[src_col]

    # Quality score column varies by release
    for qcol in ["qc", "Quality_Score_Parent"]:
        if qcol in df.columns:
            result["quality_score"] = pd.to_numeric(df[qcol], errors="coerce")
            break

    # Convert numeric columns
    for col in ["lat", "lon", "q_mw_m2", "q_uncertainty", "elevation"]:
        if col in result.columns:
            result[col] = pd.to_numeric(result[col], errors="coerce")

    # Drop rows without core data
    result = result.dropna(subset=["lat", "lon", "q_mw_m2"]).reset_index(drop=True)

    if min_quality is not None and "quality_score" in result.columns:
        result = result[result["quality_score"] <= float(min_quality)].reset_index(drop=True)

    return result


def resample_heat_flow_to_grid(hf_df: pd.DataFrame, grid_size: float = 2.0) -> pd.DataFrame:
    """Resample point heat flow measurements to a grid by taking the median per cell.

    Uses median (robust to outliers) rather than mean.
    """
    hf_df = hf_df.copy()
    hf_df["lat_bin"] = np.floor(hf_df["lat"] / grid_size) * grid_size + grid_size / 2
    hf_df["lon_bin"] = np.floor(hf_df["lon"] / grid_size) * grid_size + grid_size / 2

    grouped = hf_df.groupby(["lat_bin", "lon_bin"]).agg(
        median_q=("q_mw_m2", "median"),
        mean_q=("q_mw_m2", "mean"),
        std_q=("q_mw_m2", "std"),
        n_measurements=("q_mw_m2", "count"),
    ).reset_index()
    return grouped


def _safe_float(s: str) -> float:
    """Convert string to float, returning NaN on failure."""
    try:
        return float(s.strip()) if s.strip() else np.nan
    except (ValueError, AttributeError):
        return np.nan


# ---------------------------------------------------------------------------
# 7.  PB2002 Plate Boundary Model (Bird, 2003)
# ---------------------------------------------------------------------------

# Full plate name lookup
PB2002_PLATE_NAMES = {
    "AF": "Africa", "AN": "Antarctica", "AP": "Altiplano", "AR": "Arabia",
    "AS": "Aegean Sea", "AT": "Anatolia", "AU": "Australia", "BH": "Birds Head",
    "BR": "Balmoral Reef", "BS": "Banda Sea", "BU": "Burma", "CA": "Caribbean",
    "CL": "Caroline", "CO": "Cocos", "CR": "Conway Reef", "EA": "Easter",
    "EU": "Eurasia", "FT": "Futuna", "GP": "Galapagos", "IN": "India",
    "JF": "Juan de Fuca", "JZ": "Juan Fernandez", "KE": "Kermadec",
    "MA": "Mariana", "MN": "Manus", "MO": "Maoke", "MS": "Molucca Sea",
    "NA": "North America", "NB": "North Bismarck", "ND": "North Andes",
    "NH": "New Hebrides", "NI": "Niuafo'ou", "NZ": "Nazca", "OK": "Okhotsk",
    "ON": "Okinawa", "PA": "Pacific", "PM": "Panama", "PS": "Philippine Sea",
    "RI": "Rivera", "SA": "South America", "SB": "South Bismarck",
    "SC": "Scotia", "SL": "Shetland", "SO": "Somalia", "SS": "Solomon Sea",
    "SU": "Sunda", "SW": "Sandwich", "TI": "Timor", "TO": "Tonga",
    "WL": "Woodlark", "YA": "Yangtze",
}

# Boundary type descriptions and tectonic classification
PB2002_BOUNDARY_TYPES = {
    "SUB": "subduction",
    "OCB": "oceanic convergent",
    "CCB": "continental convergent",
    "OTF": "oceanic transform",
    "CTF": "continental transform",
    "OSR": "oceanic spreading ridge",
    "CRB": "continental rift",
}

# Higher-level tectonic setting categories
PB2002_SETTING_MAP = {
    "SUB": "subduction",
    "OCB": "subduction",
    "CCB": "convergent",
    "OTF": "transform",
    "CTF": "transform",
    "OSR": "spreading",
    "CRB": "rift",
}


def load_pb2002_plates(filepath: str | Path = None) -> dict[str, np.ndarray]:
    """Load PB2002 plate polygons.

    Returns dict mapping plate code (e.g. 'NA') to Nx2 array of (lon, lat) vertices.
    """
    if filepath is None:
        filepath = EXTERNAL_DIR / "PB2002_plates.dig.txt"
    filepath = Path(filepath)

    plates = {}
    current_plate = None
    vertices = []

    with open(filepath) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("***"):
                # End of segment — but plates.dig doesn't use this marker
                continue
            # Check if this is a plate header (short alphanumeric code, no comma)
            if "," not in line and len(line) <= 5 and line.isalpha():
                if current_plate and vertices:
                    plates[current_plate] = np.array(vertices)
                current_plate = line
                vertices = []
            else:
                # Coordinate line: lon,lat
                parts = line.split(",")
                if len(parts) == 2:
                    try:
                        lon = float(parts[0])
                        lat = float(parts[1])
                        vertices.append([lon, lat])
                    except ValueError:
                        continue

    if current_plate and vertices:
        plates[current_plate] = np.array(vertices)

    return plates


def load_pb2002_boundaries(filepath: str | Path = None) -> pd.DataFrame:
    """Load PB2002 boundary segments with their types from the steps file.

    Returns DataFrame with columns: plate_pair, lon1, lat1, lon2, lat2,
    length_km, azimuth, boundary_type, tectonic_setting
    """
    if filepath is None:
        filepath = EXTERNAL_DIR / "PB2002_steps.dat.txt"
    filepath = Path(filepath)

    records = []
    with open(filepath) as f:
        for line in f:
            parts = line.split()
            if len(parts) < 15:
                continue
            try:
                plate_pair = parts[1].lstrip(":")
                lon1 = float(parts[2])
                lat1 = float(parts[3])
                lon2 = float(parts[4])
                lat2 = float(parts[5])
                length_km = float(parts[6])
                azimuth = float(parts[7])

                # Boundary type is the last field, strip : and *
                btype = parts[-1].strip(":").strip("*")
                if btype not in PB2002_BOUNDARY_TYPES:
                    btype = "unknown"

                setting = PB2002_SETTING_MAP.get(btype, "unknown")

                records.append({
                    "plate_pair": plate_pair,
                    "lon1": lon1, "lat1": lat1,
                    "lon2": lon2, "lat2": lat2,
                    "length_km": length_km,
                    "azimuth": azimuth,
                    "boundary_type": btype,
                    "tectonic_setting": setting,
                })
            except (ValueError, IndexError):
                continue

    return pd.DataFrame(records)


def classify_tectonic_setting(lat: float, lon: float,
                              boundaries_df: pd.DataFrame,
                              radius_deg: float = 2.0) -> str:
    """Classify a point's tectonic setting based on nearest PB2002 boundary.

    Finds all boundary segments within radius_deg of the point and returns
    the most common tectonic setting. Returns 'intraplate' if no boundaries
    are found within the radius.
    """
    # Use midpoint of each segment
    mid_lon = (boundaries_df["lon1"] + boundaries_df["lon2"]) / 2
    mid_lat = (boundaries_df["lat1"] + boundaries_df["lat2"]) / 2

    # Quick rectangular filter (approximate, fast)
    mask = (
        (np.abs(mid_lat - lat) < radius_deg) &
        (np.abs(mid_lon - lon) < radius_deg * max(1.0, 1.0 / np.cos(np.radians(lat + 0.01))))
    )
    nearby = boundaries_df[mask]

    if len(nearby) == 0:
        return "intraplate"

    # Return most common setting
    return nearby["tectonic_setting"].mode().iloc[0]


def classify_grid_tectonic_settings(grid_lats: np.ndarray, grid_lons: np.ndarray,
                                     boundaries_df: pd.DataFrame = None,
                                     radius_deg: float = 2.0) -> list[str]:
    """Classify an array of grid cell centers by their tectonic setting.

    Returns list of setting strings: 'subduction', 'convergent', 'transform',
    'spreading', 'rift', or 'intraplate'.
    """
    if boundaries_df is None:
        boundaries_df = load_pb2002_boundaries()

    settings = []
    for lat, lon in zip(grid_lats, grid_lons):
        settings.append(classify_tectonic_setting(lat, lon, boundaries_df, radius_deg))
    return settings

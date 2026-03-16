#!/usr/bin/env bash
# download_external_data.sh — Download supplementary datasets for Memory of the Earth
#
# Usage: bash scripts/download_external_data.sh
#
# Idempotent: skips files that already exist. Safe to re-run.
#
# NOTE: OCC UIC injection well data (used in NB04) must be downloaded manually from:
#   https://oklahoma.gov/occ/divisions/oil-gas/oil-gas-data.html
# Download annual injection files (occ_uic_2006_2010.xls, occ_uic_2011.xlsx through
# occ_uic_2024.xlsx) and place them in data/external/.

set -euo pipefail

DATA_DIR="data/external"
mkdir -p "$DATA_DIR"

ERRORS=0

download() {
    local url="$1"
    local dest="$2"
    local label="$3"

    if [[ -f "$dest" ]]; then
        echo "[skip] $label — already exists"
        return 0
    fi

    echo "[download] $label ..."
    if curl -fSL --retry 3 --retry-delay 5 -o "$dest" "$url"; then
        echo "[ok] $label"
    else
        echo "[WARN] $label — download failed (curl exit $?)"
        rm -f "$dest"
        ERRORS=$((ERRORS + 1))
    fi
}

# ──────────────────────────────────────────────────────────────────────
# 1. Global CMT Catalog (NDK format) — used in NB01
# ──────────────────────────────────────────────────────────────────────
echo ""
echo "=== Global CMT Catalog ==="

download \
    "https://www.ldeo.columbia.edu/~gcmt/projects/CMT/catalog/jan76_dec20.ndk.gz" \
    "$DATA_DIR/gcmt_jan76_dec20.ndk.gz" \
    "GCMT full catalog 1976-2020 (gzipped)"

if [[ -f "$DATA_DIR/gcmt_jan76_dec20.ndk.gz" && ! -f "$DATA_DIR/gcmt_jan76_dec20.ndk" ]]; then
    echo "[unzip] Decompressing GCMT catalog ..."
    gunzip -k "$DATA_DIR/gcmt_jan76_dec20.ndk.gz"
fi

download \
    "https://www.ldeo.columbia.edu/~gcmt/projects/CMT/catalog/NEW_QUICK/qcmt.ndk" \
    "$DATA_DIR/gcmt_quick_full.ndk" \
    "GCMT quick CMT solutions"

if [[ -f "$DATA_DIR/gcmt_jan76_dec20.ndk" && -f "$DATA_DIR/gcmt_quick_full.ndk" && ! -f "$DATA_DIR/gcmt_combined.ndk" ]]; then
    echo "[combine] Merging GCMT catalogs ..."
    cat "$DATA_DIR/gcmt_jan76_dec20.ndk" "$DATA_DIR/gcmt_quick_full.ndk" \
        > "$DATA_DIR/gcmt_combined.ndk"
fi

# ──────────────────────────────────────────────────────────────────────
# 2. GSRM v2.1 Strain Rate Model — used in NB01, NB03
# ──────────────────────────────────────────────────────────────────────
echo ""
echo "=== GSRM v2.1 Strain Rate Model ==="

download \
    "http://geodesy.unr.edu/GSRM/model/GSRM_average_strain_v2.1.txt.Z" \
    "$DATA_DIR/GSRM_average_strain_v2.1.txt.Z" \
    "GSRM v2.1 strain rate grid (compressed)"

if [[ -f "$DATA_DIR/GSRM_average_strain_v2.1.txt.Z" && ! -f "$DATA_DIR/GSRM_average_strain_v2.1.txt" ]]; then
    echo "[unzip] Decompressing GSRM strain rate grid ..."
    uncompress -k "$DATA_DIR/GSRM_average_strain_v2.1.txt.Z" 2>/dev/null \
        || gzip -d -k "$DATA_DIR/GSRM_average_strain_v2.1.txt.Z" 2>/dev/null \
        || echo "[WARN] Could not decompress GSRM file — try manually: uncompress $DATA_DIR/GSRM_average_strain_v2.1.txt.Z"
fi

# ──────────────────────────────────────────────────────────────────────
# 3. IHFC Global Heat Flow Database (2024 Release) — used in NB01
# ──────────────────────────────────────────────────────────────────────
echo ""
echo "=== IHFC Global Heat Flow Database (2024) ==="

download \
    "https://datapub.gfz.de/download/10.5880.FIDGEO.2024.014-VEueRf/GHFBD-R2024_v.2026-03.zip" \
    "$DATA_DIR/GHFDB-R2024.zip" \
    "IHFC GHFDB 2024 release (ZIP)"

if [[ -f "$DATA_DIR/GHFDB-R2024.zip" ]]; then
    echo "[unzip] Extracting GHFDB ..."
    unzip -o -q "$DATA_DIR/GHFDB-R2024.zip" -d "$DATA_DIR/"
fi

# ──────────────────────────────────────────────────────────────────────
# 4. PB2002 Plate Boundary Model — used in NB01, NB02, NB03, NB05
# ──────────────────────────────────────────────────────────────────────
echo ""
echo "=== PB2002 Plate Boundary Model ==="

download \
    "http://peterbird.name/oldFTP/PB2002/PB2002_boundaries.dig.txt" \
    "$DATA_DIR/PB2002_boundaries.dig.txt" \
    "PB2002 plate boundaries"

download \
    "http://peterbird.name/oldFTP/PB2002/PB2002_plates.dig.txt" \
    "$DATA_DIR/PB2002_plates.dig.txt" \
    "PB2002 plate polygons"

download \
    "http://peterbird.name/oldFTP/PB2002/PB2002_steps.dat.txt" \
    "$DATA_DIR/PB2002_steps.dat.txt" \
    "PB2002 boundary segments"

# ──────────────────────────────────────────────────────────────────────
# 5. SCEDC Southern California Catalog — used in NB05
# ──────────────────────────────────────────────────────────────────────
echo ""
echo "=== SCEDC Southern California Catalog ==="

for START in 2000 2005 2010 2015 2020; do
    END=$((START + 5))
    DEST="$DATA_DIR/scedc_${START}_${END}.txt"
    download \
        "https://service.scedc.caltech.edu/fdsnws/event/1/query?starttime=${START}-01-01&endtime=${END}-01-01&minmagnitude=0.5&format=text" \
        "$DEST" \
        "SCEDC catalog ${START}-${END}"
done

if [[ -f "$DATA_DIR/scedc_2000_2005.txt" && ! -f "$DATA_DIR/scedc_combined.txt" ]]; then
    echo "[combine] Merging SCEDC catalog chunks ..."
    head -1 "$DATA_DIR/scedc_2000_2005.txt" > "$DATA_DIR/scedc_combined.txt"
    for f in "$DATA_DIR"/scedc_*_*.txt; do
        tail -n +2 "$f" >> "$DATA_DIR/scedc_combined.txt"
    done
fi

# ──────────────────────────────────────────────────────────────────────
# 6. USGS Oklahoma Full Catalog — used in NB04
#    Paginated monthly to avoid the 20,000-event API limit
# ──────────────────────────────────────────────────────────────────────
echo ""
echo "=== USGS Oklahoma Full Catalog (M0.1+, paginated) ==="

OK_DEST="$DATA_DIR/ogs_usgs_oklahoma_full.csv"
if [[ -f "$OK_DEST" ]]; then
    echo "[skip] Oklahoma full catalog — already exists"
else
    OK_TMP_DIR="$DATA_DIR/oklahoma_full_tmp"
    mkdir -p "$OK_TMP_DIR"
    FIRST=1

    for YEAR in $(seq 2000 2025); do
        for MONTH in $(seq 1 12); do
            MONTH_PAD=$(printf "%02d" "$MONTH")
            CHUNK="$OK_TMP_DIR/ok_${YEAR}_${MONTH_PAD}.csv"

            if [[ -f "$CHUNK" ]]; then
                continue
            fi

            # Compute end date (first day of next month)
            if [[ "$MONTH" -eq 12 ]]; then
                END_YEAR=$((YEAR + 1))
                END_MONTH="01"
            else
                END_YEAR=$YEAR
                END_MONTH=$(printf "%02d" $((MONTH + 1)))
            fi

            URL="https://earthquake.usgs.gov/fdsnws/event/1/query?format=csv&starttime=${YEAR}-${MONTH_PAD}-01&endtime=${END_YEAR}-${END_MONTH}-01&minlatitude=33.5&maxlatitude=37.5&minlongitude=-100.0&maxlongitude=-94.5&minmagnitude=0.1&orderby=time-asc"

            if curl -fsSL --retry 3 --retry-delay 5 -o "$CHUNK" "$URL"; then
                LINES=$(wc -l < "$CHUNK" | tr -d ' ')
                echo "[ok] Oklahoma ${YEAR}-${MONTH_PAD}: $((LINES - 1)) events"
            else
                echo "[WARN] Oklahoma ${YEAR}-${MONTH_PAD} — download failed"
                rm -f "$CHUNK"
                ERRORS=$((ERRORS + 1))
            fi

            sleep 0.5
        done
    done

    # Combine monthly chunks into single CSV
    echo "[combine] Merging Oklahoma monthly chunks ..."
    FIRST=1
    for f in "$OK_TMP_DIR"/ok_*.csv; do
        if [[ "$FIRST" -eq 1 ]]; then
            cat "$f" > "$OK_DEST"
            FIRST=0
        else
            tail -n +2 "$f" >> "$OK_DEST"
        fi
    done

    TOTAL=$(( $(wc -l < "$OK_DEST" | tr -d ' ') - 1 ))
    echo "[ok] Oklahoma full catalog: $TOTAL total events"
fi

# ──────────────────────────────────────────────────────────────────────
# 7. OCC Injection Well Volumes — manual download required
# ──────────────────────────────────────────────────────────────────────
echo ""
echo "=== OCC Injection Well Volumes ==="
echo "[NOTE] OCC UIC injection data must be downloaded manually."
echo "       Visit: https://oklahoma.gov/occ/divisions/oil-gas/oil-gas-data.html"
echo "       Download annual injection files (occ_uic_2006_2010.xls, occ_uic_2011.xlsx"
echo "       through occ_uic_2024.xlsx) and place them in $DATA_DIR/"

# ──────────────────────────────────────────────────────────────────────
# Summary
# ──────────────────────────────────────────────────────────────────────
echo ""
echo "========================================"
if [[ "$ERRORS" -gt 0 ]]; then
    echo "Done with $ERRORS warning(s). Check output above."
else
    echo "All downloads completed successfully."
fi
echo "========================================"

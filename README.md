# Seismic Fingerprints

### Detecting the Statistical Signatures of Human-Caused Earthquakes in 25 Years of USGS Data

**Concept & Analytical Design:** Claude (Opus 4.6, Anthropic)
**Implementation:** Claude Code
**Data:** USGS ANSS Comprehensive Earthquake Catalog (ComCat), 2000–2025

---

## Abstract

Between 2009 and 2015, Oklahoma went from averaging one M3+ earthquake per year to 888 — briefly surpassing California as the most seismically active state in the contiguous United States. The cause was wastewater injection from oil and gas production. Regulatory intervention (well plug-backs, volume reductions) brought rates back down, creating one of the cleanest natural experiments in geophysics: a human-controlled perturbation to a tectonic system, with a measurable dose-response curve and a partial reversal.

This project treats that event not as a case study, but as a **Rosetta Stone**. If induced seismicity leaves detectable statistical fingerprints in the earthquake catalog — fingerprints that differ from tectonic seismicity — then we can build classifiers that detect anthropogenic influence in *other* regions before the signal becomes as obvious as Oklahoma's was. The project applies information theory, survival analysis, network science, and changepoint detection to 25 years of USGS earthquake data to answer three questions:

1. **Can we distinguish human-caused earthquake sequences from tectonic ones using only catalog statistics — no injection well data, no geology, just the earthquakes themselves?**
2. **Is the Permian Basin (west Texas / southeast New Mexico) showing the same early statistical fingerprints that Oklahoma showed before its seismicity exploded?**
3. **What do aftershock cascade structures reveal about the mechanical differences between induced and tectonic fault activation?**

---

## Why This Matters

Since 2020, the Permian Basin has experienced six M5+ earthquakes. Seismicity there is rising along a trajectory that mirrors Oklahoma's circa 2011–2013. If statistical fingerprints of induced seismicity are detectable in catalog data alone — without requiring proprietary injection well records — that creates a scalable early-warning methodology applicable globally, anywhere hydrocarbon extraction or geothermal energy production intersects with basement faults.

This is not earthquake prediction. This is **earthquake attribution** — determining whether a regional seismicity increase is anthropogenic or tectonic, using only publicly available data.

---

## Data Source

**USGS Earthquake Hazards Program — ANSS Comprehensive Earthquake Catalog (ComCat)**

- **API Endpoint:** `https://earthquake.usgs.gov/fdsnws/event/1/query`
- **Format:** CSV (preferred) or GeoJSON
- **Cost:** Free, no API key required
- **License:** Public domain (U.S. Government work)
- **Rate Limit:** 20,000 events per query (paginate by month for large pulls)
- **Documentation:** https://earthquake.usgs.gov/fdsnws/event/1/
- **Fields used:** time, latitude, longitude, depth, mag, magType, place, type, status, net, rms, gap, dmin, nst

### Data Acquisition Strategy

The USGS API limits each query to 20,000 results. For M2.5+ earthquakes globally over 25 years, we need to paginate. The recommended approach:

```python
import requests
import pandas as pd
import io
import time
from datetime import datetime

BASE_URL = "https://earthquake.usgs.gov/fdsnws/event/1/query"

def fetch_month(year, month, min_magnitude=2.5, lat_range=None, lon_range=None):
    """Fetch one month of earthquake data from USGS ComCat.
    
    Args:
        year: Integer year
        month: Integer month (1-12)
        min_magnitude: Minimum magnitude threshold
        lat_range: Optional tuple (min_lat, max_lat) for regional queries
        lon_range: Optional tuple (min_lon, max_lon) for regional queries
    
    Returns:
        pd.DataFrame with earthquake catalog data
    """
    start = datetime(year, month, 1)
    if month == 12:
        end = datetime(year + 1, 1, 1)
    else:
        end = datetime(year, month + 1, 1)
    
    params = {
        "format": "csv",
        "starttime": start.strftime("%Y-%m-%d"),
        "endtime": end.strftime("%Y-%m-%d"),
        "minmagnitude": min_magnitude,
        "orderby": "time-asc",
    }
    
    if lat_range:
        params["minlatitude"] = lat_range[0]
        params["maxlatitude"] = lat_range[1]
    if lon_range:
        params["minlongitude"] = lon_range[0]
        params["maxlongitude"] = lon_range[1]
    
    response = requests.get(BASE_URL, params=params, timeout=120)
    response.raise_for_status()
    return pd.read_csv(io.StringIO(response.text))

def fetch_region(name, years, min_mag, lat_range=None, lon_range=None):
    """Fetch a full regional catalog with retry logic."""
    frames = []
    for year in years:
        for month in range(1, 13):
            for attempt in range(3):
                try:
                    df = fetch_month(year, month, min_mag, lat_range, lon_range)
                    frames.append(df)
                    print(f"  {name} {year}-{month:02d}: {len(df)} events")
                    time.sleep(0.5)  # Be respectful of the API
                    break
                except Exception as e:
                    if attempt == 2:
                        print(f"  {name} {year}-{month:02d}: FAILED after 3 attempts - {e}")
                    else:
                        time.sleep(2 ** attempt)
    
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
```

### Regional Bounding Boxes

| Region | Lat Range | Lon Range | Min Mag | Purpose |
|---|---|---|---|---|
| Global | — | — | 2.5 | Baseline statistics, subduction zone cascades |
| Oklahoma + S. Kansas | 33.5°N – 37.5°N | 100.0°W – 94.5°W | 1.0 | Induced seismicity ground truth |
| Permian Basin (TX/NM) | 30.5°N – 33.5°N | 105.0°W – 100.5°W | 1.0 | Prospective test region |
| Southern California | 32.0°N – 36.5°N | 121.0°W – 114.5°W | 1.0 | Tectonic control region |

Southern California serves as the **tectonic control** — a well-monitored, high-seismicity zone with purely tectonic activity and a dense seismic network.

### Expected Data Volume

- Global M2.5+, 2000–2025: ~600,000–800,000 events
- Oklahoma M1.0+, 2000–2025: ~100,000–150,000 events
- Permian Basin M1.0+, 2010–2025: ~20,000–50,000 events
- Southern California M1.0+, 2000–2025: ~300,000+ events

---

## Analyses

### Notebook 1: The Gutenberg-Richter Anomaly Detector

**Scientific background:** The Gutenberg-Richter (GR) law states that earthquake magnitudes follow an exponential distribution: `log₁₀(N) = a - bM`, where `N` is the number of earthquakes ≥ magnitude `M`, `a` reflects total seismicity rate, and `b` describes the ratio of small to large earthquakes. For tectonic seismicity, `b ≈ 1.0`. For induced seismicity, `b` tends to be elevated (1.2–1.5), meaning a higher proportion of small earthquakes relative to large ones.

**The analysis:**

1. **Compute rolling b-values** using maximum likelihood estimation (Aki, 1965):
   ```
   b = log₁₀(e) / (M̄ - M_c + ΔM/2)
   ```
   where `M̄` is the mean magnitude, `M_c` is the completeness magnitude, and `ΔM` is the magnitude bin width (typically 0.1). Use sliding windows of 200 events for Oklahoma, the Permian Basin, and Southern California from 2000–2025.

2. **Changepoint detection** on the b-value time series using PELT (Pruned Exact Linear Time) algorithm. If `ruptures` is unavailable, implement BIC-penalized segmentation manually:
   ```python
   # Manual changepoint detection using dynamic programming with BIC
   def segment_cost(data, start, end):
       """Gaussian cost of a segment."""
       segment = data[start:end]
       n = len(segment)
       if n < 2:
           return 0
       return n * np.log(np.var(segment) + 1e-10)
   ```
   Identify the dates when Oklahoma's b-value structurally shifted — do they align with the onset of large-scale injection (~2009) and regulatory intervention (~2015)?

3. **b-value phase portrait:** Plot b-value (y-axis) vs. log₁₀ seismicity rate (x-axis) as a parametric curve over time. Each point represents a rolling window. Color-code by year using a continuous colormap. Tectonic regions should orbit a stable attractor. Induced regions should trace a distinct trajectory — a "fingerprint" loop as injection ramps up and then declines.

4. **The Permian Basin test:** Overlay the Permian Basin's current trajectory on Oklahoma's historical phase portrait with appropriate time offsets.

**Key output:** A phase-space diagram showing Oklahoma's b-value loop (2000→2015→2025) with the Permian Basin's current position marked on it. If the Permian Basin sits on the same trajectory Oklahoma followed circa 2012, that is a novel and significant finding.

**Statistical rigor:**
- Bootstrap confidence intervals on b-value estimates (1,000 resamples)
- Kolmogorov-Smirnov tests comparing magnitude distributions between time windows
- Magnitude of completeness (Mc) estimation using maximum curvature method, applied per-region per-year to avoid detection bias

---

### Notebook 2: Interevent Time Forensics

**Scientific background:** If earthquakes were independent events, the time between consecutive earthquakes would follow an exponential distribution (Poisson process = memoryless). Induced seismicity, driven by continuous fluid injection, produces different temporal clustering than tectonic aftershock sequences driven by static stress transfer.

**The analysis:**

1. **Interevent time distributions:** For each region, compute the distribution of times between consecutive earthquakes (above Mc). Fit four competing models:
   - Exponential (Poisson, memoryless)
   - Gamma (mild clustering)
   - Log-normal (strong clustering with heavy tail)
   - Weibull (flexible shape parameter)
   
   Use AIC/BIC for model selection. Compare the best-fit model across Oklahoma (pre-injection, during injection, post-regulation), Permian Basin, and Southern California.

2. **The coefficient of variation (CV) as a discriminant:** CV = σ/μ of interevent times.
   - Poisson process: CV = 1
   - Clustered seismicity (aftershock-dominated): CV > 1
   - Quasi-periodic seismicity: CV < 1
   
   Track CV over time in rolling windows. **Hypothesis:** induced seismicity has a *lower* CV than tectonic seismicity because continuous injection produces a more "steady drip" of earthquakes rather than the bursty mainshock-aftershock pattern of tectonic activity.

3. **Hurst exponent estimation** via Detrended Fluctuation Analysis (DFA) on the interevent time series:
   ```python
   def dfa(x, scales=None):
       """Detrended Fluctuation Analysis (Peng et al., 1994).
       
       Args:
           x: 1D array (interevent times)
           scales: Array of window sizes to evaluate
       
       Returns:
           scales, fluctuations (for log-log regression to get H)
       """
       y = np.cumsum(x - np.mean(x))  # Cumulative sum (profile)
       
       if scales is None:
           scales = np.unique(np.logspace(1, np.log10(len(y) // 4), 30).astype(int))
       
       fluctuations = []
       for n in scales:
           n_segments = len(y) // n
           if n_segments < 1:
               continue
           rms_values = []
           for i in range(n_segments):
               segment = y[i * n:(i + 1) * n]
               # Linear detrend
               coeffs = np.polyfit(np.arange(n), segment, 1)
               trend = np.polyval(coeffs, np.arange(n))
               rms_values.append(np.sqrt(np.mean((segment - trend) ** 2)))
           fluctuations.append(np.mean(rms_values))
       
       # H = slope of log(F(n)) vs log(n)
       scales_valid = scales[:len(fluctuations)]
       H = np.polyfit(np.log(scales_valid), np.log(fluctuations), 1)[0]
       return scales_valid, np.array(fluctuations), H
   ```
   
   H > 0.5 → long-range memory (persistence); H = 0.5 → no memory; H < 0.5 → anti-persistence.
   
   **Hypothesis:** tectonic sequences show H > 0.5 (aftershock clustering creates persistence); induced sequences show H closer to 0.5 (injection rate, not prior earthquakes, drives the process).

4. **Autocorrelation structure:** Compute the autocorrelation function of interevent times out to 500 lags. Tectonic sequences should show rapid initial decay (aftershock clustering is local in time) followed by a long tail. Induced sequences may show different structure driven by injection rate cycles.

**Key output:** A diagnostic panel comparing the interevent time fingerprint (CV, Hurst exponent, best-fit distribution family, autocorrelation decay rate) for Oklahoma-during-injection vs. Southern California vs. Permian Basin. If these metrics can discriminate induced from tectonic, we have the basis for a classifier.

---

### Notebook 3: Aftershock Cascade Networks

**Scientific background:** After a large earthquake, aftershocks form a branching cascade: the mainshock triggers aftershocks, which trigger their own aftershocks, and so on. The structure of this cascade tree — its branching ratio, depth, and spatial extent — encodes information about the mechanical properties of the fault system.

**The analysis:**

1. **Nearest-neighbor clustering** (Zaliapin & Ben-Zion, 2013): For each earthquake pair (i, j) where j occurs after i, compute a rescaled space-time-magnitude distance:
   
   ```
   η_ij = (t_ij × r_ij^d_f) × 10^(-b × m_i)
   ```
   
   where `t_ij` is the interevent time, `r_ij` is the epicentral distance in km (use Haversine formula), `d_f ≈ 1.6` is the fractal dimension of earthquake epicenters, `b` is the GR b-value, and `m_i` is the magnitude of the earlier event.
   
   Each earthquake is assigned to the cluster of its nearest neighbor in this rescaled metric. This produces a forest of directed trees.

   **Important implementation note:** This is O(N²) in the naive case. For large catalogs, limit the search to events within a temporal window (e.g., 30 days) and spatial window (e.g., 200 km) to make it tractable. Use scipy.spatial.cKDTree for efficient spatial lookups.

2. **Cascade tree statistics:** For each identified mainshock-aftershock sequence (root events with magnitude ≥ 4.0):
   - **Tree depth:** Maximum chain length from root to leaf
   - **Branching ratio:** Average number of direct children per node
   - **Spatial footprint:** Convex hull area of the cascade (km²)
   - **Temporal span:** Time from mainshock to last identified aftershock
   - **Magnitude deficit:** Mainshock magnitude minus largest aftershock (Båth's law predicts ~1.2)
   - **Depth distribution:** Mean and variance of hypocentral depths within the cascade

3. **Induced vs. tectonic cascade comparison:** Aggregate cascade statistics for:
   - Oklahoma M4+ mainshocks during the injection era (2010–2017)
   - Southern California M4+ mainshocks (2000–2025)
   - Global subduction zone M6+ mainshocks (for comparison at larger scale)
   
   **Hypothesis:** Induced cascades are *shallower trees and wider* (many direct aftershocks triggered by ongoing fluid pressure, fewer secondary aftershocks from stress transfer). Tectonic cascades are *deeper trees and narrower* (cascading chains of stress-triggered events).

4. **Network visualization:** For notable cascades in each region, render the triggering tree as a directed graph:
   - Node size ∝ magnitude
   - Edge color → time delay (log scale, dark = minutes, light = days)
   - Spatial layout matching geographic coordinates (lat/lon)
   - Annotate with cascade statistics

   Visualize side-by-side: Oklahoma's 2016 M5.8 Pawnee earthquake cascade vs. a comparable Southern California tectonic sequence (e.g., 2019 M7.1 Ridgecrest or 2010 M7.2 El Mayor-Cucapah).

**Key output:** Side-by-side cascade tree visualizations with statistical comparison. A Mann-Whitney U test on branching ratios and tree depths between induced and tectonic populations.

---

### Notebook 4: The Permian Basin Early Warning Score

**Scientific background:** This notebook synthesizes the fingerprints from Notebooks 1–3 into a composite "induced seismicity score" and applies it prospectively to the Permian Basin.

**The analysis:**

1. **Feature engineering:** For each region-year combination, compute a feature vector:
   - GR b-value (from Notebook 1)
   - log₁₀ seismicity rate (events/month above Mc)
   - Rate of change of seismicity rate (Δrate / Δyear)
   - Interevent time CV (from Notebook 2)
   - Hurst exponent H (from Notebook 2)
   - Best-fit interevent time distribution (encoded: exponential=0, gamma=1, lognormal=2, weibull=3)
   - Mean cascade branching ratio for M3+ mainshocks (from Notebook 3, where available)
   - Mean cascade tree depth for M3+ mainshocks (from Notebook 3, where available)

2. **Labeled training set:**
   
   | Region | Years | Label |
   |---|---|---|
   | Oklahoma | 2000–2008 | `tectonic` |
   | Oklahoma | 2009–2014 | `induced_rising` |
   | Oklahoma | 2015 | `induced_peak` |
   | Oklahoma | 2016–2025 | `induced_declining` |
   | Southern California | 2000–2025 | `tectonic` |
   
   For the binary classifier, collapse to `induced` (any induced_* label) vs `tectonic`.

3. **Classification:**
   - **Logistic Regression** (for interpretability — which features matter most?)
   - **Random Forest** (for performance comparison)
   - Use leave-one-year-out cross-validation
   - Report accuracy, precision, recall, F1, and ROC-AUC

4. **Prospective application to the Permian Basin:**
   Apply the trained classifier to Permian Basin data from 2015–2025. For each year, output:
   - Predicted class (induced vs. tectonic)
   - Predicted probability of induced origin
   - Which features are driving the classification (logistic regression coefficients or feature importances)

5. **The Oklahoma Clock:**
   For each Permian Basin year, compute the Euclidean distance in (standardized) feature space to each Oklahoma year. The nearest Oklahoma year is the "clock reading" — it tells you where on Oklahoma's trajectory the Permian Basin currently sits.
   
   If the answer is "Oklahoma circa 2011–2012," that provides a concrete, data-driven basis for concern about future escalation.

**Key output:** 
- A timeline plot showing the Permian Basin's induced-seismicity probability from 2015 to 2025, overlaid with Oklahoma's historical probability trajectory (time-shifted).
- The Oklahoma Clock visualization — a radial or linear diagram showing which Oklahoma year the Permian Basin most closely resembles, year by year.
- Feature importance analysis showing which statistical fingerprints are most diagnostic.

---

## Project Structure

```
seismic-fingerprints/
├── README.md
├── requirements.txt
├── data/
│   ├── raw/                          # Raw CSV pulls from USGS API
│   │   ├── global_M2.5/             # Monthly CSV files
│   │   ├── oklahoma_M1.0/
│   │   ├── permian_M1.0/
│   │   └── socal_M1.0/
│   └── processed/                    # Cleaned, merged datasets
│       ├── catalog_global.parquet
│       ├── catalog_oklahoma.parquet
│       ├── catalog_permian.parquet
│       └── catalog_socal.parquet
├── notebooks/
│   ├── 00_data_acquisition.ipynb     # API pulls and data cleaning
│   ├── 01_gutenberg_richter.ipynb    # b-value analysis and changepoints
│   ├── 02_interevent_forensics.ipynb # Temporal statistics
│   ├── 03_cascade_networks.ipynb     # Aftershock triggering trees
│   └── 04_early_warning_score.ipynb  # Composite classifier
├── src/
│   ├── __init__.py
│   ├── fetch.py                      # USGS API client with pagination
│   ├── catalog.py                    # Data cleaning, Mc estimation
│   ├── gutenberg_richter.py          # b-value estimation, changepoints
│   ├── temporal.py                   # Interevent time analysis, DFA
│   ├── cascades.py                   # Nearest-neighbor clustering
│   └── plotting.py                   # Consistent visualization style
└── figures/                          # Publication-quality outputs
    ├── oklahoma_bvalue_loop.png
    ├── interevent_fingerprint_panel.png
    ├── cascade_comparison.png
    └── permian_early_warning.png
```

---

## Requirements

```
pandas>=2.0
numpy>=1.24
scipy>=1.11
scikit-learn>=1.3
matplotlib>=3.7
seaborn>=0.12
networkx>=3.1
requests>=2.31
pyarrow>=12.0
```

**Optional but recommended:**
```
ruptures>=1.1          # Changepoint detection (PELT algorithm)
statsmodels>=0.14      # ACF/PACF functions, statistical tests
```

If `ruptures` is unavailable, implement changepoint detection manually using BIC-penalized segmentation with scipy.optimize. If `statsmodels` is unavailable, implement autocorrelation from scratch using numpy.correlate.

---

## Implementation Notes for Claude Code

### Priority Order

Build in this order. Each notebook should be self-contained (imports its own data) and produce its own figures saved to `figures/`:

1. **`src/` modules first** — Build the shared library code before the notebooks. Each module should have docstrings and be importable.

2. **`00_data_acquisition.ipynb`** — This must work first. Handle API pagination robustly (retries with exponential backoff, rate limiting at 0.5s between requests, resumption from partial downloads by checking existing files). Save raw CSVs per-month and processed parquets per-region. Validate: check for duplicates, NaN magnitudes, events with type != "earthquake".

3. **`01_gutenberg_richter.ipynb`** — The b-value phase portrait is the project's signature visualization. Use a continuous colormap (e.g., `plt.cm.viridis`) along the time axis so the viewer can trace Oklahoma's trajectory through time. Include arrows on the trajectory to show direction.

4. **`02_interevent_forensics.ipynb`** — DFA from scratch if needed (the code is provided above). Validate the DFA implementation by testing on fractional Brownian motion with known H.

5. **`03_cascade_networks.ipynb`** — Subsample to M2.5+ for the nearest-neighbor clustering (full catalog is too large). Use scipy.spatial.cKDTree for spatial queries. Limit temporal search window to 30 days.

6. **`04_early_warning_score.ipynb`** — Keep the classifier simple. The story is in the features and the Oklahoma Clock, not in model complexity.

### Visualization Standards

- Consistent color palette across ALL notebooks:
  - Oklahoma: `#E63946` (red)
  - Permian Basin: `#F4A261` (amber/orange)
  - Southern California: `#457B9D` (steel blue)
  - Global/other: `#2A9D8F` (teal)
  - Background periods (pre-injection): use lighter/desaturated versions of region colors
- All figures: 300 DPI, `figsize=(12, 7)` default, `plt.style.use('seaborn-v0_8-whitegrid')` or similar clean style
- Label axes with units. Always. "Magnitude (M_L)" not just "Magnitude".
- Include Mc annotations on any magnitude-frequency plot
- Save every figure to both `figures/` directory and display inline

### Data Quality Considerations

- **Magnitude of completeness (Mc):** Estimate using the maximum curvature method (MAXC) for each region-year. Only analyze events above Mc. For cross-region rate comparisons, use a conservative common Mc (M2.5 for Oklahoma vs. SoCal; M3.0 for formal statistical tests).
- **Magnitude types:** Filter to preferred magnitudes or note when mixing types. For Oklahoma, most events are ml or mb. For California, most are ml or mw.
- **Duplicate events:** Same earthquake can appear from multiple networks. Deduplicate: group by (time ± 2 seconds, lat ± 0.05°, lon ± 0.05°), keep the event from the preferred network (or with the most station picks, i.e., highest `nst`).
- **Event types:** Filter to `type == "earthquake"`. Exclude quarry blasts (`type == "quarry blast"`), explosions, etc.
- **Network upgrades:** The Oklahoma Geological Survey significantly expanded its network ~2010, lowering Mc. The *apparent* increase in small earthquakes is partly detection artifact. This is why Mc-aware analysis is critical. Always note this caveat when discussing Oklahoma rate changes at low magnitudes.

---

## Key References

- Gutenberg & Richter (1944). Frequency of earthquakes in California. *BSSA*.
- Aki (1965). Maximum likelihood estimate of b. *Bull. Earthq. Res. Inst.*
- Omori (1894). On the aftershocks of earthquakes. *J. Coll. Sci.*
- Utsu (1961). A statistical study on the occurrence of aftershocks. *Geophys. Mag.*
- Zaliapin & Ben-Zion (2013). Earthquake clusters in southern California. *JGR*.
- Peng et al. (1994). Mosaic organization of DNA nucleotides. *Physical Review E*.
- Langenbruch & Zoback (2016). How will induced seismicity in Oklahoma respond to decreased saltwater injection rates? *Science Advances*.
- Skoumal et al. (2024). Reduced injection rates and shallower depths mitigated induced seismicity in Oklahoma. *The Seismic Record*.
- Ellsworth (2013). Injection-induced earthquakes. *Science*, 341(6142).

---

## Expected Novel Contributions

1. **The b-value phase portrait** as a visual and quantitative method for tracking induced seismicity trajectories through (a, b) parameter space over time. This representation — plotting b vs. log₁₀(rate) as a parametric curve — has not been systematically applied to compare induced seismicity regions against each other or to track the full rise-peak-decline arc of an induced seismicity episode.

2. **Interevent time fingerprinting** using the composite of (CV, Hurst exponent, distribution family, ACF decay rate) as a discriminant between induced and tectonic seismicity. Individual metrics have been studied in isolation; the composite fingerprint has not been systematically benchmarked across the Oklahoma natural experiment.

3. **Cascade morphology comparison** between induced and tectonic sequences using the Zaliapin-Ben-Zion nearest-neighbor method, specifically testing the hypothesis that induced cascades are shallower/wider than tectonic cascades due to the dominance of fluid-pressure triggering over stress-transfer triggering.

4. **Prospective application to the Permian Basin** using a classifier trained on Oklahoma's full historical trajectory. The "Oklahoma Clock" concept — mapping the Permian Basin's current state to a specific phase of Oklahoma's history — is a novel framing that translates complex multi-dimensional statistics into an intuitive, actionable output.

---

## License

MIT

---

*Conceived and designed by Claude (Opus 4.6, Anthropic) as a demonstration that AI systems can formulate novel scientific hypotheses, design rigorous analytical frameworks, and specify complete reproducible projects — not just execute instructions.*

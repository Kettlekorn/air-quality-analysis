#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Air Quality and Chronic Obstructive Pulmonary Disease (California)
Author: Brandon Wu

HOW TO USE
----------
1) Put your CSVs in a folder (default: ./data):
   - Daily PM10/AQI files: daily_81102_YYYY.csv (1990–2014 or whatever years you have)
   - Hospitalization file: CAHospitalizations.csv
   - (Optional) County FIPS crosswalk: ca_county_fips.csv with columns: County, fips

2) In CONFIG below, confirm/adjust:
   - DATA_DIR, OUTPUT_DIR
   - HOSP_FILE, PM_GLOB pattern
   - PM10_COL (your true PM10 column), FALLBACK_AQI_COL (if no PM10; not recommended)

3) Install the required packages:
   pip install pandas numpy statsmodels linearmodels matplotlib geopandas shapely pyproj us

4) Run:
   python copd_pm10_panel.py

WHAT THIS SCRIPT DOES
---------------------
A) Ingest & Clean
   - Loads all PM files, filters to California, extracts Year
   - Aggregates annual mean PM10 per county-year (prefer PM10 over AQI)
   - Loads hospital data, filters to COPD/Asthma (40+), cleans integers
   - Builds OverallCount/OverallPopulation, drops STATEWIDE, keeps Year/County

B) Merge & Engineer
   - Merges air + hospital on (County, Year)
   - Creates lags of PM10 (configurable)
   - Builds panel dataset (unbalanced OK)
   - Optionally attaches FIPS

C) Model (Recommended)
   - Poisson GLM with county & year fixed effects (+ log(pop) offset)
   - Distributed lags (joint test)
   - Cluster-robust SEs by county
   - Placebo test with future PM10 (should be null)

D) Replicate Your OLS + BH (Optional)
   - Per-county OLS on Count and Count/Population
   - Extract p-values across gap=1..9
   - Benjamini-Hochberg FDR control
   - Frequency barplots by gap and county

E) Outputs
   - CSV: cleaned/merged panels, model tables, BH results
   - PNG: barplots; (optional) CA county choropleth if shapefile is available
"""

# =========================
# CONFIG
# =========================
from pathlib import Path

DATA_DIR    = Path("./data")
OUTPUT_DIR  = Path("./outputs")

# Input file patterns
PM_GLOB     = "daily_81102_*.csv"     # pattern for daily air files in DATA_DIR
HOSP_FILE   = "CAHospitalizations.csv"

# Column names in PM files
PM10_COL            = "Arithmetic.Mean"   # <-- set to your true PM10 concentration column
FALLBACK_AQI_COL    = "AQI"               # only used if PM10_COL missing (not recommended)

# Hospitalization filter
COPD_FILTER_TEXT    = "COPD or Asthma in Older Adults (Age 40+)"

# Lags to include in the panel model
PM_LAGS             = [1, 2]              # e.g., 1-2 year lag (edit as you like)

# Optional FIPS crosswalk CSV in DATA_DIR with columns: County, fips
FIPS_CROSSWALK      = "ca_county_fips.csv"  # set to None if not available

# Replicate OLS + BH section?
RUN_PER_COUNTY_OLS_BH = True

# Plotting options
MAKE_MAP = False  # requires county shapefile; set True if you have it

# =========================
# IMPORTS
# =========================
import sys
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore")


# =========================
# UTILS
# =========================
def ensure_dirs():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

def safe_numeric_series(s):
    """Strip non-digits and cast to float; empty -> NaN."""
    s = s.astype(str).str.replace(r'[^0-9]', '', regex=True)
    s = s.replace('', np.nan).astype(float)
    return s

def benjamini_hochberg(df, p_col="PValue", alpha=0.05):
    """Return a copy of df with BH-adjusted p-values and significance flag."""
    out = df.copy()
    p = out[p_col].values
    p = np.where(np.isfinite(p), p, np.nan)
    mask = np.isfinite(p)
    if mask.sum() == 0:
        out["p_adj"] = np.nan
        out["significant"] = False
        return out
    rej, p_adj, _, _ = multipletests(p[mask], alpha=alpha, method="fdr_bh")
    out["p_adj"] = np.nan
    out.loc[mask, "p_adj"] = p_adj
    out["significant"] = False
    out.loc[mask, "significant"] = rej
    return out

def print_header(title):
    print("\n" + "="*len(title))
    print(title)
    print("="*len(title))


# =========================
# A) LOAD & CLEAN AIR DATA
# =========================
def load_pm10_year(file_path: Path) -> pd.DataFrame:
    df = pd.read_csv(file_path)
    df.columns = df.columns.str.replace(" ", ".", regex=False)
    # Filter to California
    if "State.Name" not in df.columns:
        raise ValueError(f"'State.Name' missing in {file_path.name}")
    df = df[df["State.Name"] == "California"].copy()

    # Extract Year from Date.Local
    if "Date.Local" not in df.columns:
        raise ValueError(f"'Date.Local' missing in {file_path.name}")
    df["Year"] = pd.to_datetime(df["Date.Local"], errors="coerce").dt.year

    # Choose PM10 metric
    col = PM10_COL if PM10_COL in df.columns else None
    if col is None:
        if FALLBACK_AQI_COL in df.columns:
            warnings.warn(f"{file_path.name}: {PM10_COL} not found; using {FALLBACK_AQI_COL} (not recommended).")
            col = FALLBACK_AQI_COL
        else:
            raise ValueError(f"Neither {PM10_COL} nor {FALLBACK_AQI_COL} found in {file_path.name}.")

    # Aggregate annual mean PM10 per county-year
    if "County.Name" not in df.columns:
        raise ValueError(f"'County.Name' missing in {file_path.name}")
    out = (df.groupby(["Year", "County.Name"])
             .agg(PM10=(col, "mean"))
             .reset_index())
    out = out.rename(columns={"County.Name": "County"})
    return out

def build_air_panel():
    files = sorted((DATA_DIR).glob(PM_GLOB))
    if not files:
        raise FileNotFoundError(f"No files matching {PM_GLOB} in {DATA_DIR}")
    air = pd.concat([load_pm10_year(f) for f in files], ignore_index=True)
    air = air.dropna(subset=["Year", "County", "PM10"]).copy()
    air["Year"] = air["Year"].astype(int)
    return air


# ==============================
# B) LOAD & CLEAN HOSPITAL DATA
# ==============================
def load_hospital():
    f = DATA_DIR / HOSP_FILE
    if not f.exists():
        raise FileNotFoundError(f"Hospitalization file not found: {f}")
    hosp = pd.read_csv(f)

    # Filter for COPD/Asthma (40+)
    if "PQIDescription" not in hosp.columns:
        raise ValueError("'PQIDescription' column missing in hospital data.")
    h = hosp.loc[hosp["PQIDescription"].eq(COPD_FILTER_TEXT)].copy()

    # Clean numeric fields
    needed = ["Count_ICD10","Population_ICD10","Count_ICD9","Population_ICD9"]
    for col in needed:
        if col not in h.columns:
            raise ValueError(f"'{col}' missing in hospital data.")
        h[col] = safe_numeric_series(h[col])

    h["OverallCount"] = h["Count_ICD10"].fillna(h["Count_ICD9"])
    h["OverallPopulation"] = h["Population_ICD10"].fillna(h["Population_ICD9"])

    # Keep essentials, drop STATEWIDE and invalid pop
    must = ["Year","County","OverallCount","OverallPopulation"]
    for m in must:
        if m not in h.columns:
            raise ValueError(f"'{m}' missing after cleaning hospital data.")
    h = h[(h["OverallPopulation"].notna()) & (h["OverallPopulation"] > 0)]
    if "County" not in h.columns:
        raise ValueError("'County' missing in hospital data.")
    h = h[h["County"] != "STATEWIDE"]
    h = h[must].copy()
    h["Year"] = h["Year"].astype(int)

    return h


# ======================================
# C) MERGE, LAG FEATURES, ATTACH FIPS
# ======================================
def attach_fips(df):
    if not FIPS_CROSSWALK:
        df["fips"] = np.nan
        return df
    cross = DATA_DIR / FIPS_CROSSWALK
    if not cross.exists():
        warnings.warn(f"FIPS crosswalk not found: {cross}; proceeding without FIPS.")
        df["fips"] = np.nan
        return df
    x = pd.read_csv(cross, dtype={"fips": str})
    if not {"County","fips"}.issubset(x.columns):
        warnings.warn("FIPS crosswalk must have columns: County, fips")
        df["fips"] = np.nan
        return df
    out = df.merge(x[["County","fips"]], on="County", how="left")
    return out

def make_lags(df, var="PM10", group="County", lags=(1,2)):
    df = df.sort_values([group, "Year"]).copy()
    for L in lags:
        df[f"{var}_lag{L}"] = df.groupby(group)[var].shift(L)
    return df

def build_panel():
    air = build_air_panel()
    hosp = load_hospital()

    merged = hosp.merge(air, on=["County","Year"], how="inner")
    merged = make_lags(merged, var="PM10", group="County", lags=PM_LAGS)
    merged["log_pop"] = np.log(merged["OverallPopulation"])
    merged = attach_fips(merged)

    # Save intermediate
    merged.to_csv(OUTPUT_DIR / "merged_panel.csv", index=False)
    return merged


# ======================================
# D) PANEL POISSON GLM WITH FE + LAGS
# ======================================
def fit_panel_poisson(df, lags=PM_LAGS):
    # Drop rows missing any lag
    need = [f"PM10_lag{L}" for L in lags]
    d = df.dropna(subset=need + ["OverallCount","log_pop"]).copy()

    # Build formula: count ~ lags + C(Year) + C(County), offset = log_pop
    rhs = " + ".join([f"PM10_lag{L}" for L in lags] + ["C(Year)","C(County)"])
    formula = f"OverallCount ~ {rhs}"

    print_header("Fitting Poisson GLM with County/Year Fixed Effects and Distributed Lags")
    print("Formula:", formula)
    model = smf.glm(
        formula=formula,
        data=d,
        family=sm.families.Poisson(),
        offset=d["log_pop"]
    ).fit(cov_type="cluster", cov_kwds={"groups": d["County"]})

    # Save summary
    with open(OUTPUT_DIR / "panel_poisson_summary.txt", "w") as f:
        f.write(model.summary().as_text())

    # Wald joint test for all lag coefficients
    if len(lags) > 1:
        constraints = " = 0, ".join([f"PM10_lag{L}" for L in lags]) + " = 0"
        wtest = model.wald_test(constraints)
    else:
        wtest = model.wald_test(f"PM10_lag{lags[0]} = 0")

    # Save wald test
    with open(OUTPUT_DIR / "panel_poisson_lag_wald_test.txt", "w") as f:
        f.write(str(wtest))

    print(model.summary())
    print_header("Joint Wald Test for Lag Block")
    print(wtest)

    # Export coefficients
    coefs = model.params.rename("coef").to_frame()
    coefs["se"] = model.bse
    coefs["pval"] = model.pvalues
    coefs.to_csv(OUTPUT_DIR / "panel_poisson_coefs.csv")

    return model


# ======================================
# E) PLACEBO: FUTURE PM10 SHOULD BE NULL
# ======================================
def placebo_future_lag(df, years_ahead=1):
    # Create lead variable: future PM10 (should not predict present outcomes)
    d = df.sort_values(["County","Year"]).copy()
    d[f"PM10_lead{years_ahead}"] = d.groupby("County")["PM10"].shift(-years_ahead)
    d = d.dropna(subset=[f"PM10_lead{years_ahead}","OverallCount","log_pop"]).copy()

    formula = f"OverallCount ~ PM10_lead{years_ahead} + C(Year) + C(County)"
    print_header(f"Placebo Test: Poisson GLM with Future PM10 (+{years_ahead})")
    print("Formula:", formula)

    model = smf.glm(
        formula=formula,
        data=d,
        family=sm.families.Poisson(),
        offset=d["log_pop"]
    ).fit(cov_type="cluster", cov_kwds={"groups": d["County"]})

    with open(OUTPUT_DIR / f"placebo_lead{years_ahead}_summary.txt", "w") as f:
        f.write(model.summary().as_text())

    print(model.summary())
    return model


# ===================================================
# F) OPTIONAL: REPLICATE PER-COUNTY OLS + BH (R-STYLE)
# ===================================================
from typing import List, Dict, Tuple

def per_county_ols(df, counties: List[str], gaps=range(1,10)) -> pd.DataFrame:
    """
    Replicates your R approach:
      - For each county
      - For each gap in 1..9
      - Fit OLS: (a) OverallCount ~ PM10_lag{gap}
                (b) (OverallCount/Population) ~ PM10_lag{gap}
      - Extract p-values
    Note: This is for comparison; Poisson GLM with offsets is preferable for counts.
    """
    rows = []
    for c in counties:
        cd = df[df["County"] == c].sort_values("Year").copy()

        # Precompute lags up to max(gaps)
        max_gap = max(gaps)
        for L in range(1, max_gap+1):
            col = f"PM10_lag{L}"
            if col not in cd.columns:
                cd[col] = cd["PM10"].shift(L)

        # Build rate safely
        cd["Rate"] = cd["OverallCount"] / cd["OverallPopulation"]

        for gap in gaps:
            lag_col = f"PM10_lag{gap}"
            sub = cd.dropna(subset=[lag_col, "OverallCount", "Rate"])

            if len(sub) < 5:
                continue

            # (a) Count model
            try:
                m1 = smf.ols(f"OverallCount ~ {lag_col}", data=sub).fit()
                p1 = m1.pvalues.get(lag_col, np.nan)
            except Exception:
                p1 = np.nan

            # (b) Rate model
            try:
                m2 = smf.ols(f"Rate ~ {lag_col}", data=sub).fit()
                p2 = m2.pvalues.get(lag_col, np.nan)
            except Exception:
                p2 = np.nan

            rows.append({"County": c, "GapPeriod": gap, "Model": "Count", "PValue": p1})
            rows.append({"County": c, "GapPeriod": gap, "Model": "Count/Population", "PValue": p2})

    out = pd.DataFrame(rows)
    out.to_csv(OUTPUT_DIR / "per_county_ols_raw_pvals.csv", index=False)
    return out

def plot_gap_hist(df_sig: pd.DataFrame, model_label: str, fname: str):
    counts = df_sig["GapPeriod"].value_counts().sort_index()
    plt.figure()
    counts.plot(kind="bar")
    plt.xlabel("Gap Period")
    plt.ylabel("Frequency")
    plt.title(f"Significant P-Values ({model_label})")
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / fname, dpi=160)
    plt.close()

def plot_county_freq(df_sig: pd.DataFrame, model_label: str, fname: str):
    counts = df_sig["County"].value_counts()
    plt.figure(figsize=(10,5))
    counts.plot(kind="bar")
    plt.ylabel("Frequency")
    plt.title(f"{model_label}: Significant by County")
    plt.xticks(rotation=75, ha="right")
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / fname, dpi=160)
    plt.close()


# =========================
# MAIN
# =========================
def main():
    ensure_dirs()

    # Ingest & clean
    print_header("Building Air Panel (annual county PM10)")
    air = build_air_panel()
    air.to_csv(OUTPUT_DIR / "air_annual_pm10.csv", index=False)
    print(air.head())

    print_header("Loading Hospital Data (COPD 40+)")
    hosp = load_hospital()
    hosp.to_csv(OUTPUT_DIR / "hospital_clean.csv", index=False)
    print(hosp.head())

    print_header("Merging, Lagging, FIPS")
    panel = build_panel()
    print(panel.head())

    # Recommended: Panel Poisson GLM with FE + Distributed Lags
    model = fit_panel_poisson(panel, lags=PM_LAGS)

    # Placebo test (future PM10 should be null)
    _ = placebo_future_lag(panel, years_ahead=1)

    # Optional replication of per-county OLS + BH (like your R analysis)
    if RUN_PER_COUNTY_OLS_BH:
        print_header("Replicating Per-County OLS + BH (Comparison Only)")
        counties = sorted(panel["County"].dropna().unique().tolist())
        raw = per_county_ols(panel, counties=counties, gaps=range(1,10))

        # Separate families (Count vs Count/Population), apply BH within each
        raw_count = raw[raw["Model"] == "Count"].copy()
        raw_rate  = raw[raw["Model"] == "Count/Population"].copy()

        bh_count = benjamini_hochberg(raw_count, p_col="PValue", alpha=0.05)
        bh_rate  = benjamini_hochberg(raw_rate,  p_col="PValue", alpha=0.05)

        bh_count.to_csv(OUTPUT_DIR / "per_county_ols_bh_count.csv", index=False)
        bh_rate.to_csv(OUTPUT_DIR / "per_county_ols_bh_rate.csv", index=False)

        # Significant only
        sig_count = bh_count[bh_count["significant"]]
        sig_rate  = bh_rate[bh_rate["significant"]]

        # Barplots like your R code
        if not sig_count.empty:
            plot_gap_hist(sig_count, "Count", "sig_gap_hist_count.png")
            plot_county_freq(sig_count, "Count", "sig_by_county_count.png")
        if not sig_rate.empty:
            plot_gap_hist(sig_rate, "Count/Population", "sig_gap_hist_rate.png")
            plot_county_freq(sig_rate, "Count/Population", "sig_by_county_rate.png")

        print("Per-county OLS + BH finished. See outputs/ for CSVs and PNGs.")

    # (Optional) Mapping block placeholder
    if MAKE_MAP:
        print_header("Map Rendering Skipped (set MAKE_MAP=True and provide shapefile to enable).")

    print_header("DONE")
    print(f"Outputs written to: {OUTPUT_DIR.resolve()}")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print("\n[ERROR]", e)
        sys.exit(1)


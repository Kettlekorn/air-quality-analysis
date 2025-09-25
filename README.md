# Air Quality Analysis Across California Counties - Python & R Implementation

Author: Brandon Wu  

---

## Overview
In this project I examine the relationship between ambient air pollution (PM10 concentrations) and hospitalization rates for Chronic Obstructive Pulmonary Disease (COPD) across California counties. This work replicates and extends my original R analysis by implementing the full workflow in Python, including data cleaning, merging, lag construction, regression modeling, and multiple-testing adjustment.

The Python version emphasizes panel data methods and count models, while still allowing replication of the per-county OLS with Benjamini–Hochberg corrections from the R project. For transparency, I’ve included the original R Markdown (`AQICOPDProject.Rmd`) as an archival artifact. Supporting/knitted files are intentionally omitted to keep the repository lean; the R version covered a smaller county subset as an early prototype, and the Python pipeline generalizes it with a unified panel count model
---

## Data Sources
- **Air Quality (daily):** `daily_81102_YYYY.csv` files (1990–2014), each containing daily air quality measures.  
  - Required fields: `State.Name`, `Date.Local`, `County.Name`.  
  - Preferred pollutant column: `Arithmetic.Mean` for PM10 concentration.  
  - If unavailable, the pipeline falls back to the `AQI` field.
- **Hospitalizations:** `CAHospitalizations.csv` with counts and populations for COPD/Asthma among adults aged 40+.  
  - Fields cleaned and used: `Year`, `County`, `Count_ICD9/10`, `Population_ICD9/10`.  
- **(Optional) Crosswalk:** `ca_county_fips.csv` with `County, fips` for mapping and consistency checks.

---

## Methodology

### 1. Data Preparation
- **Air quality:**  
  - Filter to California records.  
  - Derive `Year` from `Date.Local`.  
  - Compute annual mean PM10 by `County, Year`.  
- **Hospitalizations:**  
  - Restrict to COPD/Asthma cases in older adults (age 40+).  
  - Remove non-numeric entries from ICD9/ICD10 fields.  
  - Construct `OverallCount` and `OverallPopulation`.  
  - Exclude `STATEWIDE` entries.  
- **Merge:**  
  - Join by `County, Year`.  
  - Generate lagged PM10 exposures (`PM10_lag1`, `PM10_lag2`, etc.).  
  - Compute `log_pop = log(OverallPopulation)` for offset modeling.

### 2. Regression Models
- **Panel Poisson GLM (recommended):**  
  - Outcome: COPD hospitalization counts.  
  - Exposure: lagged PM10 levels.  
  - Controls: county and year fixed effects.  
  - Offset: log of population.  
  - Inference: cluster-robust standard errors at the county level.  
  - Provides coefficient estimates and a joint Wald test for lag terms.  
- **Placebo check:**  
  - Uses *future* PM10 as a predictor to confirm that spurious associations are not driving results.  
- **Optional replication of R workflow:**  
  - Fits OLS separately for each county across lag values 1–9.  
  - Evaluates significance for both count and crude rate models.  
  - Applies the Benjamini–Hochberg procedure to control false discovery rate.  
  - Produces barplots of significant associations by lag period and county.

---

## Outputs
All outputs are written to the `./outputs` directory:
- **Tables**  
  - `air_annual_pm10.csv` — county–year PM10 means  
  - `hospital_clean.csv` — cleaned COPD counts and populations  
  - `merged_panel.csv` — merged dataset with lags  
  - `panel_poisson_coefs.csv` — coefficients, SEs, p-values  
  - `panel_poisson_summary.txt` — full regression summary  
  - `panel_poisson_lag_wald_test.txt` — joint lag significance test  
  - `per_county_ols_bh_count.csv` / `per_county_ols_bh_rate.csv` (if replication enabled)
- **Figures**  
  - Histograms of significant lag results by gap period  
  - Barplots of significant results by county  

---

## Usage

### Requirements
- Python 3.9+ recommended

## Data acquisition (as of 2025)

This project does not auto-download data. Place the following files in `./data/` before running the script.

### 1) Air quality (PM10, parameter 81102)
- Primary source (pre-generated annual daily files):
  https://aqs.epa.gov/aqsweb/airdata/download_files.html
  • Download the “Daily” files for parameter **81102** (PM10), one archive per year:
    `daily_81102_1990.zip`, …, `daily_81102_2014.zip`
  • Unzip each to obtain CSVs (e.g., `daily_81102_1990.csv`).
  • Expected columns: `State.Name`, `Date.Local`, `County.Name`, and preferably `Arithmetic.Mean` (PM10). If `Arithmetic.Mean` is absent, the pipeline falls back to `AQI`.

- Alternative access (interactive tool):
  https://www.epa.gov/outdoor-air-quality-data/download-daily-data

### 2) Hospitalizations (COPD/Asthma, age 40+)
- California Health & Human Services (CHHS) Open Data Portal:
  https://data.chhs.ca.gov/dataset/rates-of-preventable-hospitalizations-for-selected-medical-conditions-by-county
  • Save the county-level table as `CAHospitalizations.csv`.
  • The script filters `PQIDescription == "COPD or Asthma in Older Adults (Age 40+)"` and uses:
    `Year`, `County`, `Count_ICD10`, `Population_ICD10`, `Count_ICD9`, `Population_ICD9`.

### 3) Measure definition (context/reference)
- AHRQ Prevention Quality Indicator 05 (COPD or Asthma in Older Adults) technical specifications:
  https://qualityindicators.ahrq.gov/measures/PQI_TechSpec
  • Use the most recent PDF linked on that page for coding/measure details.

### Suggested directory layout
./data/
├── daily_81102_1990.csv
├── daily_81102_1991.csv
├── …
├── daily_81102_2014.csv
├── CAHospitalizations.csv
└── (optional) ca_county_fips.csv

### Installation
pip install pandas numpy statsmodels linearmodels matplotlib geopandas shapely pyproj us

### Running the analysis
python copd_pm10_panel.py

> Place all input CSVs in ./data/. Outputs will be written to ./outputs/.

---

## Notes (step by step)

1. Air data prep: I load the daily `daily_81102_YYYY.csv` files, keep only California rows, parse `Date.Local` into `Year`, and compute a county–year mean for PM10. I prefer `Arithmetic.Mean`; if it is missing I fall back to `AQI`.

2. Hospitalization cleaning: I filter to “COPD or Asthma in Older Adults (Age 40+)”. I strip non-digits from ICD9/ICD10 fields, build `OverallCount` and `OverallPopulation`, drop `STATEWIDE`, and keep `Year` and `County`.

3. Merge and lags: I merge on `County, Year`, then create lagged PM10 exposures (e.g., `PM10_lag1`, `PM10_lag2`). I also compute `log_pop = log(OverallPopulation)` for use as an offset.

4. Main model (panel Poisson GLM): I model counts with county and year fixed effects and an offset for population. Outcome = `OverallCount`, exposure = lagged PM10, controls = `C(County)` and `C(Year)`, offset = `log_pop`. Inference uses cluster-robust SEs by county. I also run a joint Wald test on the block of lag coefficients.

5. Placebo check: I repeat the model using future PM10 (a lead) to confirm it does not predict current hospitalizations.

6. Optional replication of my R workflow: If enabled, I fit per-county OLS across lag values (1–9) for both counts and crude rates, apply Benjamini–Hochberg FDR, and export barplots for significant results.

7. Outputs I review: 
   - `panel_poisson_summary.txt` and `panel_poisson_coefs.csv` for effect sizes and p-values  
   - `panel_poisson_lag_wald_test.txt` for the joint significance of lags  

8. Interpretation notes: PM10 concentration is more interpretable than AQI. Poisson with an offset is appropriate for counts. County fixed effects absorb time-invariant local differences; year fixed effects absorb statewide shocks (policy, coding, wildfire years). The placebo helps guard against spurious correlations.

---

## Repository contents and R implementation

This repository contains the Python pipeline and the original R implementation, along with the supplementary outputs generated during analysis. The folder already includes the processed CSVs (e.g., `merged_panel.csv`, BH-adjusted p-value tables) and figures used in the report, so the project is fully browsable without rerunning the code.

### R implementation (scope and purpose)
- The R notebook in `r/` documents the early version of the analysis and final narrative report. It focuses on per-county OLS regressions with Benjamini–Hochberg correction, consistent with the study design described in the report.  
- **Scope note.** The R implementation was intentionally limited to a smaller geographic scope

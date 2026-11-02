# Air Pollution, Green Space, and Neurodevelopmental Outcomes

R scripts for studying associations of wildfire PM2.5, ambient PM2.5, NO2, ozone, and green space with childhood neurodevelopmental outcomes using Cox models and exposure windows.

## Analysis files

| Script | Main analysis |
| --- | --- |
| `Wildfire_ASD_JAMA_Pediatrics_code.R` | Wildfire and background PM2.5; ASD and ADHD; 33 prenatal weeks |
| `PM_NDD_ES&T.txt` | Ambient PM2.5 and temperature; seven outcomes; 12 preconception weeks, 37 prenatal weeks and 48 postnatal months |
| `Gaseous_NDD_complete_annual.R` | NO2 and ozone; eight outcomes, including any neurodevelopmental disorder; 12 preconception weeks, 37 prenatal weeks and four postnatal years |
| `GreenSpace_NDD_complete.R` | NDVI for preconception, pregnancy and postnatal windows; separate and mutually adjusted models, subgroups, VIF, tree-canopy, age-restricted and IPCW analyses |

## Requirements and inputs

Use R with compatible versions of `dlnm` and `survival`; `splines` is included with R.

```r
install.packages(c("dlnm", "survival"))
```

Inputs must be prepared RDS datasets with one row per eligible child, validated outcomes, linked exposures and required covariates. Cohort construction, eligibility exclusions and exposure linkage precede these scripts. Follow-up durations are in days from birth. Use `Deliv_Date` and `season_conception`.

- **Wildfire:** configure `input_file`. The full workflow requires weekly `gsPM_1`–`gsPM_37` and `gnsPM_1`–`gnsPM_37`, ASD/ADHD endpoints and the fields listed in Section 0. 
- **Ambient PM2.5:** configure outcome-specific `CFG$input_files` and verify `CFG$column_maps`. Supply actual age-aligned monthly PM2.5 and temperature, in wide format or through `monthly_long_file`. RH remains four annual values. 
- **NO2/O3:** configure `CFG$input_files`, or load `dat`, `dat_26`, `dat_30` and `dat_32`. Shorter-window sensitivities require the corresponding eligible cohorts or a prepared cohort before prenatal-completeness selection. 
- **Green space:** configure the input path and field names in Step 1 of `GreenSpace_NDD_complete.R`. Provide one NDVI summary per child for each exposure window (preconception, pregnancy and postnatal), plus outcome indicators, follow-up, covariates, and diagnosis age for the age-restricted analyses. 
The scripts define their variable lists and model settings. PM2.5 and gaseous inputs require unique `patient_birth_ID` values.

## Running an analysis

1. Edit paths and settings inside wildfire Section 0 or the Step 1 configuration of the PM2.5, gaseous and green-space scripts.
2. Choose a new or empty output directory.
3. Run each analysis in a separate R session:

```bash
Rscript Wildfire_ASD_JAMA_Pediatrics_code.R
Rscript PM_NDD_ES&T.txt
Rscript Gaseous_NDD_complete_annual.R
Rscript GreenSpace_NDD_complete.R
```

Alternatively, use `source("filename")` in R. Use this approach when gaseous input objects are already loaded. Relative paths use the working directory. Sourcing recreates `CFG`, so edit configuration entries inside the script. PM2.5 and gaseous `run_*` switches control optional analyses.

## Main settings and interpretation

- **Wildfire:** natural exposure spline with 3 df, lag spline with 2 df; background PM2.5 is linear with 2 lag df. HRs compare one calculated IQR with zero.
- **Ambient PM2.5:** main early-exposure lag df are 6/4/4 for PM2.5/temperature/RH; early subgroups and postnatal models use 4/4/4. Main exposure functions are linear. PM2.5 HRs compare 19 versus 9 µg/m³.
- **NO2/O3:** linear exposure functions and 4 lag df; HRs per 10 ppb, with co-pollutant adjustment.
- **Green space:** NDVI models use a fixed 0.12 increase; each exposure window is modeled separately and then mutually adjusted for the other windows. The script also calculates VIF across the three NDVI windows, subgroup models, optional tree-canopy models, diagnosis-age restrictions at 3 and 8 years, and stabilized censoring weights capped at the 99th percentile.
- Subgroup, sensitivity and IPCW analyses are included. 
## Outputs and checks

Outputs include HR tables with 95% confidence intervals, descriptive summaries, subgroup/sensitivity results, exposure mappings, IPCW diagnostics and plots. All scripts save `sessionInfo.txt`.

For PM2.5 and NO2/O3, start with `all_HRs.csv` and `analysis_status.csv`; review every failure or warning. For wildfire, start with `overall_table.csv` and `sensitivity_results.csv` and inspect the run log. The green-space script saves `all_HRs.csv`, NDVI summaries, VIF results and settings.

Check sample sizes, event counts, exposure ordering, convergence and weight diagnostics before interpreting results. Retain script versions, settings and input provenance; keep participant-level data within the authorized environment.

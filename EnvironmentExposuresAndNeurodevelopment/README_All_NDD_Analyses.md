# Air Pollution and Neurodevelopmental Outcomes

R scripts for studying associations of wildfire PM2.5, ambient PM2.5, NO2 and ozone with childhood neurodevelopmental outcomes using Cox models and distributed lag exposure terms.

## Analysis files

| Script | Main analysis |
| --- | --- |
| `Wildfire_ASD_JAMA_Pediatrics_code.R` | Wildfire and background PM2.5; ASD and ADHD; 33 prenatal weeks |
| `PM_NDD_ES&T.R` | Ambient PM2.5 and temperature; seven outcomes; 12 preconception weeks, 37 prenatal weeks and 48 postnatal months |
| `Gaseous_NDD_complete_annual.R` | NO2 and ozone; seven outcomes; 12 preconception weeks, 37 prenatal weeks and four postnatal years |

## Requirements and inputs

Use R with compatible versions of `dlnm`, `survival` and `splines`. The latter is included with R.

```r
install.packages(c("dlnm", "survival"))
```

Inputs must be prepared RDS datasets with one row per eligible child, validated outcomes, linked exposures and required covariates. Cohort construction, eligibility exclusions and exposure linkage precede these scripts. Follow-up durations are in days from birth. Use `Deliv_Date` and `season_conception`; verify dates, units and exposure chronology.

- **Wildfire:** configure `input_file`. The full workflow requires weekly `gsPM_1`–`gsPM_37` and `gnsPM_1`–`gnsPM_37`, ASD/ADHD endpoints and the fields listed in Section 0. Exposure columns retain descending suffix order; check `exposure_map.csv` before interpreting gestational timing.
- **Ambient PM2.5:** configure outcome-specific `CFG$input_files` and verify `CFG$column_maps`. Supply actual age-aligned monthly PM2.5 and temperature, in wide format or through `monthly_long_file`. RH remains four annual values. Repeating annual averages across months is not valid monthly input.
- **NO2/O3:** configure `CFG$input_files`, or load `dat`, `dat_26`, `dat_30` and `dat_32`. Shorter-window sensitivities require the corresponding eligible cohorts or a prepared cohort before prenatal-completeness selection. Supply the verified administrative end and exit-reason coding for disenrollment IPCW.

The scripts define their complete variable lists and model settings. PM2.5 and gaseous inputs require unique `patient_birth_ID` values.

## Running an analysis

1. Edit paths and settings inside wildfire Section 0 or the PM2.5/gaseous Step 1 configuration.
2. Choose a new or empty output directory.
3. Run each analysis in a separate R session:

```bash
Rscript Wildfire_ASD_JAMA_Pediatrics_code.R
Rscript PM_NDD_ES&t.R
Rscript Gaseous_NDD_complete_annual.R
```

Alternatively, use `source("filename")` in R. Use this approach when gaseous input objects are already loaded. Relative paths use the working directory. Sourcing recreates `CFG`, so edit configuration entries inside the script. PM2.5 and gaseous `run_*` switches control optional analyses.

## Main settings and interpretation

- **Wildfire:** natural exposure spline with 3 df, lag spline with 2 df; background PM2.5 is linear with 2 lag df. HRs compare one calculated IQR with zero.
- **Ambient PM2.5:** main early-exposure lag df are 6/4/4 for PM2.5/temperature/RH; early subgroups and postnatal models use 4/4/4. Main exposure functions are linear. PM2.5 HRs compare 19 versus 9 µg/m³.
- **NO2/O3:** linear exposure functions and 4 lag df; HRs per 10 ppb, with co-pollutant adjustment.
- Subgroup, sensitivity and IPCW analyses are included. Both postnatal scripts default to follow-up from birth with fixed exposure histories, and diagnoses after the timeframe of exposure.


## Outputs and checks

Outputs include HR tables with 95% confidence intervals, descriptive summaries, subgroup/sensitivity results, exposure mappings, IPCW diagnostics and plots. All scripts save `sessionInfo.txt`.

For PM2.5 and NO2/O3, start with `all_HRs.csv` and `analysis_status.csv`; review every failure or warning. For wildfire, start with `overall_table.csv` and `sensitivity_results.csv` and inspect the run log.

Check sample sizes, event counts, exposure ordering, convergence and weight diagnostics before interpreting results. Retain script versions, settings and input provenance; keep participant-level data within the authorized environment.

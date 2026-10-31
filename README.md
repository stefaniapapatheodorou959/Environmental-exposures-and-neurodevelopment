# Wildfire_ASD_JAMA_Pediatrics
Code and README files
Prenatal Wildfire PM2.5 and ASD/ADHD
R code for evaluating associations between prenatal wildfire PM2.5 exposure and autism spectrum disorder (ASD) and attention-deficit/hyperactivity disorder (ADHD). The analysis uses Cox proportional-hazards models with distributed lag nonlinear exposure terms and includes descriptive, subgroup, and sensitivity analyses.

Related publication: Papatheodorou et al., JAMA Pediatrics. doi:10.1001/jamapediatrics.2026.5045.
Software requirements
The analysis requires R and the following packages:
- dlnm: exposure and lag-response functions and predictions.
- survival: survival models and Kaplan–Meier estimates.
- splines: natural spline functions; included with R.
Install the external packages in a compatible R environment:
install.packages(c("dlnm", "survival"))
The script records the R and package versions in sessionInfo.txt. Use a consistent software environment when comparing results across runs. The AIC model grid and censoring-weight calculations can require substantial computation and memory.
Input data
The input is an RDS data frame containing one row per eligible mother–child observation. Cohort eligibility, mother–child linkage, prenatal exposure assignment, outcome ascertainment, and follow-up endpoints must be established before running the analysis.

Data component	Required information
Outcomes	Event_1_MAIN for ASD and Event_3_MAIN for ADHD, coded 0/1
Follow-up	FU_Event_1_Correct and FU_Event_3_Correct, measured in days from birth to the corresponding outcome or follow-up endpoint
Weekly exposures	Numeric wildfire concentrations gsPM_1–`gsPM_37` and background concentrations gnsPM_1–`gnsPM_37`, in µg/m³
Individual characteristics	Maternal age, race/ethnicity, conception season, binary parity proxy, prenatal folate dispensing, substance use indicators, and child sex
Contextual characteristics	Income, education, housing, population composition, population density, county BMI, hospital distance, urbanicity, and pregnancy temperature
Geography and clustering	birth_yr, FIPS, ZIP_Index, state, and a maternal grouping variable mother_cluster
Descriptive measures	Pregnancy-average wildfire and background PM2.5 and Infant_length_enrollment
Calendar dates	Birth date and outcome-specific diagnosis dates for the ICD-period analysis
Follow-up reasons	asd_exit_reason and adhd_exit_reason for censoring-weight estimation


Step 0 lists the source column names. Use consistent names throughout the script. Binary fields must contain 0, 1, or missing values. Infant_Gender_Male uses 0 for female and 1 for male. Calendar dates must have a verified date representation; numeric dates require conversion using their correct origin.

race_factor accepts White, Black, Hispanic, Asian/Pacific, and Other/Unknown. If this field is absent, Step 1 derives it from the five mutually exclusive race/ethnicity indicators named in that section.

Follow-up reasons must use these exact values:

Value	Meaning
outcome	The outcome occurred at the end of observed follow-up
loss_to_followup	Observable follow-up ended through loss of coverage or follow-up
administrative	Follow-up reached the prespecified study end
death	Death ended follow-up before the outcome

Running the analysis
1. Place the script and input data in the approved analysis environment.
2. Set the working directory and edit the paths in Step 0:
   input_file <- "data/data_wild2_excl.rds"
   output_dir <- "results/wildfire_analysis"
3. Use a new or empty output directory and run the script from beginning to end:
   source("Wildfire_ASD_ADHD_FINAL_CODE.txt")
The same script can be run from a terminal:
Rscript Wildfire_ASD_ADHD_FINAL_CODE.txt
Relative paths are resolved from the working directory. The script creates the output directory and stops if it already contains files. All analysis sections run sequentially.
Analysis workflow
Steps	Purpose
0–5	    Validate inputs; define race/ethnicity, regions, exposure histories, IQRs, and outcomes
6–8	    Fit primary ASD and ADHD models and extract cumulative hazard ratios
9–14	Fit sex-, race/ethnicity-, and region-specific models
15–19	Produce overall and subgroup plots and summaries
20   	Summarize participant characteristics, follow-up, and age-8 incidence
21  	Compare exposure and lag spline degrees of freedom using AIC
22–25	Evaluate 28- and 37-week windows, diagnoses after ages 3 and 5, linear exposure response, and minimal adjustment
26	    Apply inverse-probability-of-censoring weights (IPCW)
27–29	Evaluate maternal clustering, ICD diagnosis periods, and E-values
30	    Save combined sensitivity results, spline settings, and software versions


The primary model covers 33 prenatal weeks. Wildfire PM2.5 uses a natural spline with 3 degrees of freedom for exposure response and 2 degrees of freedom for the lag response; background PM2.5 uses a linear exposure term and a 2-df lag spline. Baseline hazards are stratified by birth year and county, with ZIP-clustered robust standard errors. The full adjustment set is defined in full_covariates.

Each outcome model uses complete cases for its required variables and weekly exposures. Sample sizes may therefore differ across outcomes and analyses. The exposure matrices use descending column suffixes: gsPM_33 corresponds to modeled lag 1 in the primary model. Confirm the correspondence between these suffixes and gestational timing before interpreting lag-response plots.

Cumulative hazard ratios compare one IQR with zero exposure at every included week. The IQR is calculated from the pooled weekly concentrations. Subgroup models use the primary exposure contrast; alternative exposure windows have their own IQRs. Diagnosis-age sensitivity analyses exclude early cases while retaining birth as the time origin for included observations.

IPCW uses outcome-specific models of loss to follow-up and time-varying stabilized weights evaluated immediately before outcome event times. The denominator includes prenatal exposures, the full adjustment set, child sex, and poverty; both censoring models use birth-year/county strata. Administrative endpoints and death are not treated as loss-to-follow-up events. Interpretation requires conditional independent censoring, adequate positivity, and appropriate censoring models. Robust standard errors condition on the estimated weights.
Output files
All files are written to output_dir.
Files: Contents
cohort_counts.csv, exposure_contrasts.csv, exposure_map.csv: Cohort counts, exposure comparisons, and weekly column ordering
overall_table.csv, sex_table.csv, race_table.csv, region_table.csv, overall_and_subgroup_results.csv: Primary and subgroup estimates
table1_long.csv, table1_threshold.csv, followup_metrics.csv, incidence_age8.csv: Descriptive summaries; age-8 estimates require sufficient observed follow-up
AIC_ASD.csv, AIC_ADHD.csv: Spline model comparisons and any recorded fitting errors
alternative_windows.csv, diagnosis_age_sensitivity.csv, linear_sensitivity.csv, minimal_model.csv: Exposure-window, diagnosis-age, and model-specification sensitivity results
main_vs_ipcw_comparison.csv, IPCW_weights_*.csv, IPCW_endpoints_*.csv, IPCW_settings_*.txt:	Matched unweighted/weighted estimates and censoring diagnostics
maternal_clustering.csv, ICD_transition.csv, E_values.csv, sensitivity_results.csv:	Additional sensitivity analyses and combined estimates
ER_*.png, LAG_*.png:	Cumulative exposure-response and weekly lag-response plots
basis_settings.txt, sessionInfo.txt: Spline definitions and software environment


Model tables report the analysis sample size, outcome events, person-periods, exposure contrast, hazard ratio (HR), and 95% confidence limits (LCL, UCL). Person-period counts can exceed participant counts when follow-up is split. IPCW weight diagnostics summarize risk-set memberships rather than unique children. The script does not generate the national annual exposure summary in Table 2 of the publication.
Reviewing a run
Check cohort counts, model-specific sample sizes, event totals, exposure ordering, convergence messages, and censoring-weight diagnostics before interpreting results. Resolve missing fields, inconsistent endpoint dates, nonfinite estimates, or singular models at their source. A successful script run alone does not establish reproduction of published estimates.

Keep participant-level data and grouping identifiers within the authorized data environment. The script exports aggregate results and does not save participant records or fitted model objects. Apply the relevant disclosure requirements before sharing outputs.
# Environmental Exposures and Neurodevelopment

This repository contains  R analysis sample code related to environmental exposures during pregnancy and childhood neurodevelopmental outcomes. The analyses use Medicaid mother–child cohort data and examine air pollution and green space.

## Analyses

The repository includes code for:

- Prenatal wildfire and background PM2.5 in relation to ASD and ADHD (PMCID: PMC13711229).
- Ambient PM2.5 with monthly postnatal exposure analyses (PMCID: PMC13674074).
- NO2 and ozone with annual postnatal exposure analyses (PMID: 42480828).
- NDVI during preconception, pregnancy, and postnatal periods (PMCID: PMC12380237).

See the individual scripts and their configuration sections for the outcomes, exposure windows, model settings, and required variables.

## Data and use

The analytic datasets are not included. Scripts require access to the appropriately prepared cohort data and any linked exposure files in an authorized research environment. Do not upload participant-level or otherwise restricted data to this repository.

These scripts do not perform every upstream data-preparation step. Cohort construction, eligibility exclusions, outcome ascertainment, and exposure linkage should be completed and checked before running the analyses. Review each script’s configuration and confirm that the variable names, dates, units, exposure order, and censoring definitions match your data.

## Requirements

Use a compatible R installation. Depending on the script, install:

```r
install.packages(c("survival", "dlnm"))
```

The `splines` package is included with R.

## Running the scripts

1. Download or clone the repository.
2. Open the relevant script and update the input paths, output directory, and variable mappings in its configuration section.
3. Run the script from R or from a terminal with `Rscript`.

For example:

```bash
Rscript path/to/analysis_script.R
```

Outputs and checks are described in each script. Review sample sizes, event counts, warnings, model convergence, exposure coding, and IPCW diagnostics before interpreting results.

## Reproducibility

The scripts support the analyses but do not, by themselves, guarantee reproduction of published estimates. Exact replication depends on the original cohort, exposure-linkage files, variable definitions, exclusions, software versions, and analysis settings.

## Funding

This work was supported by the National Institute of Environmental Health Sciences, National Institutes of Health (R01 ES034038).

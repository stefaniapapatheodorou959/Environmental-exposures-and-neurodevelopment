############################################################
## PRENATAL WILDFIRE PM2.5, ASD AND ADHD
## Papatheodorou et al. JAMA Pediatrics.
## DOI: 10.1001/jamapediatrics.2026.5045
## Overall, subgroup and sensitivity analyses.
############################################################

library(dlnm)
library(survival)
library(splines)

############################################################
## 0) READ DATA AND CHECK REQUIRED FIELDS
############################################################
## One row per eligible mother-child observation. Cohort assembly,
## exposure linkage and outcome ascertainment precede this script.
## Edit these two paths for the approved analysis environment.

input_file <- "data/data_wild2_excl.rds"
output_dir <- "results/wildfire_analysis"

data_wild2_excl <- readRDS(input_file)
if (!is.data.frame(data_wild2_excl) || !nrow(data_wild2_excl) ||
    anyDuplicated(names(data_wild2_excl))) stop("Invalid input data frame.")
if (dir.exists(output_dir) &&
    length(list.files(output_dir, all.files = TRUE, no.. = TRUE))) {
  stop("Use an empty output directory.")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

require_columns <- function(data, fields) {
  missing <- setdiff(fields, names(data))
  if (length(missing)) stop("Missing fields: ", paste(missing, collapse = ", "))
}
required_fields <- c(
  "FU_Event_1_Correct", "Event_1_MAIN", "FU_Event_3_Correct", "Event_3_MAIN",
  "birth_yr", "FIPS", "ZIP_Index", "Age", "season_conception", "UA",
  "repeated_pregnancy", "pop_density_30sec_NASA", "BlkPct", "HspPct",
  "nearest_hospital_km", "temp", "Infant_Gender_Male", "pct_owner_occ",
  "mean_bmi", "Alcohol_abuse", "Smoking", "Substance_abuse", "Poor_nutrition",
  "prenatal_folate", "medhouseholdincome", "education", "poverty", "state",
  "mother_cluster", "Deliv_Date", "Event_1_diagdate", "Event_3_diagdate",
  "Infant_length_enrollment", "pregnancy_wildfire_mean", "pregnancy_background_mean",
  "asd_exit_reason", "adhd_exit_reason"
)
require_columns(data_wild2_excl, required_fields)
data_wild2_excl$.row_index <- seq_len(nrow(data_wild2_excl))

binary_fields <- c("repeated_pregnancy", "prenatal_folate", "Smoking",
                   "Alcohol_abuse", "Substance_abuse", "Poor_nutrition",
                   "Infant_Gender_Male", "Event_1_MAIN", "Event_3_MAIN")
for (v in binary_fields) {
  x <- data_wild2_excl[[v]]
  if (is.factor(x)) x <- as.character(x)
  if (is.character(x)) x <- trimws(x)
  if (is.logical(x)) x <- as.integer(x)
  if ((!is.numeric(x) && !is.character(x)) ||
      any(!is.na(x) & !x %in% c(0, 1))) stop("Expected 0/1/NA in ", v)
  data_wild2_excl[[v]] <- as.integer(x)
}
for (v in c("season_conception", "UA")) {
  data_wild2_excl[[v]] <- factor(data_wild2_excl[[v]])
}

############################################################
## 1) DEFINE RACE FACTOR
############################################################

race_levels_all <- c("White", "Black", "Hispanic", "Asian/Pacific", "Other/Unknown")
if (!"race_factor" %in% names(data_wild2_excl)) {
  race_fields <- c("Race_WHITE", "Race_BLACK_OR_AFRICAN_AMERICAN",
                   "Race_HISPANIC_OR_LATINO", "Race_ASIAN_OR_OTHER_PACIFIC",
                   "Race_UNKNOWN_OTHER")
  require_columns(data_wild2_excl, race_fields)
  race_data <- data_wild2_excl[race_fields]
  if (!all(vapply(race_data, function(x) is.numeric(x) &&
                  all(is.na(x) | x %in% c(0, 1)), logical(1)))) {
    stop("Race indicators must be numeric 0/1/NA.")
  }
  complete <- complete.cases(race_data)
  if (any(rowSums(race_data[complete, , drop = FALSE]) != 1)) {
    stop("Resolve overlapping or unclassified race indicators.")
  }
  race <- rep(NA_character_, nrow(data_wild2_excl))
  for (j in seq_along(race_fields)) {
    race[which(complete & race_data[[j]] == 1)] <- race_levels_all[j]
  }
  data_wild2_excl$race_factor <- race
}
race <- as.character(data_wild2_excl$race_factor)
if (any(!is.na(race) & !race %in% race_levels_all)) stop("Unexpected race category.")
data_wild2_excl$race_factor <- factor(race, levels = race_levels_all)
race_levels <- setdiff(race_levels_all, "Other/Unknown")

############################################################
## 2) US REGIONS
############################################################

northeast_states <- c("CT", "ME", "MA", "NH", "RI", "VT", "NJ", "NY", "PA")
midwest_states <- c("IL", "IN", "MI", "OH", "WI", "IA", "KS", "MN", "MO", "NE", "ND", "SD")
south_states <- c("DE", "DC", "FL", "GA", "MD", "NC", "SC", "VA", "WV",
                  "AL", "KY", "MS", "TN", "AR", "LA", "OK", "TX")
west_states <- c("AK", "AZ", "CA", "CO", "HI", "ID", "MT", "NV", "NM", "OR", "UT", "WA", "WY")
region_states <- list(NORTHEAST = northeast_states, MIDWEST = midwest_states,
                      SOUTH = south_states, WEST = west_states)
state <- toupper(trimws(as.character(data_wild2_excl$state)))
if ("STATECODE" %in% names(data_wild2_excl)) {
  fallback <- toupper(trimws(as.character(data_wild2_excl$STATECODE)))
  valid <- unlist(region_states, use.names = FALSE)
  replace <- !state %in% valid & fallback %in% valid
  state[replace] <- fallback[replace]
}
data_wild2_excl$state <- state
data_wild2_excl$region <- NA_character_
for (region in names(region_states)) {
  data_wild2_excl$region[state %in% region_states[[region]]] <- region
}
region_levels <- names(region_states)
data_wild2_excl$region <- factor(data_wild2_excl$region, levels = region_levels)

############################################################
## 3) DEFINE REVERSED EXPOSURE MATRICES ONCE
############################################################
## Select columns by name rather than position. Preserve the original
## reversed order: gsPM_33 is modeled lag 1 in the 33-week analysis.
## Confirm suffix chronology against the exposure data dictionary.

build_exposures <- function(data, lag_max) {
  wild_names <- paste0("gsPM_", lag_max:1)
  bg_names <- paste0("gnsPM_", lag_max:1)
  require_columns(data, c(wild_names, bg_names))
  if (!all(vapply(data[c(wild_names, bg_names)], is.numeric, logical(1)))) {
    stop("Weekly exposure columns must be numeric.")
  }
  wild <- as.matrix(data[wild_names])
  bg <- as.matrix(data[bg_names])
  for (x in list(wild, bg)) {
    if (all(is.na(x)) || any(!is.finite(x) & !is.na(x)) || any(x < 0, na.rm = TRUE)) {
      stop("Invalid weekly exposure concentrations.")
    }
  }
  list(wild = wild, bg = bg)
}
exposures <- build_exposures(data_wild2_excl, 33L)
revexp_s3 <- exposures$wild
revexp_ds3 <- exposures$bg
rm(exposures)

############################################################
## 4) GLOBAL IQR AND EXPOSURE CONTRAST
############################################################

IQR_wild <- IQR(as.vector(revexp_s3), na.rm = TRUE)
IQR_bg <- IQR(as.vector(revexp_ds3), na.rm = TRUE)
if (!is.finite(IQR_wild) || IQR_wild <= 0 || !is.finite(IQR_bg) || IQR_bg <= 0) {
  stop("Exposure IQRs must be finite and positive.")
}
write.csv(data.frame(exposure = c("wildfire", "background"), reference = 0,
                     comparison = c(IQR_wild, IQR_bg)),
          file.path(output_dir, "exposure_contrasts.csv"), row.names = FALSE)
write.csv(data.frame(lag = 1:33, wildfire_column = colnames(revexp_s3),
                     background_column = colnames(revexp_ds3)),
          file.path(output_dir, "exposure_map.csv"), row.names = FALSE)

############################################################
## 5) OUTCOMES AND FOLLOW-UP
############################################################

outcomes <- data.frame(
  name = c("ASD", "ADHD"),
  fu = c("FU_Event_1_Correct", "FU_Event_3_Correct"),
  ev = c("Event_1_MAIN", "Event_3_MAIN"),
  exit = c("asd_exit_reason", "adhd_exit_reason"),
  diagnosis_date = c("Event_1_diagdate", "Event_3_diagdate")
)
for (i in seq_len(nrow(outcomes))) {
  t <- data_wild2_excl[[outcomes$fu[i]]]
  if (!is.numeric(t) || any(!is.na(t) & (!is.finite(t) | t <= 0))) {
    stop("Invalid follow-up in ", outcomes$fu[i])
  }
}
cohort_counts <- data.frame(
  metric = c("Cohort N", "ASD events", "ADHD events"),
  observed = c(nrow(data_wild2_excl), sum(data_wild2_excl$Event_1_MAIN, na.rm = TRUE),
               sum(data_wild2_excl$Event_3_MAIN, na.rm = TRUE)),
  published = c(893212, 8357, 22302)
)
print(cohort_counts)
write.csv(cohort_counts, file.path(output_dir, "cohort_counts.csv"), row.names = FALSE)

############################################################
## 6) EXTRACT THE EXACT CUMULATIVE HR AND FORMAT RESULTS
############################################################
extract_cumHR <- function(predLG) {
  c(HR = unname(predLG$allRRfit[1]), LCL = unname(predLG$allRRlow[1]),
    UCL = unname(predLG$allRRhigh[1]))
}
result_row <- function(fit, outcome, analysis) {
  data.frame(outcome = outcome, analysis = analysis, N = fit$n,
    events = fit$events, person_periods = fit$periods,
    exposure = c("wildfire", "background"), reference = 0,
    IQR = c(fit$IQR_wild, fit$IQR_bg),
    HR = c(fit$cumHR["HR"], fit$bgHR["HR"]),
    LCL = c(fit$cumHR["LCL"], fit$bgHR["LCL"]),
    UCL = c(fit$cumHR["UCL"], fit$bgHR["UCL"]), row.names = NULL)
}
save_table <- function(rows, filename) {
  table <- do.call(rbind, rows)
  write.csv(table, file.path(output_dir, filename), row.names = FALSE)
  table
}

############################################################
## 7) HELPER: OVERALL MODEL
############################################################

full_covariates <- c(
  "Age", "race_factor", "repeated_pregnancy", "season_conception",
  "prenatal_folate", "Substance_abuse", "Smoking", "Alcohol_abuse", "Poor_nutrition",
  "medhouseholdincome", "education", "pct_owner_occ", "BlkPct", "HspPct",
  "pop_density_30sec_NASA", "mean_bmi", "nearest_hospital_km", "UA", "temp"
)

analytic_sample <- function(data, exp_wild_mat, exp_bg_mat, fu_var, ev_var,
                            covariates = full_covariates, cluster_var = "ZIP_Index") {
  needed <- unique(c(fu_var, ev_var, "birth_yr", "FIPS", cluster_var, covariates))
  require_columns(data, c(needed, ".row_index"))
  if (anyDuplicated(data$.row_index) || anyNA(data$.row_index) ||
      any(data$.row_index < 1 | data$.row_index > nrow(exp_wild_mat)) ||
      nrow(exp_wild_mat) != nrow(exp_bg_mat)) stop("Exposure rows do not align.")
  keep <- complete.cases(data[needed]) &
    complete.cases(exp_wild_mat[data$.row_index, , drop = FALSE]) &
    complete.cases(exp_bg_mat[data$.row_index, , drop = FALSE])
  dat <- droplevels(data[keep, , drop = FALSE])
  numeric_fields <- needed[vapply(dat[needed], is.numeric, logical(1))]
  if (any(vapply(dat[numeric_fields], function(x) any(!is.finite(x)), logical(1)))) {
    stop("Nonfinite model covariate.")
  }
  dat
}

fit_overall <- function(fu_var, ev_var, data, exp_wild_mat, exp_bg_mat,
                        IQR_wild, IQR_bg, covariates = full_covariates,
                        exposure_df = 3L, lag_df = 2L, linear = FALSE,
                        cluster_var = "ZIP_Index", intervals = NULL, predictions = TRUE) {
  dat <- analytic_sample(data, exp_wild_mat, exp_bg_mat, fu_var, ev_var,
                         covariates, cluster_var)
  if (is.null(intervals)) {
    dat$.start <- 0
    dat$.stop <- dat[[fu_var]]
    dat$.event <- dat[[ev_var]]
    dat$.weight <- 1
  } else {
    require_columns(intervals, c("row_index", "tstart", "tstop", "event", "weight"))
    it <- intervals[intervals$row_index %in% dat$.row_index, , drop = FALSE]
    dat <- dat[match(it$row_index, dat$.row_index), , drop = FALSE]
    dat$.start <- it$tstart; dat$.stop <- it$tstop
    dat$.event <- it$event; dat$.weight <- it$weight
  }
  if (!nrow(dat) || !sum(dat$.event)) stop("No eligible observations or outcome events.")
  if (any(!is.finite(as.matrix(dat[c(".start", ".stop", ".event", ".weight")]))) ||
      any(dat$.start < 0 | dat$.stop <= dat$.start | dat$.weight <= 0) ||
      any(!dat$.event %in% c(0, 1))) stop("Invalid survival intervals.")

  ## Construct each basis once per child, then expand for split follow-up.
  rows <- unique(dat$.row_index)
  lag_max <- ncol(exp_wild_mat)
  if (lag_max != ncol(exp_bg_mat)) stop("Exposure windows differ.")
  cb_wild <- crossbasis(
    exp_wild_mat[rows, , drop = FALSE], lag = c(1, lag_max),
    argvar = if (linear) list(fun = "lin") else list(fun = "ns", df = exposure_df),
    arglag = list(fun = "ns", df = lag_df)
  )
  cb_bg <- crossbasis(
    exp_bg_mat[rows, , drop = FALSE], lag = c(1, lag_max),
    argvar = list(fun = "lin"), arglag = list(fun = "ns", df = 2L)
  )
  expand_basis <- function(cb, index) {
    attributes_before <- attributes(cb)
    ans <- cb[index, , drop = FALSE]
    for (nm in setdiff(names(attributes_before), c("dim", "dimnames"))) {
      attr(ans, nm) <- attributes_before[[nm]]
    }
    ans
  }
  index <- match(dat$.row_index, rows)
  cb_wild <- expand_basis(cb_wild, index)
  cb_bg <- expand_basis(cb_bg, index)

  formula <- as.formula(paste(
    "Surv(.start, .stop, .event) ~ cb_wild + cb_bg +",
    "strata(birth_yr) + strata(FIPS) +", paste(covariates, collapse = " + ")
  ))
  mod <- coxph(formula, data = dat, weights = .weight,
               cluster = dat[[cluster_var]], robust = TRUE, ties = "efron",
               na.action = na.fail, singular.ok = FALSE,
               x = FALSE, y = FALSE, model = FALSE)
  if (any(!is.finite(coef(mod)))) stop("Nonfinite model coefficients.")
  if (!predictions) return(list(AIC = AIC(mod)))

  predLG <- crosspred(cb_wild, model = mod, cen = 0, at = IQR_wild)
  predBG <- crosspred(cb_bg, model = mod, cen = 0, at = IQR_bg)
  predER <- crosspred(cb_wild, model = mod, cen = 0,
                      at = seq(0, max(exp_wild_mat, na.rm = TRUE), length.out = 100))
  cumHR <- extract_cumHR(predLG)
  bgHR <- extract_cumHR(predBG)
  if (any(!is.finite(c(cumHR, bgHR)))) stop("Nonfinite model predictions.")

  list(predER = predER, predLG = predLG, cumHR = cumHR, bgHR = bgHR,
       n = length(rows), events = sum(dat$.event), periods = nrow(dat),
       IQR_wild = IQR_wild, IQR_bg = IQR_bg, AIC = AIC(mod),
       basis = list(wildfire = attributes(cb_wild)[c("argvar", "arglag", "lag")],
                    background = attributes(cb_bg)[c("argvar", "arglag", "lag")]))
}

############################################################
## 8) RUN OVERALL ASD + ADHD
############################################################

res_overall <- list()
overall_rows <- list()
for (i in seq_len(nrow(outcomes))) {
  nm <- outcomes$name[i]
  cat("Overall model:", nm, "\n")
  res_overall[[nm]] <- fit_overall(
    outcomes$fu[i], outcomes$ev[i], data_wild2_excl,
    revexp_s3, revexp_ds3, IQR_wild, IQR_bg
  )
  overall_rows[[nm]] <- result_row(res_overall[[nm]], nm, "Overall")
}
overall_table <- save_table(overall_rows, "overall_table.csv")
print(overall_table)

############################################################
## 9) SEX-STRATIFIED SETUP
############################################################

sex_levels <- c(0, 1)
sex_labels <- c("Female", "Male")

############################################################
## 10) HELPER: SEX-STRATIFIED MODEL
############################################################
fit_sex_strat <- function(fu_var, ev_var, sex_value) {
  keep <- which(data_wild2_excl$Infant_Gender_Male == sex_value)
  fit_overall(fu_var, ev_var, data_wild2_excl[keep, , drop = FALSE],
              revexp_s3, revexp_ds3, IQR_wild, IQR_bg)
}

############################################################
## 11) RUN SEX-STRATIFIED ASD + ADHD
############################################################

results_sex <- list()
sex_rows <- list()
for (i in seq_len(nrow(outcomes))) {
  for (j in seq_along(sex_levels)) {
    key <- paste(outcomes$name[i], sex_labels[j], sep = "_")
    results_sex[[key]] <- fit_sex_strat(outcomes$fu[i], outcomes$ev[i], sex_levels[j])
    sex_rows[[key]] <- result_row(results_sex[[key]], outcomes$name[i], sex_labels[j])
  }
}
sex_table <- save_table(sex_rows, "sex_table.csv")

############################################################
## 12) HELPER: RACE- AND REGION-STRATIFIED MODELS
############################################################

fit_group_strat <- function(fu_var, ev_var, group_var, group_value) {
  keep <- which(data_wild2_excl[[group_var]] == group_value)
  covariates <- if (group_var == "race_factor") {
    setdiff(full_covariates, "race_factor")
  } else full_covariates
  fit_overall(fu_var, ev_var, data_wild2_excl[keep, , drop = FALSE],
              revexp_s3, revexp_ds3, IQR_wild, IQR_bg, covariates = covariates)
}

############################################################
## 13) RACE-STRATIFIED ASD + ADHD
############################################################

results_race <- list()
race_rows <- list()
for (i in seq_len(nrow(outcomes))) {
  for (rc in race_levels) {
    key <- paste(outcomes$name[i], make.names(rc), sep = "_")
    results_race[[key]] <- fit_group_strat(outcomes$fu[i], outcomes$ev[i], "race_factor", rc)
    race_rows[[key]] <- result_row(results_race[[key]], outcomes$name[i], rc)
  }
}
race_table <- save_table(race_rows, "race_table.csv")

############################################################
## 14) REGION-STRATIFIED ASD + ADHD
############################################################

results_region <- list()
region_rows <- list()
for (i in seq_len(nrow(outcomes))) {
  for (rg in region_levels) {
    key <- paste(outcomes$name[i], rg, sep = "_")
    results_region[[key]] <- fit_group_strat(outcomes$fu[i], outcomes$ev[i], "region4", rg)
    region_rows[[key]] <- result_row(results_region[[key]], outcomes$name[i], rg)
  }
}
region_table <- save_table(region_rows, "region_table.csv")

############################################################
## 15) GRAPHS: OVERALL
############################################################


plot_results <- function(results, label, columns = 2L) {
  panels <- length(results)
  rows <- ceiling(panels / columns)
  for (type in c("ER", "LAG")) {
    png(file.path(output_dir, paste0(type, "_", label, ".png")),
        width = 650 * columns, height = 500 * rows, res = 150)
    tryCatch({
      par(mfrow = c(rows, columns), mar = c(4, 4, 3, 1))
      for (nm in names(results)) {
        fit <- results[[nm]]
        if (type == "ER") {
          plot(fit$predER, "overall", xlab = "Wildfire PM2.5 (ug/m3)",
               ylab = "HR", main = gsub("_", " ", nm), ci = "area")
          abline(v = fit$IQR_wild, lty = 3)
        } else {
          plot(fit$predLG, "slices", var = fit$IQR_wild,
               xlab = "Modeled prenatal lag (weeks)", ylab = "HR",
               main = gsub("_", " ", nm), ci = "area")
        }
        abline(h = 1, lty = 2)
      }
    }, finally = dev.off())
  }
}
plot_results(res_overall, "ASD_ADHD_overall")

############################################################
## 16) GRAPHS: SEX-STRATIFIED
############################################################

plot_results(results_sex, "ASD_ADHD_sex_stratified")

############################################################
## 17) GRAPHS: RACE-STRATIFIED
############################################################

plot_results(results_race, "ASD_ADHD_race_stratified", columns = 4L)

############################################################
## 18) GRAPHS: REGION-STRATIFIED
############################################################

plot_results(results_region, "ASD_ADHD_region_stratified", columns = 4L)

############################################################
## 19) SAVE OVERALL AND SUBGROUP SUMMARY
############################################################

model_summary <- rbind(overall_table, sex_table, race_table, region_table)
write.csv(model_summary, file.path(output_dir, "overall_and_subgroup_results.csv"), row.names = FALSE)

############################################################
## 20) DESCRIBE PARTICIPANTS, FOLLOW-UP AND AGE-8 INCIDENCE
############################################################
descriptives <- function(d,outdir) {
  # Pregnancy-average exposures for descriptive summaries.
  threshold <- median(d$pregnancy_wildfire_mean,na.rm=TRUE)
  groups <- list(Total=seq_len(nrow(d)),Low=which(d$pregnancy_wildfire_mean<=threshold),
                High=which(d$pregnancy_wildfire_mean>threshold))
  continuous <- c('Age','medhouseholdincome','education','nearest_hospital_km',
                  'pregnancy_wildfire_mean','pregnancy_background_mean','poverty')
  categorical <- c('Infant_Gender_Male','race_factor','region4','UA')
  rows <- list()
  for(g in names(groups)) {
    a <- d[groups[[g]],,drop=FALSE]
    rows[[length(rows)+1L]] <- data.frame(group=g,variable='N',level='',n=nrow(a),missing=0,
                                        mean=NA_real_,sd=NA_real_,percent=NA_real_)
    for(v in continuous) {
      x<-a[[v]]
      rows[[length(rows)+1L]] <- data.frame(group=g,variable=v,level='',n=sum(!is.na(x)),
        missing=sum(is.na(x)),mean=mean(x,na.rm=TRUE),sd=sd(x,na.rm=TRUE),percent=NA_real_)
    }
    for(v in categorical) for(l in unique(as.character(d[[v]][!is.na(d[[v]])]))) {
      x<-a[[v]];n<-sum(x==l,na.rm=TRUE)
      rows[[length(rows)+1L]] <- data.frame(group=g,variable=v,level=l,n=n,
        missing=sum(is.na(x)),mean=NA_real_,sd=NA_real_,percent=100*n/sum(!is.na(x)))
    }
  }
  write.csv(do.call(rbind,rows),file.path(outdir,'table1_long.csv'),row.names=FALSE)
  write.csv(data.frame(pregnancy_mean_median=threshold),file.path(outdir,'table1_threshold.csv'),row.names=FALSE)
  follows<-list();km<-list()
  for(i in seq_len(nrow(outcomes))) {
    o <- outcomes$name[i]
    tt<-d[[outcomes$fu[i]]];ev<-d[[outcomes$ev[i]]]
    for(status in c(0,1)) {
      x<-tt[which(ev==status)]/365.25
      follows[[length(follows)+1L]]<-data.frame(outcome=o,status=status,metric='analysis_time',
        n=sum(!is.na(x)),mean_years=mean(x,na.rm=TRUE),sd_years=sd(x,na.rm=TRUE),median_years=median(x,na.rm=TRUE))
      x<-d$Infant_length_enrollment[which(ev==status)]/365.25
      follows[[length(follows)+1L]]<-data.frame(outcome=o,status=status,metric='total_enrollment',
        n=sum(!is.na(x)),mean_years=mean(x,na.rm=TRUE),sd_years=sd(x,na.rm=TRUE),median_years=median(x,na.rm=TRUE))
    }
    for(sample in 'full_cohort') {
      k<-!is.na(tt)&!is.na(ev)
      z<-data.frame(time=tt[k],event=ev[k])
      sf<-survfit(Surv(time,event)~1,data=z)
      # No extrapolation to age 8 if not supported by observed follow-up.
      ss<-summary(sf,times=8*365.25,extend=FALSE)
      if(length(ss$surv)) km[[length(km)+1L]]<-data.frame(outcome=o,sample=sample,
        age=8,incidence=1-ss$surv,lower=1-ss$upper,upper=1-ss$lower,n_risk=ss$n.risk)
    }
  }
  write.csv(do.call(rbind,follows),file.path(outdir,'followup_metrics.csv'),row.names=FALSE)
  if(length(km)) write.csv(do.call(rbind,km),file.path(outdir,'incidence_age8.csv'),row.names=FALSE)
}

descriptives(data_wild2_excl, output_dir)

############################################################
## 21) AIC COMPARISON OF SPLINE DEGREES OF FREEDOM
############################################################

for (i in seq_len(nrow(outcomes))) {
  dat <- analytic_sample(data_wild2_excl, revexp_s3, revexp_ds3,
                         outcomes$fu[i], outcomes$ev[i])
  aic_rows <- list()
  for (exposure_df in 2:6) {
    for (lag_df in 2:6) {
      fit <- tryCatch(
        fit_overall(outcomes$fu[i], outcomes$ev[i], dat, revexp_s3, revexp_ds3,
                    IQR_wild, IQR_bg, exposure_df = exposure_df, lag_df = lag_df,
                    predictions = FALSE),
        error = function(e) e)
      failed <- inherits(fit, "error")
      aic_rows[[length(aic_rows) + 1L]] <- data.frame(
        outcome = outcomes$name[i], N = nrow(dat), exposure_df = exposure_df, lag_df = lag_df,
        AIC = if (failed) NA_real_ else fit$AIC,
        error = if (failed) conditionMessage(fit) else "")
    }
  }
  aic_table <- do.call(rbind, aic_rows)
  if (!any(is.finite(aic_table$AIC))) stop("All AIC models failed for ", outcomes$name[i])
  aic_table$delta_AIC <- aic_table$AIC - min(aic_table$AIC, na.rm = TRUE)
  write.csv(aic_table, file.path(output_dir, paste0("AIC_", outcomes$name[i], ".csv")), row.names = FALSE)
}

############################################################
## 22) ALTERNATIVE EXPOSURE WINDOWS: 28 AND 37 WEEKS
############################################################

results_lag <- list()
lag_rows <- list()
for (lag_max in c(28L, 37L)) {
  exp_window <- build_exposures(data_wild2_excl, lag_max)
  iqr_wild <- IQR(as.vector(exp_window$wild), na.rm = TRUE)
  iqr_bg <- IQR(as.vector(exp_window$bg), na.rm = TRUE)
  if (!is.finite(iqr_wild) || iqr_wild <= 0 || !is.finite(iqr_bg) || iqr_bg <= 0) {
    stop("Invalid exposure IQR for the ", lag_max, "-week window.")
  }
  for (i in seq_len(nrow(outcomes))) {
    key <- paste(outcomes$name[i], lag_max, "weeks", sep = "_")
    results_lag[[key]] <- fit_overall(outcomes$fu[i], outcomes$ev[i], data_wild2_excl,
      exp_window$wild, exp_window$bg, iqr_wild, iqr_bg)
    lag_rows[[key]] <- result_row(results_lag[[key]], outcomes$name[i], paste0(lag_max, "_weeks"))
  }
}
lag_table <- save_table(lag_rows, "alternative_windows.csv")
plot_results(results_lag, "alternative_windows")
rm(exp_window)

############################################################
## 23) DIAGNOSES AFTER AGES 3 AND 5 YEARS
############################################################
results_age <- list()
age_rows <- list()
for (i in seq_len(nrow(outcomes))) {
  t <- data_wild2_excl[[outcomes$fu[i]]]
  event <- data_wild2_excl[[outcomes$ev[i]]]
  for (age in c(3, 5)) {
    keep <- which(event == 0 | (event == 1 & t > age * 365.25))
    key <- paste(outcomes$name[i], "after_age", age, sep = "_")
    results_age[[key]] <- fit_overall(outcomes$fu[i], outcomes$ev[i],
      data_wild2_excl[keep, , drop = FALSE], revexp_s3, revexp_ds3, IQR_wild, IQR_bg)
    age_rows[[key]] <- result_row(results_age[[key]], outcomes$name[i], paste0("After_age_", age))
  }
}
age_table <- save_table(age_rows, "diagnosis_age_sensitivity.csv")
plot_results(results_age, "diagnosis_age_sensitivity")

############################################################
## 24) LINEAR EXPOSURE-RESPONSE MODEL
############################################################

results_linear <- list()
linear_rows <- list()
for (i in seq_len(nrow(outcomes))) {
  nm <- outcomes$name[i]
  results_linear[[nm]] <- fit_overall(outcomes$fu[i], outcomes$ev[i], data_wild2_excl,
    revexp_s3, revexp_ds3, IQR_wild, IQR_bg, linear = TRUE)
  linear_rows[[nm]] <- result_row(results_linear[[nm]], nm, "Linear")
}
linear_table <- save_table(linear_rows, "linear_sensitivity.csv")
plot_results(results_linear, "linear_sensitivity")

############################################################
## 25) MINIMALLY ADJUSTED MODEL
############################################################

minimal_covariates <- c("Age", "race_factor", "season_conception",
                        "medhouseholdincome", "UA", "temp")
results_minimal <- list()
minimal_rows <- list()
for (i in seq_len(nrow(outcomes))) {
  nm <- outcomes$name[i]
  results_minimal[[nm]] <- fit_overall(outcomes$fu[i], outcomes$ev[i], data_wild2_excl,
    revexp_s3, revexp_ds3, IQR_wild, IQR_bg, covariates = minimal_covariates)
  minimal_rows[[nm]] <- result_row(results_minimal[[nm]], nm, "Minimal")
}
minimal_table <- save_table(minimal_rows, "minimal_model.csv")
plot_results(results_minimal, "minimal_model")

############################################################
## 26) IPCW SENSITIVITY ANALYSIS FOR OVERALL CENSORING
############################################################

estimate_ipcw <- function(data, time, event, denominator_covariates) {
  # 1. Prepare one record per child and the birth-year/county strata.
  needed <- unique(c('.row_index', time, event, 'birth_yr', 'FIPS', denominator_covariates))
  stopifnot(all(needed %in% names(data)), !anyNA(data[needed]),
            !anyDuplicated(data$.row_index), all(is.finite(data[[time]])),
            all(data[[time]] > 0), all(data[[event]] %in% c(0, 1)))
  d <- data
  d$.time <- d[[time]]
  d$.censored <- 1L - d[[event]]
  d$.stratum <- interaction(d$birth_yr, d$FIPS, drop = TRUE, lex.order = TRUE)

  # 2. Fit numerator and denominator models for the censoring hazard.
  fit_censor <- function(predictors) {
    f <- reformulate(c('strata(.stratum)', predictors),
                     response = 'Surv(.time, .censored)')
    withCallingHandlers(
      coxph(f, data = d, ties = 'efron', model = TRUE, na.action = na.fail,
            singular.ok = FALSE),
      warning = function(w) stop('Censoring model: ', conditionMessage(w)))
  }
  numerator <- denominator <- NULL
  if (any(d$.censored == 1L)) {
    numerator <- fit_censor(character(0))
    denominator <- fit_censor(denominator_covariates)
  }

  # 3. Split follow-up at outcome event times within each outcome stratum.
  subjects <- data.frame(row_index = d$.row_index, position = seq_len(nrow(d)),
    tstop = d$.time, event = d[[event]], stratum = d$.stratum)
  pieces <- split(subjects, subjects$stratum, drop = TRUE)
  total_rows <- 0
  for (k in seq_along(pieces)) {
    z <- pieces[[k]]
    cuts <- sort(unique(z$tstop[z$event == 1L]))
    total_rows <- total_rows + sum(findInterval(z$tstop, cuts) + !z$tstop %in% cuts)
    if (total_rows > 20000000) stop('Risk-set expansion exceeds 20 million rows; review available memory.')
    z <- if (length(cuts)) survSplit(Surv(tstop, event) ~ ., data = z, cut = cuts,
                                    start = 'tstart') else transform(z, tstart = 0)
    z$used <- z$tstop %in% cuts
    pieces[[k]] <- z
  }
  intervals <- do.call(rbind, pieces)
  stopifnot(sum(intervals$event) == sum(d[[event]]), any(intervals$used))

  # 4. Calculate log G(t-) using each child's censoring stratum.
  # Centered baseline hazards and fitted linear predictors use the same reference.
  log_G <- function(fit) {
    if (is.null(fit)) return(rep(0, nrow(intervals)))
    h <- basehaz(fit, centered = TRUE)
    labels <- as.character(fit$model[['strata(.stratum)']])
    if (!'strata' %in% names(h)) h$strata <- unique(labels)
    curves <- split(h, as.character(h$strata), drop = TRUE)
    rows <- split(seq_len(nrow(intervals)), labels[intervals$position])
    lp <- if (length(coef(fit))) fit$linear.predictors else rep(0, nrow(d))
    stopifnot(all(is.finite(lp)), all(names(rows) %in% names(curves)))
    answer <- numeric(nrow(intervals))
    for (s in names(rows)) {
      j <- rows[[s]]
      curve <- curves[[s]]
      before <- findInterval(intervals$tstop[j], curve$time, left.open = TRUE)
      H0 <- c(0, curve$hazard)[before + 1L]
      answer[j] <- -exp(log(H0) + lp[intervals$position[j]])
    }
    if (any(!is.finite(answer))) stop('Invalid censoring survival; inspect model fit.')
    answer
  }
  log_den <- log_G(denominator)
  intervals$weight <- exp(log_G(numerator) - log_den)
  if (any(!is.finite(intervals$weight) | intervals$weight <= 0)) {
    stop('Invalid weights; inspect censoring models and positivity.')
  }

  # 5. Summarize weights at outcome risk sets. ESS is calculated within each set.
  used <- intervals$used
  w <- intervals$weight[used]
  riskset <- paste(intervals$stratum[used], intervals$tstop[used], sep = ':')
  scaled <- w / ave(w, riskset, FUN = max)
  sums <- rowsum(cbind(n = 1, w = scaled, w2 = scaled^2), riskset, reorder = FALSE)
  ess_fraction <- sums[, 'w']^2 / (sums[, 'n'] * sums[, 'w2'])
  q <- quantile(w, c(0, .5, .95, .99, 1), names = FALSE)
  diagnostics <- data.frame(N = nrow(d), events = sum(d[[event]]),
    censored = sum(d$.censored), risk_set_memberships = length(w), mean = mean(w),
    CV = if (length(w) > 1L) sd(w / mean(w)) else 0,
    min = q[1], median = q[2], p95 = q[3], p99 = q[4], max = q[5],
    min_G = exp(min(log_den[used])), proportion_G_below_001 = mean(log_den[used] < log(.01)),
    min_riskset_ESS_fraction = min(ess_fraction),
    median_riskset_ESS_fraction = median(ess_fraction))
  if (any(log_den[used] < log(.01))) warning('Some G(t-) values are below 0.01; review positivity.')
  list(intervals = intervals[c('row_index', 'tstart', 'tstop', 'event', 'weight')],
       diagnostics = diagnostics)
}

# 6. Fit the unweighted and weighted DLNMs in the same complete-case sample.
results_compare <- comparison_rows <- list()
for (i in seq_len(nrow(outcomes))) {
  nm <- outcomes$name[i]
  censor_covariates <- c(full_covariates, 'Infant_Gender_Male', 'poverty')
  dat <- analytic_sample(data_wild2_excl, revexp_s3, revexp_ds3,
    outcomes$fu[i], outcomes$ev[i], covariates = censor_covariates)
  cw <- crossbasis(revexp_s3[dat$.row_index, , drop = FALSE], lag = c(1, 33),
    argvar = list(fun = 'ns', df = 3), arglag = list(fun = 'ns', df = 2))
  cb <- crossbasis(revexp_ds3[dat$.row_index, , drop = FALSE], lag = c(1, 33),
    argvar = list(fun = 'lin'), arglag = list(fun = 'ns', df = 2))
  wild_terms <- paste0('.ipcw_w', seq_len(ncol(cw)))
  bg_terms <- paste0('.ipcw_b', seq_len(ncol(cb)))
  censor_data <- cbind(dat, setNames(as.data.frame(cw), wild_terms),
                      setNames(as.data.frame(cb), bg_terms))
  ipcw <- estimate_ipcw(censor_data, outcomes$fu[i], outcomes$ev[i],
                       c(censor_covariates, wild_terms, bg_terms))
  write.csv(ipcw$diagnostics, file.path(output_dir, paste0('IPCW_weights_', nm, '.csv')),
            row.names = FALSE)
  results_compare[[paste0(nm, '_Unweighted')]] <- fit_overall(
    outcomes$fu[i], outcomes$ev[i], dat, revexp_s3, revexp_ds3, IQR_wild, IQR_bg)
  results_compare[[paste0(nm, '_IPCW')]] <- fit_overall(
    outcomes$fu[i], outcomes$ev[i], dat, revexp_s3, revexp_ds3, IQR_wild, IQR_bg,
    intervals = ipcw$intervals)
  for (type in c('Unweighted', 'IPCW')) {
    key <- paste(nm, type, sep = '_')
    comparison_rows[[key]] <- result_row(results_compare[[key]], nm, paste0('IPCW_', type))
  }
  rm(censor_data, dat, cw, cb, ipcw)
}

# 7. Save the HR comparison and the DLNM exposure-response and lag-response plots.
comparison_table <- save_table(comparison_rows, 'main_vs_ipcw_comparison.csv')
plot_results(results_compare, 'unweighted_vs_IPCW')

############################################################
## 27) MATERNAL CLUSTERING FOR REPEATED PREGNANCIES
############################################################

results_maternal <- list()
maternal_rows <- list()
for (i in seq_len(nrow(outcomes))) {
  dat <- analytic_sample(data_wild2_excl, revexp_s3, revexp_ds3,
    outcomes$fu[i], outcomes$ev[i], cluster_var = "mother_cluster")
  dat <- analytic_sample(dat, revexp_s3, revexp_ds3, outcomes$fu[i], outcomes$ev[i])
  for (cluster in c("ZIP_Index", "mother_cluster")) {
    key <- paste(outcomes$name[i], cluster, sep = "_")
    results_maternal[[key]] <- fit_overall(outcomes$fu[i], outcomes$ev[i], dat,
      revexp_s3, revexp_ds3, IQR_wild, IQR_bg, cluster_var = cluster)
    maternal_rows[[key]] <- result_row(results_maternal[[key]], outcomes$name[i], cluster)
  }
}
maternal_table <- save_table(maternal_rows, "maternal_clustering.csv")
plot_results(results_maternal, "maternal_clustering")

############################################################
## 28) ICD-9 AND ICD-10 PERIODS
############################################################
## Split follow-up at October 1, 2015 on each child's age time scale.
## Children observed across the transition contribute person-time to both eras.
## Check that the split preserves total follow-up and outcome events.

as_verified_date <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, "POSIXt")) return(as.Date(x, tz = "UTC"))
  if (is.numeric(x)) stop("Convert numeric dates using their verified origin before analysis.")
  out <- as.Date(x)
  if (any(!is.na(x) & is.na(out))) stop("Unparseable calendar date.")
  out
}
icd_intervals <- function(dat, fu_var, ev_var, diagnosis_var) {
  birth <- as_verified_date(dat$Deliv_Date2)
  diagnosis <- as_verified_date(dat[[diagnosis_var]])
  time <- dat[[fu_var]]; event <- dat[[ev_var]]
  if (anyNA(birth) || any(event == 1 & is.na(diagnosis))) stop("Missing required ICD dates.")
  cases <- which(event == 1)
  if (any(abs(as.numeric(diagnosis[cases] - birth[cases]) - time[cases]) > 0.01)) {
    stop("Diagnosis dates and follow-up times disagree.")
  }
  transition <- as.Date("2015-10-01")
  if (any(event == 1 & diagnosis == transition, na.rm = TRUE)) {
    stop("An event occurs on the ICD boundary; reconcile the within-day timing before splitting.")
  }
  boundary <- as.numeric(transition - birth)
  pre <- data.frame(row_index = dat$.row_index, tstart = 0, tstop = pmin(time, boundary),
    event = as.integer(event == 1 & diagnosis < transition), weight = 1, era = "ICD9")
  post <- data.frame(row_index = dat$.row_index, tstart = pmax(0, boundary), tstop = time,
    event = as.integer(event == 1 & diagnosis >= transition), weight = 1, era = "ICD10")
  intervals <- rbind(pre, post)
  intervals$event[is.na(intervals$event)] <- 0L
  intervals <- intervals[intervals$tstop > intervals$tstart, , drop = FALSE]
  duration <- rowsum(intervals$tstop - intervals$tstart, intervals$row_index, reorder = FALSE)
  events <- rowsum(intervals$event, intervals$row_index, reorder = FALSE)
  index <- match(as.character(dat$.row_index), rownames(duration))
  if (anyNA(index) || any(abs(duration[index, 1] - time) > 1e-7) ||
      any(events[index, 1] != event)) stop("ICD split does not preserve follow-up/events.")
  intervals
}
results_icd <- list()
icd_rows <- list()
for (i in seq_len(nrow(outcomes))) {
  dat <- analytic_sample(data_wild2_excl, revexp_s3, revexp_ds3, outcomes$fu[i], outcomes$ev[i])
  intervals <- icd_intervals(dat, outcomes$fu[i], outcomes$ev[i], outcomes$diagnosis_date[i])
  for (era in c("ICD9", "ICD10")) {
    key <- paste(outcomes$name[i], era, sep = "_")
    era_intervals <- intervals[intervals$era == era, , drop = FALSE]
    results_icd[[key]] <- fit_overall(outcomes$fu[i], outcomes$ev[i], dat,
      revexp_s3, revexp_ds3, IQR_wild, IQR_bg, intervals = era_intervals)
    icd_rows[[key]] <- result_row(results_icd[[key]], outcomes$name[i], era)
  }
}
icd_table <- save_table(icd_rows, "ICD_transition.csv")
plot_results(results_icd, "ICD_transition")

############################################################
## 29) E-VALUES FOR THE PRIMARY WILDFIRE ESTIMATES
############################################################
## Use the unrounded HR as a risk-ratio approximation for rare outcomes.
## Use the confidence bound nearest 1; a CI containing 1 gives E-value 1.

evalue <- function(hr, lcl, ucl) {
  if (any(!is.finite(c(hr, lcl, ucl))) || min(hr, lcl, ucl) <= 0 || lcl > hr || hr > ucl) {
    stop("Invalid HR or confidence interval.")
  }
  f <- function(x) { r <- max(x, 1/x); r + sqrt(r * (r - 1)) }
  bound <- if (lcl <= 1 && ucl >= 1) 1 else if (hr >= 1) lcl else ucl
  c(estimate = unname(f(hr)), CI = unname(f(bound)))
}
evalue_rows <- list()
for (nm in outcomes$name) {
  hr <- res_overall[[nm]]$cumHR
  value <- evalue(hr["HR"], hr["LCL"], hr["UCL"])
  evalue_rows[[nm]] <- data.frame(outcome = nm, HR = unname(hr["HR"]),
    LCL = unname(hr["LCL"]), UCL = unname(hr["UCL"]),
    E_value = unname(value["estimate"]), E_value_CI = unname(value["CI"]))
}
evalue_table <- save_table(evalue_rows, "E_values.csv")

############################################################
## 30) SAVE SENSITIVITY SUMMARY AND SOFTWARE VERSIONS
############################################################
## Export aggregate tables and spline settings. Participant records,
## identifiers, dates and fitted model objects are not written to disk.

sensitivity_summary <- rbind(lag_table, age_table, linear_table, minimal_table,
                              comparison_table, maternal_table, icd_table)
write.csv(sensitivity_summary, file.path(output_dir, "sensitivity_results.csv"), row.names = FALSE)
all_results <- list(overall = res_overall, sex = results_sex, race = results_race,
  region = results_region, windows = results_lag, age = results_age,
  linear = results_linear, minimal = results_minimal, ipcw = results_compare,
  maternal = results_maternal, icd = results_icd)
basis_settings <- lapply(all_results, function(group) lapply(group, function(fit) fit$basis))
dput(basis_settings, file = file.path(output_dir, "basis_settings.txt"))
capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
cat("Analysis completed. Aggregate tables and plots are in:", output_dir, "\n")

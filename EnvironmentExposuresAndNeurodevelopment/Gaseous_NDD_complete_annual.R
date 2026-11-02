## 1) LOAD PACKAGES AND SET INPUTS

library(dlnm)
library(survival)
library(splines)

CFG <- list(
  input_files = list(),
  input_objects = c(main = "dat", `26` = "dat_26", `30` = "dat_30", `32` = "dat_32"),
  unrestricted_input_file = NULL,
  output_dir = "results_Gaseous_NDD_steps",
  study_end = as.Date(NA_character_),
  birth_years = 2001:2013,
  outcomes = c("ASD", "ADHD", "LDif", "DSLD", "DCD", "IntD", "BDis", "AnyNDD"),
  pollutants = c("no2", "ozone"),
  main_prenatal_weeks = 37L,
  sensitivity_weeks = c(26L, 30L, 32L),
  pre_lag_df = c(ap = 4L, temp = 4L, rh = 4L),
  post_lag_df = c(ap = 4L, temp = 4L, rh = 4L),
  increment_ppb = 10,
  id_column = "patient_birth_ID",
  enrollment_column = "Infant_length_enrollment.x",
  conception_date_column = NULL,
  death_date_column = NULL,
  exit_reason_column = NULL,
  loss_to_followup_values = c("disenrollment", "loss_to_followup"),
  other_exit_reason_values = c("death", "administrative", "end_of_study"),
  censor_mode = "disenrollment",
  categorical_covariates = c("season_conception", "Alcohol_abuse_yn", "Smoking_yn",
    "Substance_abuse_yn", "Poor_nutrition_yn", "prenatal_folate_yn"),
  numeric_covariates = c("Age", "BMI", "ndvi", "medhouseholdincome", "pct_owner_occ",
    "medianhousevalue", "education", "poverty", "pop_density_30sec_NASA", "smoke_rate",
    "mean_bmi", "nearest_hospital_km", "pct_amb_visit"),
  column_map = list(),
  missing_exposure_codes = 999,
  use_prefilled_annual_gases = TRUE,
  post_adjust_early_copollutants = FALSE,
  run_subgroups = TRUE,
  subgroup_levels = list(infant_sex_cat = c("Male", "Female"),
    race_cat = c("White", "Black/African American", "Hisp/Latino"),
    population_density_tertile = c("Low", "Medium", "High")),
  population_density_cutpoints = NULL,
  run_sensitivity = TRUE,
  run_ipcw = TRUE,
  ipcw_periods = c("preconception_prenatal", "postnatal"),
  ipcw_truncate_quantile = 0.99,
  ipcw_max_rows = 20000000,
  run_aic = FALSE,
  aic_lag_df = 2:6,
  save_plots = TRUE
)

event_numbers <- c(ASD = 1L, ADHD = 3L, LDif = 5L, DSLD = 6L,
                  DCD = 7L, IntD = 9L, BDis = 12L, AnyNDD = 22L)
outcome_labels <- c(ASD = "Autism spectrum disorder", ADHD = "ADHD",
  LDif = "Learning difficulty", DSLD = "Developmental speech or language disorder",
  DCD = "Developmental coordination disorder", IntD = "Intellectual disability",
  BDis = "Behavioral disorder", AnyNDD = "Any neurodevelopmental disorder")
input_environment <- environment()

## 2) DEFINE INPUT CHECKS AND READ THE APPROPRIATE COHORT

need_columns <- function(d, columns) {
  missing <- setdiff(columns, names(d))
  if (length(missing)) stop("Missing columns: ", paste(missing, collapse = ", "))
}

as_num <- function(x, label) {
  if (is.factor(x)) x <- as.character(x)
  y <- suppressWarnings(as.numeric(x))
  if (any(!is.na(x) & !is.finite(y))) stop("Invalid numeric values in ", label)
  y
}

as_date <- function(x, label) {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, "POSIXt")) return(as.Date(x))
  if (is.numeric(x)) stop("Convert ", label, " to Date using its verified origin.")
  y <- as.Date(as.character(x))
  if (any(!is.na(x) & is.na(y))) stop("Invalid dates in ", label)
  y
}

read_cohort <- function(key, cfg) {
  path <- cfg$input_files[[key]]
  if (is.null(path)) path <- cfg$unrestricted_input_file
  if (!is.null(path)) {
    if (!file.exists(path)) stop("Input file does not exist: ", path)
    d <- readRDS(path)
  } else {
    object <- unname(cfg$input_objects[key])
    if (length(object) != 1L || is.na(object) ||
        !exists(object, envir = input_environment, inherits = TRUE))
      stop("Provide CFG$input_files[['", key, "']] or its prepared input object. ",
           "For 26/30/32-week analyses use the corresponding eligible cohort, ",
           "or CFG$unrestricted_input_file before prenatal-completeness selection.")
    d <- get(object, envir = input_environment, inherits = TRUE)
  }
  if (!is.data.frame(d) || !nrow(d) || anyDuplicated(names(d)))
    stop("Input must be a nonempty prepared cohort with unique column names.")
  need_columns(d, c(cfg$id_column, "Deliv_Date", "birth_yr", "FIPS", "ZIP_Index"))
  if (anyNA(d[[cfg$id_column]]) || anyDuplicated(d[[cfg$id_column]]))
    stop("Expected one row per child; resolve missing or duplicated IDs.")
  d$birth_yr <- as_num(d$birth_yr, "birth_yr")
  d <- d[!is.na(d$birth_yr) & d$birth_yr %in% cfg$birth_years, , drop = FALSE]
  if (!nrow(d)) stop("No births in CFG$birth_years.")
  d
}

## 3) PREPARE CONCEPTION SEASON, COVARIATES AND STRATA

prepare_cohort <- function(d, cfg, density_cutpoints = NULL) {
  d$Deliv_Date <- as_date(d$Deliv_Date, "Deliv_Date")
  if (anyNA(d$Deliv_Date) || any(as.integer(format(d$Deliv_Date, "%Y")) != d$birth_yr))
    stop("Birth dates are missing or disagree with birth_yr.")
  if (!"season_conception" %in% names(d)) {
    if ("season_conception_cat" %in% names(d)) d$season_conception <- d$season_conception_cat
    else {
      if (is.null(cfg$conception_date_column))
        stop("Provide season_conception, or set CFG$conception_date_column to the verified prenatal time-origin date column.")
      need_columns(d, cfg$conception_date_column)
      origin <- as_date(d[[cfg$conception_date_column]], cfg$conception_date_column)
      if (any(origin > d$Deliv_Date, na.rm = TRUE)) stop("Prenatal origin is after birth.")
      month <- as.integer(format(origin, "%m"))
      d$season_conception <- c("Winter", "Spring", "Summer", "Autumn")[(month %% 12L) %/% 3L + 1L]
    }
  }
  if (!"race_cat" %in% names(d)) {
    columns <- c("Race_WHITE", "Race_BLACK_OR_AFRICAN_AMERICAN", "Race_HISPANIC_OR_LATINO",
                 "Race_ASIAN_OR_OTHER_PACIFIC", "Race_UNKNOWN_OTHER")
    need_columns(d, columns)
    x <- vapply(columns, function(v) as_num(d[[v]], v), numeric(nrow(d)))
    if (any(!is.na(x) & !x %in% 0:1)) stop("Race indicators must be 0/1/NA.")
    ok <- complete.cases(x)
    if (any(rowSums(x[ok, , drop = FALSE]) != 1L)) stop("Resolve overlapping or unclassified race indicators.")
    d$race_cat <- NA_character_
    d$race_cat[ok] <- c("White", "Black/African American", "Hisp/Latino", "Asian/Pacific",
                       "Unknown/Other")[max.col(x[ok, , drop = FALSE], ties.method = "first")]
  }
  if (!"infant_sex_cat" %in% names(d)) {
    columns <- c("Infant_Gender_Male", "Infant_Gender_Female", "Infant_Gender_Unknown")
    need_columns(d, columns)
    x <- vapply(columns, function(v) as_num(d[[v]], v), numeric(nrow(d)))
    if (any(!is.na(x) & !x %in% 0:1)) stop("Sex indicators must be 0/1/NA.")
    ok <- complete.cases(x)
    if (any(rowSums(x[ok, , drop = FALSE]) != 1L)) stop("Resolve overlapping or unclassified sex indicators.")
    d$infant_sex_cat <- NA_character_
    d$infant_sex_cat[ok] <- c("Male", "Female", "Unknown")[max.col(x[ok, , drop = FALSE], ties.method = "first")]
  }
  aliases <- c(Alcohol_abuse_yn = "Alcohol_abuse", Smoking_yn = "Smoking",
    Substance_abuse_yn = "Substance_abuse", Poor_nutrition_yn = "Poor_nutrition",
    prenatal_folate_yn = "prenatal_folate")
  for (v in names(aliases)) if (!v %in% names(d) && aliases[[v]] %in% names(d))
    d[[v]] <- d[[aliases[[v]]]]
  need_columns(d, c(cfg$categorical_covariates, cfg$numeric_covariates))
  for (v in cfg$numeric_covariates) d[[v]] <- as_num(d[[v]], v)
  for (v in c(cfg$categorical_covariates, "race_cat", "infant_sex_cat", "FIPS", "ZIP_Index"))
    d[[v]] <- factor(d[[v]])
  if ("region" %in% names(d)) d$region <- factor(d$region)
  d$.id <- seq_len(nrow(d))
  d$.zip <- d$ZIP_Index
  d$.stratum <- interaction(d$FIPS, d$birth_yr, d$race_cat, d$infant_sex_cat,
                            drop = TRUE, lex.order = TRUE)
  if (!is.null(density_cutpoints)) {
    if (length(density_cutpoints) != 2L || any(!is.finite(density_cutpoints)) ||
        density_cutpoints[1] >= density_cutpoints[2]) stop("Population-density tertile cutpoints must be two distinct numbers.")
    d$population_density_tertile <- cut(d$pop_density_30sec_NASA,
      breaks = c(-Inf, density_cutpoints, Inf), labels = c("Low", "Medium", "High"), right = TRUE)
  }
  d
}

## 4) BUILD CHRONOLOGICAL EXPOSURE MATRICES AND RECORD THEIR MAPPING

exposure_columns <- function(d, variable, period, npren, cfg) {
  key <- paste(variable, period, sep = "_")
  override <- cfg$column_map[[key]]
  index <- switch(period, precon = 12:1, prenatal = seq_len(npren), postnatal = 0:3)
  if (!is.null(override)) {
    if (period == "prenatal" && length(override) >= npren) override <- override[seq_len(npren)]
    if (!is.character(override) || length(override) != length(index)) stop("Invalid chronological column map for ", key)
    need_columns(d, override)
    return(override)
  }
  prefix <- if (variable == "temp") c("temp", "tmpC") else variable
  columns <- vapply(index, function(i) {
    candidates <- if (period == "postnatal") {
      standard <- c(paste0(prefix, "_postnatal_yr.", i), paste0(prefix, "Year.", i),
                    paste0(prefix, ".Year", i))
      filled <- if (variable %in% c("no2", "ozone") && cfg$use_prefilled_annual_gases)
        paste0(variable, ".Year", i, "_fill") else character()
      c(filled, standard)
    } else paste0(prefix, "_", period, "_wk.", i)
    found <- intersect(candidates, names(d))
    if (!length(found)) stop("Missing exposure column: ", paste(candidates, collapse = " or "))
    found[1L]
  }, character(1))
  if (anyDuplicated(columns)) stop("Duplicated exposure columns for ", key)
  columns
}

exposure_matrices <- function(d, npren, cfg, include_post = TRUE) {
  matrices <- list(); audit <- list()
  periods <- c("precon", "prenatal", if (include_post) "postnatal")
  for (v in c("no2", "ozone", "pm25", "temp", "rh")) for (period in periods) {
    key <- paste(v, period, sep = "_")
    cols <- exposure_columns(d, v, period, npren, cfg)
    x <- vapply(cols, function(z) as_num(d[[z]], z), numeric(nrow(d)))
    missing_code <- x %in% cfg$missing_exposure_codes
    dim(missing_code) <- dim(x)
    x[missing_code] <- NA_real_
    if (v %in% c("no2", "ozone", "pm25") && any(x < 0, na.rm = TRUE))
      stop("Negative pollutant values in ", key)
    if (v == "rh" && any(x < 0 | x > 100, na.rm = TRUE)) stop("RH must be in percent, 0-100.")
    matrices[[key]] <- x
    biological <- switch(period, precon = -(12:1), prenatal = seq_len(npren), postnatal = 0:3)
    coordinate <- if (period == "prenatal") biological - 1L else biological
    audit[[key]] <- data.frame(variable = v, period = period, position = seq_along(cols),
      column = cols, biological_week_or_year = biological, model_coordinate = coordinate,
      missing_n = colSums(is.na(x)), sentinel_n = colSums(missing_code))
  }
  list(matrices = matrices, audit = do.call(rbind, audit))
}

complete_early <- function(ex) {
  Reduce(`&`, lapply(ex$matrices[!grepl("_postnatal$", names(ex$matrices))], complete.cases))
}

## 5) DEFINE OUTCOME-SPECIFIC FOLLOW-UP AND CENSORING

prepare_followup <- function(d, outcome, cfg) {
  number <- event_numbers[[outcome]]
  flag <- paste0("Event_", number, "_MAIN_correct")
  duration <- paste0("Event_", number, "_FU_correct")
  enrollment_column <- cfg$enrollment_column
  if (!enrollment_column %in% names(d) && enrollment_column == "Infant_length_enrollment.x" &&
      "Infant_length_enrollment" %in% names(d)) enrollment_column <- "Infant_length_enrollment"
  need_columns(d, c(flag, duration, enrollment_column))
  event <- as_num(d[[flag]], flag)
  time <- as_num(d[[duration]], duration)
  enrollment <- as_num(d[[enrollment_column]], enrollment_column)
  if (anyNA(event) || any(!event %in% 0:1) || anyNA(enrollment) || any(enrollment < 0) ||
      any(event == 1 & (is.na(time) | time <= 0))) stop("Invalid corrected outcome or enrollment duration for ", outcome)
  admin <- rep(Inf, nrow(d))
  known_admin <- length(cfg$study_end) == 1L && !is.na(cfg$study_end)
  if (known_admin) admin <- as.numeric(cfg$study_end - d$Deliv_Date)
  if (any(admin <= 0)) stop("Study end must follow every included birth date.")
  death <- rep(Inf, nrow(d))
  if (!is.null(cfg$death_date_column)) {
    need_columns(d, cfg$death_date_column)
    date <- as_date(d[[cfg$death_date_column]], cfg$death_date_column)
    death[!is.na(date)] <- as.numeric(date[!is.na(date)] - d$Deliv_Date[!is.na(date)])
    if (any(death < 0)) stop("Death date precedes birth.")
  }
  end <- pmin(enrollment, admin, death)
  observed <- event == 1 & !is.na(time) & time <= end
  d$.event <- as.integer(observed)
  d$.stop <- ifelse(observed, time, end)
  d$.start <- 0
  d$.censor <- NA_integer_
  definition <- "administrative_end_not_supplied"
  if (known_admin) {
    early <- !observed & end < admin
    if (cfg$censor_mode == "all_early_exits") {
      d$.censor <- as.integer(early)
      definition <- "all_nonadministrative_exits_including_death"
    } else if (cfg$censor_mode == "disenrollment") {
      if (is.null(cfg$exit_reason_column)) definition <- "disenrollment_reasons_not_supplied"
      else {
        need_columns(d, cfg$exit_reason_column)
        reason <- as.character(d[[cfg$exit_reason_column]])
        d$.censor <- as.integer(early & reason %in% cfg$loss_to_followup_values & enrollment < death)
        recognized <- reason %in% c(cfg$loss_to_followup_values, cfg$other_exit_reason_values)
        d$.censor[early & (!recognized | is.na(reason) | !nzchar(trimws(reason)))] <- NA_integer_
        definition <- "specified_disenrollment_reasons"
      }
    } else stop("Unknown CFG$censor_mode.")
  }
  attr(d, "followup_audit") <- data.frame(outcome = outcome, input_n = nrow(d),
    input_events = sum(event), observed_events = sum(d$.event),
    events_outside_followup = sum(event == 1 & !observed), nonpositive_followup = sum(d$.stop <= 0),
    administrative_end = if (known_admin) as.character(cfg$study_end) else NA_character_,
    censor_definition = definition)
  attr(d, "censor_definition") <- definition
  d
}

## 6) PREPARE COMPLETE-CASE MODELS AND CO-POLLUTANT ADJUSTMENT

make_analysis <- function(d, ex, outcome, pollutant, period, npren, cfg) {
  if (!pollutant %in% c("no2", "ozone")) stop("Unknown pollutant.")
  early <- function(v) cbind(ex$matrices[[paste0(v, "_precon")]], ex$matrices[[paste0(v, "_prenatal")]])
  post <- function(v) ex$matrices[[paste0(v, "_postnatal")]]
  copollutants <- c("pm25", if (pollutant == "no2") "ozone" else "no2")
  covariates <- c(cfg$categorical_covariates, cfg$numeric_covariates)
  is_pre <- period == "preconception_prenatal"
  if (is_pre) {
    X <- list(ap = early(pollutant), temp = early("temp"), rh = early("rh"))
    lag <- c(-12L, npren - 1L)
    df <- cfg$pre_lag_df
    d$.copollutant_pm25 <- rowMeans(early(copollutants[1]))
    d$.copollutant_gas <- rowMeans(early(copollutants[2]))
  } else {
    if (period != "postnatal") stop("Unknown exposure period.")
    X <- list(ap = post(pollutant), temp = post("temp"), rh = post("rh"))
    if (any(vapply(X, ncol, integer(1)) != 4L)) stop("Postnatal models require four annual values.")
    lag <- c(0L, 3L)
    df <- cfg$post_lag_df
    d$.copollutant_pm25 <- rowMeans(post(copollutants[1]))
    d$.copollutant_gas <- rowMeans(post(copollutants[2]))
    d$.main_early_mean <- rowMeans(early(pollutant))
    covariates <- c(covariates, ".main_early_mean")
    if (cfg$post_adjust_early_copollutants) {
      d$.copollutant_pm25_early <- rowMeans(early(copollutants[1]))
      d$.copollutant_gas_early <- rowMeans(early(copollutants[2]))
      covariates <- c(covariates, ".copollutant_pm25_early", ".copollutant_gas_early")
    }
  }
  covariates <- c(covariates, ".copollutant_pm25", ".copollutant_gas")
  required <- c(covariates, ".id", ".zip", ".stratum", ".start", ".stop", ".event")
  need_columns(d, required)
  keep <- complete_early(ex) & complete.cases(d[required]) & d$.stop > d$.start
  for (x in X) keep <- keep & complete.cases(x)
  grouping <- intersect(names(cfg$subgroup_levels), names(d))
  dat <- droplevels(d[keep, unique(c(required, ".censor", grouping)), drop = FALSE])
  X <- lapply(X, function(x) x[keep, , drop = FALSE])
  if (!nrow(dat) || !sum(dat$.event)) stop("No eligible events for ", outcome, " / ", period)
  numeric_fields <- setdiff(names(dat)[vapply(dat, is.numeric, logical(1))], ".censor")
  if (any(!is.finite(as.matrix(dat[numeric_fields])))) stop("Nonfinite covariates.")
  constant <- covariates[vapply(dat[covariates], function(x) length(unique(x)) < 2L, logical(1))]
  reference_column <- if (is_pre) 12L else 1L
  reference <- round(unname(quantile(X$ap[, reference_column], 0.25)))
  list(data = dat, X = X, lag = lag, df = df, covariates = setdiff(covariates, constant),
    constant = constant, outcome = outcome, pollutant = pollutant, period = period, npren = npren,
    reference = reference, increment = cfg$increment_ppb, input_n = nrow(d),
    excluded_n = sum(!keep), censor_definition = attr(d, "censor_definition"),
    copollutants = copollutants, post_adjust_early_copollutants = cfg$post_adjust_early_copollutants)
}

subset_analysis <- function(a, variable, value) {
  need_columns(a$data, variable)
  keep <- !is.na(a$data[[variable]]) & as.character(a$data[[variable]]) == value
  b <- a
  b$data <- droplevels(a$data[keep, , drop = FALSE])
  b$X <- lapply(a$X, function(x) x[keep, , drop = FALSE])
  if (!nrow(b$data) || !sum(b$data$.event)) stop("No eligible events in subgroup ", variable, " = ", value)
  constant <- b$covariates[vapply(b$data[b$covariates], function(x) length(unique(x)) < 2L, logical(1))]
  b$covariates <- setdiff(b$covariates, constant)
  b$constant <- unique(c(b$constant, constant))
  b$input_n <- nrow(a$data); b$excluded_n <- sum(!keep)
  b
}

## 7) DEFINE CROSS-BASES AND FIT STRATIFIED COX MODELS

make_bases <- function(a, cfg) {
  lapply(setNames(c("ap", "temp", "rh"), c("ap", "temp", "rh")), function(v) {
    x <- a$X[[v]]
    if (ncol(x) != diff(a$lag) + 1L || a$df[[v]] > ncol(x)) stop("Exposure columns and lag specification disagree.")
    crossbasis(x, lag = a$lag, argvar = list(fun = "lin"), arglag = list(fun = "ns", df = a$df[[v]]))
  })
}

expand_basis <- function(cb, index) {
  old <- attributes(cb); ans <- cb[index, , drop = FALSE]
  for (v in setdiff(names(old), c("dim", "dimnames"))) attr(ans, v) <- old[[v]]
  ans
}

fit_analysis <- function(a, cfg, intervals = NULL, predictions = TRUE) {
  basis <- make_bases(a, cfg)
  dat <- a$data
  dat$.weight <- 1
  if (!is.null(intervals)) {
    position <- match(intervals$.id, dat$.id)
    if (anyNA(position)) stop("IPCW interval IDs do not match the model cohort.")
    dat <- dat[position, , drop = FALSE]
    for (v in c(".start", ".stop", ".event", ".weight")) dat[[v]] <- intervals[[v]]
    basis <- lapply(basis, expand_basis, index = position)
  }
  cb_ap <- basis$ap; cb_temp <- basis$temp; cb_rh <- basis$rh
  formula <- as.formula(paste("Surv(.start, .stop, .event) ~", paste(c("cb_ap", "cb_temp", "cb_rh",
    a$covariates, "strata(.stratum)"), collapse = " + ")))
  model <- coxph(formula, data = dat, weights = .weight, cluster = .zip, robust = TRUE,
    ties = "efron", na.action = na.fail, singular.ok = FALSE, x = FALSE, y = FALSE, model = FALSE)
  if (any(!is.finite(coef(model)))) stop("Nonfinite model coefficients.")
  if (!predictions) return(list(AIC = AIC(model)))
  coefficient_names <- paste0("cb_ap", colnames(cb_ap))
  if (!all(coefficient_names %in% names(coef(model)))) stop("Cannot identify exposure-basis coefficients.")
  beta <- coef(model)[coefficient_names]
  covariance <- vcov(model)[coefficient_names, coefficient_names, drop = FALSE]
  predict_ap <- function(at, lag = a$lag) crosspred(cb_ap, coef = beta, vcov = covariance,
    model.link = "log", cen = a$reference, at = at, lag = lag, bylag = 1)
  prediction <- predict_ap(a$reference + a$increment)
  grid <- seq(min(a$X$ap), max(a$X$ap), length.out = 101L)
  curve <- predict_ap(sort(unique(c(grid, a$reference, a$reference + a$increment))))
  windows <- if (a$period == "postnatal") list(year_0 = c(0, 0), year_1 = c(1, 1),
    year_2 = c(2, 2), year_3 = c(3, 3), postnatal_all = c(0, 3)) else {
    w <- list(preconception = c(-12, -1), trimester_1 = c(0, 11),
      trimester_2 = c(12, min(26, a$npren - 1)), prenatal = c(0, a$npren - 1),
      prenatal_weeks_10_16 = c(9, 15), preconception_prenatal_all = a$lag)
    if (a$npren >= 28L) w$trimester_3 <- c(27, a$npren - 1L)
    w
  }
  window_table <- do.call(rbind, lapply(names(windows), function(name) {
    p <- predict_ap(a$reference + a$increment, windows[[name]])
    data.frame(window = name, coordinate_start = windows[[name]][1], coordinate_end = windows[[name]][2],
      HR = unname(p$allRRfit), LCL = unname(p$allRRlow), UCL = unname(p$allRRhigh))
  }))
  coords <- seq.int(a$lag[1], a$lag[2])
  lag_table <- data.frame(coordinate = coords,
    period = if (a$period == "postnatal") "postnatal" else ifelse(coords < 0, "preconception", "prenatal"),
    week_or_year = if (a$period == "postnatal") coords else ifelse(coords < 0, coords, coords + 1L),
    HR = as.numeric(prediction$matRRfit), LCL = as.numeric(prediction$matRRlow),
    UCL = as.numeric(prediction$matRRhigh))
  list(prediction = prediction, curve = curve, window_HRs = window_table, lag_HRs = lag_table,
    beta = beta, vcov = covariance, AIC = AIC(model), N = nrow(a$data), events = sum(a$data$.event),
    intervals = nrow(dat), formula = paste(deparse(formula), collapse = " "),
    basis = lapply(basis, function(x) attributes(x)[c("argvar", "arglag", "lag")]))
}

## 8) EXTRACT HAZARD RATIOS AND DRAW DLNM PLOTS

result_table <- function(fit, a, label) {
  cbind(data.frame(analysis = label, outcome = a$outcome, outcome_label = unname(outcome_labels[a$outcome]),
    pollutant = a$pollutant, exposure_period = a$period, prenatal_weeks = a$npren,
    N = fit$N, events = fit$events, increment_ppb = a$increment, reference_ppb = a$reference,
    lag_df = a$df[["ap"]], AIC = fit$AIC), fit$window_HRs)
}

save_plots <- function(fit, a, prefix) {
  pdf(paste0(prefix, "_plots.pdf"), width = 10, height = 4.8)
  on.exit(dev.off(), add = TRUE)
  par(mfrow = c(1, 2), mar = c(5, 4, 3, 1))
  plot(fit$curve, "overall", ci = "area", xlab = paste(a$pollutant, "(ppb)"),
    ylab = "Cumulative HR", main = paste(a$outcome, a$pollutant))
  abline(h = 1, lty = 2)
  plot(fit$prediction, "slices", var = a$reference + a$increment, ci = "area", xaxt = "n",
    xlab = if (a$period == "postnatal") "Postnatal year (0 = birth year)" else "Preconception / prenatal week",
    ylab = paste0("HR per ", a$increment, " ppb"), main = paste(a$outcome, a$pollutant))
  if (a$period == "postnatal") axis(1, at = 0:3, labels = 0:3)
  else {
    ticks <- unique(c(-12, -8, -4, 0, seq(8, a$npren - 1, by = 8), a$npren - 1))
    axis(1, at = ticks, labels = ifelse(ticks < 0, ticks, ticks + 1L))
    abline(v = -0.5, lty = 3)
  }
  abline(h = 1, lty = 2)
}

## 9) DEFINE INDEPENDENT DESCRIPTIVE TABLES AND EXPOSURE CORRELATIONS

save_descriptives <- function(d, ex, cfg, prefix) {
  eligible <- complete_early(ex)
  d <- droplevels(d[eligible, , drop = FALSE])
  e <- lapply(ex$matrices, function(x) x[eligible, , drop = FALSE])
  if (!nrow(d)) stop("No complete early exposure histories for descriptive tables.")
  for (v in c("no2", "ozone", "pm25", "temp", "rh")) {
    d[[paste0(v, "_precon_mean")]] <- rowMeans(e[[paste0(v, "_precon")]])
    d[[paste0(v, "_prenatal_mean")]] <- rowMeans(e[[paste0(v, "_prenatal")]])
    if (!is.null(e[[paste0(v, "_postnatal")]])) {
      d[[paste0(v, "_postnatal_mean")]] <- rowMeans(e[[paste0(v, "_postnatal")]])
      for (j in 0:3) d[[paste0(v, "_year_", j)]] <- e[[paste0(v, "_postnatal")]][, j + 1L]
    }
  }
  groups <- list(all = rep(TRUE, nrow(d)))
  for (v in c("no2", "ozone")) {
    x <- d[[paste0(v, "_prenatal_mean")]]
    cutoff <- median(x, na.rm = TRUE)
    groups[[paste0(v, "_below_median")]] <- !is.na(x) & x < cutoff
    groups[[paste0(v, "_at_or_above_median")]] <- !is.na(x) & x >= cutoff
    write.csv(data.frame(variable = v, prenatal_mean_median = cutoff),
      paste0(prefix, "_", v, "_median.csv"), row.names = FALSE)
  }
  numeric_columns <- unique(c(cfg$numeric_covariates, grep("_(mean|year_[0-3])$", names(d), value = TRUE)))
  categorical_columns <- unique(c(cfg$categorical_covariates, "race_cat", "infant_sex_cat",
                                  intersect(c("region", "population_density_tertile"), names(d))))
  rows <- list(); k <- 0L
  for (g in names(groups)) {
    z <- d[groups[[g]], , drop = FALSE]
    for (v in numeric_columns) {
      x <- z[[v]]; good <- x[is.finite(x)]
      q <- if (length(good)) quantile(good, c(.25, .5, .75), names = FALSE) else rep(NA_real_, 3)
      k <- k + 1L
      rows[[k]] <- data.frame(group = g, variable = v, level = "", N = nrow(z),
        observed_n = length(good), missing_n = sum(is.na(x)), n = NA_integer_, percent = NA_real_,
        mean = if (length(good)) mean(good) else NA_real_, sd = if (length(good) > 1) sd(good) else NA_real_,
        p25 = q[1], median = q[2], p75 = q[3])
    }
    for (v in categorical_columns) {
      counts <- table(as.character(z[[v]]), useNA = "no")
      for (level in names(counts)) {
        k <- k + 1L
        rows[[k]] <- data.frame(group = g, variable = v, level = level, N = nrow(z),
          observed_n = sum(counts), missing_n = sum(is.na(z[[v]])), n = unname(counts[level]),
          percent = 100 * unname(counts[level]) / sum(counts), mean = NA_real_, sd = NA_real_,
          p25 = NA_real_, median = NA_real_, p75 = NA_real_)
      }
    }
  }
  write.csv(do.call(rbind, rows), paste0(prefix, "_descriptives.csv"), row.names = FALSE)
  events <- list()
  for (outcome in cfg$outcomes) {
    z <- prepare_followup(d, outcome, cfg)
    for (g in names(groups)) {
      keep <- groups[[g]] & z$.stop > 0
      t <- z$.stop[keep] / 365.25
      events[[paste(outcome, g)]] <- data.frame(outcome = outcome, group = g, N = sum(keep),
        events = sum(z$.event[keep]), followup_mean_years = if (length(t)) mean(t) else NA_real_,
        followup_median_years = if (length(t)) median(t) else NA_real_)
    }
  }
  write.csv(do.call(rbind, events), paste0(prefix, "_outcomes.csv"), row.names = FALSE)
  for (v in c("no2", "ozone")) {
    cols <- intersect(c(paste0(v, c("_precon_mean", "_prenatal_mean", "_postnatal_mean")),
                        paste0(v, "_year_", 0:3)), names(d))
    x <- as.matrix(d[cols])
    write.csv(cor(x, method = "spearman", use = "pairwise.complete.obs"), paste0(prefix, "_", v, "_correlations.csv"))
    write.csv(crossprod(1L * !is.na(x)), paste0(prefix, "_", v, "_correlation_N.csv"))
  }
}

## 10) DEFINE TIME-VARYING IPCW AND WEIGHT DIAGNOSTICS

estimate_ipcw <- function(a,cfg) {
  d <- a$data
  if (anyNA(d$.censor)) stop("IPCW requires CFG$study_end and verified exit reasons for disenrollment weights. Check CFG$exit_reason_column and CFG$loss_to_followup_values. All-early-exit weights require an explicit CFG$censor_mode change.")
  basis <- make_bases(a,cfg)
  basis_names <- character()
  for(v in names(basis)) {
    b <- as.data.frame(unclass(basis[[v]]))
    names(b) <- paste0(".cens_",v,seq_len(ncol(b)))
    basis_names <- c(basis_names,names(b))
    d <- cbind(d,b)
  }
  cov <- c(a$covariates,basis_names)
  fit_censor <- function(terms) {
    f <- as.formula(paste("Surv(.start,.stop,.censor) ~",paste(c("strata(.stratum)",terms),collapse=" + ")))
    withCallingHandlers(
      survival::coxph(f,data=d,ties="efron",model=TRUE,na.action=na.fail,singular.ok=FALSE),
      warning=function(w) stop("Censoring model warning: ",conditionMessage(w)))
  }
  numerator <- denominator <- NULL
  if(any(d$.censor==1L)) {
    numerator <- fit_censor(character())
    denominator <- fit_censor(cov)
  }
  subjects <- d[c(".id",".stratum",".start",".stop",".event")]
  subjects$.position <- seq_len(nrow(d))
  subjects$.entry <- subjects$.start
  pieces <- split(subjects,subjects$.stratum,drop=TRUE)
  total <- 0
  for(k in seq_along(pieces)) {
    z <- pieces[[k]]; cuts <- sort(unique(z$.stop[z$.event==1]))
    total <- total + sum(1L + findInterval(z$.stop,cuts)-findInterval(z$.start,cuts))
    if(total > cfg$ipcw_max_rows) stop("IPCW expansion exceeds CFG$ipcw_max_rows. Run one outcome at a time or raise the limit after checking RAM.")
    if(length(cuts)) z <- survival::survSplit(Surv(.start,.stop,.event)~.,
      data=z,cut=cuts,start=".start",end=".stop",event=".event")
    z$.used <- z$.stop %in% cuts
    pieces[[k]] <- z
  }
  it <- do.call(rbind,pieces)
  if(sum(it$.event)!=sum(d$.event) || any(it$.stop<=it$.start)) stop("Invalid IPCW splitting.")
  log_G <- function(fit) {
    if(is.null(fit)) return(rep(0,nrow(it)))
    h <- survival::basehaz(fit,centered=TRUE)
    labels <- as.character(fit$model[["strata(.stratum)"]])
    if(!"strata" %in% names(h)) h$strata <- unique(labels)
    curves <- split(h,as.character(h$strata),drop=TRUE)
    rows <- split(seq_len(nrow(it)),labels[it$.position])
    lp <- if(length(coef(fit))) fit$linear.predictors else rep(0,nrow(d))
    if(any(!is.finite(lp)) || !all(names(rows) %in% names(curves))) stop("Censoring strata do not align.")
    ans <- numeric(nrow(it))
    for(s in names(rows)) {
      j <- rows[[s]]; curve <- curves[[s]]
      before <- findInterval(it$.stop[j],curve$time,left.open=TRUE)
      entry <- findInterval(it$.entry[j],curve$time)
      H <- pmax(0,c(0,curve$hazard)[before+1L]-c(0,curve$hazard)[entry+1L])
      ans[j] <- -exp(log(H)+lp[it$.position[j]])
    }
    if(any(!is.finite(ans))) stop("Invalid predicted censoring survival.")
    ans
  }
  log_den <- log_G(denominator)
  it$.weight <- exp(log_G(numerator)-log_den)
  if(any(!is.finite(it$.weight)|it$.weight<=0)) stop("Invalid IPCW; inspect positivity and model fit.")
  used <- it$.used
  if(!any(used)) stop("No outcome risk sets for IPCW.")
  w <- it$.weight[used]
  cap <- as.numeric(quantile(w,cfg$ipcw_truncate_quantile,names=FALSE))
  trunc <- it
  trunc$.weight <- pmin(trunc$.weight,cap)
  riskset <- paste(it$.stratum[used],it$.stop[used],sep=":")
  diagnostics_one <- function(weights,type) {
    x <- weights[used]
    scaled <- x/ave(x,riskset,FUN=max)
    sums <- rowsum(cbind(n=1,w=scaled,w2=scaled^2),riskset,reorder=FALSE)
    ess <- sums[,"w"]^2/sums[,"w2"]
    q <- quantile(x,c(0,.5,.95,.99,1),names=FALSE)
    data.frame(type=type,N=nrow(d),events=sum(d$.event),censor_events=sum(d$.censor),
      intervals=nrow(it),censor_definition=a$censor_definition,mean=mean(x),CV=sd(x)/mean(x),min=q[1],median=q[2],p95=q[3],p99=q[4],max=q[5],
      min_G=exp(min(log_den[used])),min_riskset_ESS_fraction=min(ess/sums[,"n"]),
      median_riskset_ESS_fraction=median(ess/sums[,"n"]),
      truncation_quantile=cfg$ipcw_truncate_quantile,truncation_cap=cap,
      fraction_capped=mean(w>cap))
  }
  cols <- c(".id",".start",".stop",".event",".weight")
  if(any(log_den[used]<log(.01))) warning("Some censoring survival probabilities are below 0.01; inspect IPCW diagnostics.")
  list(raw=it[cols],truncated=trunc[cols],diagnostics=rbind(
    diagnostics_one(it$.weight,"untruncated"),diagnostics_one(trunc$.weight,"truncated")))
}

## 11) DEFINE AIC COMPARISONS ON IDENTICAL ANALYTIC SAMPLES

compare_aic <- function(a, cfg) {
  candidates <- cfg$aic_lag_df[cfg$aic_lag_df <= ncol(a$X$ap)]
  grid <- expand.grid(ap = candidates, temp = candidates, rh = candidates)
  rows <- lapply(seq_len(nrow(grid)), function(i) {
    b <- a; b$df <- unlist(grid[i, ])
    warnings <- character()
    ans <- tryCatch(withCallingHandlers(fit_analysis(b, cfg, predictions = FALSE),
      warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") }),
      error = function(e) e)
    data.frame(grid[i, ], N = nrow(a$data), events = sum(a$data$.event),
      AIC = if (inherits(ans, "error")) NA_real_ else ans$AIC,
      error = if (inherits(ans, "error")) conditionMessage(ans) else "",
      warning = paste(unique(warnings), collapse = " | "))
  })
  out <- do.call(rbind, rows)
  if (!any(is.finite(out$AIC))) stop("All AIC candidates failed.")
  out$delta_AIC <- out$AIC - min(out$AIC, na.rm = TRUE)
  out
}

## 12) OPEN THE OUTPUT FOLDER AND DEFINE RESULT SAVING

stopifnot(length(CFG$outcomes) > 0L, !anyDuplicated(CFG$outcomes),
  all(CFG$outcomes %in% names(event_numbers)), all(CFG$pollutants %in% c("no2", "ozone")),
  !anyDuplicated(CFG$pollutants), CFG$main_prenatal_weeks == 37L,
  all(CFG$sensitivity_weeks %in% c(26L, 30L, 32L)), CFG$increment_ppb == 10,
  CFG$ipcw_truncate_quantile > 0, CFG$ipcw_truncate_quantile <= 1,
  all(CFG$ipcw_periods %in% c("preconception_prenatal", "postnatal")))
if (dir.exists(CFG$output_dir) && length(list.files(CFG$output_dir, all.files = TRUE, no.. = TRUE)))
  stop("Use an empty CFG$output_dir.")
dir.create(CFG$output_dir, recursive = TRUE, showWarnings = FALSE)
dput(CFG, file = file.path(CFG$output_dir, "settings_used.R"))
capture.output(sessionInfo(), file = file.path(CFG$output_dir, "sessionInfo.txt"))
message("Inputs must contain the validated outcome algorithms and completed cohort exclusions, including genetic/chromosomal exclusions.")
message("Annual postnatal models retain birth-origin follow-up and fixed years 0-3 exposure histories; some exposures may follow an early diagnosis.")
if (is.na(CFG$study_end)) message("CFG$study_end is unset: main models use prepared enrollment endpoints; IPCW cannot run until the administrative end is supplied.")
if (CFG$censor_mode == "all_early_exits") message("IPCW targets all nonadministrative exits, including death when it is not separated; use exit reasons for disenrollment-only weights.")
if (CFG$censor_mode == "disenrollment" && is.null(CFG$exit_reason_column))
  message("Set CFG$exit_reason_column and CFG$loss_to_followup_values for disenrollment IPCW; main models can run without these fields.")
results <- list(); analysis_status <- list()

save_status <- function() {
  if (length(results)) write.csv(do.call(rbind, results), file.path(CFG$output_dir, "all_HRs.csv"), row.names = FALSE)
  if (length(analysis_status)) write.csv(do.call(rbind, analysis_status),
    file.path(CFG$output_dir, "analysis_status.csv"), row.names = FALSE)
}

run_step <- function(label, code) {
  message(label)
  warnings <- character()
  ans <- tryCatch(withCallingHandlers(force(code), warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
  }), error = function(e) e)
  failed <- inherits(ans, "error")
  analysis_status[[length(analysis_status) + 1L]] <<- data.frame(analysis = label,
    status = if (failed) "FAILED" else if (length(warnings)) "WARNING" else "OK",
    message = if (failed) conditionMessage(ans) else paste(unique(warnings), collapse = " | "))
  save_status()
  if (failed) message("FAILED: ", label, ": ", conditionMessage(ans))
  invisible(ans)
}

save_fit <- function(a, label, intervals = NULL) {
  fit <- fit_analysis(a, CFG, intervals)
  table <- result_table(fit, a, label)
  prefix <- file.path(CFG$output_dir, gsub("[^A-Za-z0-9_.-]", "_", label))
  write.csv(table, paste0(prefix, "_HR.csv"), row.names = FALSE)
  write.csv(fit$lag_HRs, paste0(prefix, "_lag_HRs.csv"), row.names = FALSE)
  write.csv(data.frame(input_n = a$input_n, analyzed_n = fit$N, events = fit$events,
    excluded_n = a$excluded_n, dropped_constant = paste(a$constant, collapse = ";")),
    paste0(prefix, "_sample.csv"), row.names = FALSE)
  dput(list(formula = fit$formula, basis = fit$basis, copollutants = a$copollutants,
    post_adjust_early_copollutants = a$post_adjust_early_copollutants), file = paste0(prefix, "_model_specification.R"))
  saveRDS(fit, paste0(prefix, "_predictions.rds"))
  if (CFG$save_plots) save_plots(fit, a, prefix)
  results[[length(results) + 1L]] <<- table
  invisible(fit)
}

## 13) PREPARE THE MAIN COHORT, TERTILES AND DESCRIPTIVE TABLES

prepared <- run_step("main_cohort_preparation", {
  d <- prepare_cohort(read_cohort("main", CFG), CFG)
  ex <- exposure_matrices(d, CFG$main_prenatal_weeks, CFG)
  cuts <- CFG$population_density_cutpoints
  if (is.null(cuts)) {
    x <- d$pop_density_30sec_NASA[complete_early(ex)]
    if (!any(is.finite(x))) stop("Population density is unavailable in the main exposure cohort.")
    cuts <- unname(quantile(x, c(1/3, 2/3), na.rm = TRUE))
  }
  d <- prepare_cohort(d, CFG, cuts)
  write.csv(data.frame(lower_tertile_cutpoint = cuts[1], upper_tertile_cutpoint = cuts[2]),
    file.path(CFG$output_dir, "population_density_cutpoints.csv"), row.names = FALSE)
  write.csv(ex$audit, file.path(CFG$output_dir, "main_exposure_mapping.csv"), row.names = FALSE)
  list(data = d, exposure = ex, density_cutpoints = cuts)
})
if (inherits(prepared, "error")) stop("Main preparation failed; inspect analysis_status.csv.")
main_data <- prepared$data
main_exposure <- prepared$exposure
density_cutpoints <- prepared$density_cutpoints
rm(prepared)
run_step("main_descriptives", save_descriptives(main_data, main_exposure, CFG, file.path(CFG$output_dir, "main")))

for (outcome in CFG$outcomes) {
  followup <- run_step(paste0(outcome, "_followup"), {
    d <- prepare_followup(main_data, outcome, CFG)
    write.csv(attr(d, "followup_audit"), file.path(CFG$output_dir, paste0(outcome, "_followup_audit.csv")), row.names = FALSE)
    d
  })
  if (inherits(followup, "error")) next
  for (pollutant in CFG$pollutants) for (period in c("preconception_prenatal", "postnatal")) {

    ## 14) FIT THE ADJUSTED PRECONCEPTION, PRENATAL AND ANNUAL POSTNATAL MODELS

    label <- paste(outcome, pollutant, period, sep = "_")
    a <- run_step(paste0(label, "_prepare"),
      make_analysis(followup, main_exposure, outcome, pollutant, period, CFG$main_prenatal_weeks, CFG))
    if (inherits(a, "error")) next
    run_step(paste0(label, "_main_adjusted"), save_fit(a, paste0(label, "_main_adjusted")))

    ## 15) FIT SEX, RACE AND POPULATION-DENSITY SUBGROUP MODELS

    if (CFG$run_subgroups) for (group in names(CFG$subgroup_levels)) {
      for (value in CFG$subgroup_levels[[group]]) {
        group_label <- paste(label, group, value, sep = "_")
        run_step(group_label, {
          b <- subset_analysis(a, group, value)
          save_fit(b, group_label)
        })
      }
    }

    ## 16) COMPARE UNWEIGHTED AND IPCW MODELS ON THE SAME SAMPLE

    if (CFG$run_ipcw && period %in% CFG$ipcw_periods) run_step(paste0(label, "_IPCW"), {
      weights <- estimate_ipcw(a, CFG)
      write.csv(weights$diagnostics, file.path(CFG$output_dir, paste0(label, "_IPCW_diagnostics.csv")), row.names = FALSE)
      save_fit(a, paste0(label, "_IPCW_same_sample_unweighted"))
      save_fit(a, paste0(label, "_IPCW_untruncated"), weights$raw)
      save_fit(a, paste0(label, "_IPCW_truncated"), weights$truncated)
      if (CFG$save_plots) {
        pdf(file.path(CFG$output_dir, paste0(label, "_IPCW_weights.pdf")), width = 8, height = 4)
        par(mfrow = c(1, 2))
        hist(log10(weights$raw$.weight), breaks = 40, main = "Untruncated IPCW", xlab = "log10(weight)")
        hist(log10(weights$truncated$.weight), breaks = 40, main = "Truncated IPCW", xlab = "log10(weight)")
        dev.off()
      }
      rm(weights)
    })

    ## 17) OPTIONALLY COMPARE LAG DEGREES OF FREEDOM USING AIC

    if (CFG$run_aic) run_step(paste0(label, "_AIC"), {
      write.csv(compare_aic(a, CFG), file.path(CFG$output_dir, paste0(label, "_AIC.csv")), row.names = FALSE)
    })
    rm(a)
  }
  rm(followup)
  gc(verbose = FALSE)
}

## 18) FIT THE 26-, 30- AND 32-WEEK PRENATAL SENSITIVITY ANALYSES

if (CFG$run_sensitivity) for (npren in CFG$sensitivity_weeks) {
  key <- as.character(npren)
  sens <- run_step(paste0("cohort_", npren, "wk_preparation"), {
    d <- prepare_cohort(read_cohort(key, CFG), CFG, density_cutpoints)
    ex <- exposure_matrices(d, npren, CFG, include_post = FALSE)
    write.csv(ex$audit, file.path(CFG$output_dir, paste0("exposure_mapping_", npren, "wk.csv")), row.names = FALSE)
    list(data = d, exposure = ex)
  })
  if (inherits(sens, "error")) next
  for (outcome in CFG$outcomes) {
    followup <- run_step(paste(outcome, npren, "wk_followup", sep = "_"), prepare_followup(sens$data, outcome, CFG))
    if (inherits(followup, "error")) next
    for (pollutant in CFG$pollutants) {
      label <- paste(outcome, pollutant, "preconception_prenatal", npren, "wk_adjusted", sep = "_")
      run_step(label, {
        a <- make_analysis(followup, sens$exposure, outcome, pollutant,
          "preconception_prenatal", npren, CFG)
        save_fit(a, label)
      })
    }
  }
  rm(sens)
  gc(verbose = FALSE)
}

## 19) SAVE COMBINED RESULTS AND IDENTIFY ANY FAILED ANALYSES

save_status()
status <- do.call(rbind, analysis_status)
if (any(status$status == "FAILED"))
  stop("Some requested analyses failed. Successful results were saved; inspect analysis_status.csv.")
if (any(status$status == "WARNING")) warning("Inspect warnings in analysis_status.csv before interpreting results.")
message("Analysis finished: ", normalizePath(CFG$output_dir))

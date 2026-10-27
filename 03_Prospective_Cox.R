==============================================================================
# Script 3 (Expanded): Prospective Cox Proportional Hazards Regression
# Description: Incident psoriasis was modeled as a time-varying exposure using time-dependent Cox regression to estimate its association with four major categories of respiratory tract infection, 
each individual respiratory infection, and pneumonia-specific subgroups, with interaction assessed by likelihood ratio tests.
==============================================================================
###########Time-varying Cox regression for four major categories of respiratory infection
library(survival)

# 1. Load datasets
data2 <- read.csv("total_data_new_407208.csv")
followup_time <- read.csv("20250613data_followuptime(50w).csv")

# 2. Merge official follow-up time (124th column)
fu_colname <- names(followup_time)[124]
data2 <- merge(
  data2,
  followup_time[, c("eid", fu_colname)],
  by = "eid",
  all.x = TRUE,
  sort = FALSE
)
names(data2)[names(data2) == fu_colname] <- "fu_time"

# 3. Define key column mappings and covariates
date_baseline_col <- "p53_i0"
date_pso_col      <- "p131742"
cov_cols          <- 25:53

outcome_categories <- list(
  "Pneumonia"       = 4:10,
  "Lower_RI"        = 11:13,
  "Upper_RI"        = 14:20,
  "Influenza"       = 21:23
)

# Convert baseline and psoriasis dates
base_date <- as.Date(trimws(as.character(data2[[date_baseline_col]])), "%Y-%m-%d")
pso_date  <- as.Date(trimws(as.character(data2[[date_pso_col]])), "%Y-%m-%d")

if (!is.numeric(data2$fu_time)) {
  data2$fu_time <- as.numeric(as.character(data2$fu_time))
}

for (nm in names(data2)[cov_cols]) {
  if (is.character(data2[[nm]])) data2[[nm]] <- as.factor(data2[[nm]])
}
cov_names <- names(data2)[cov_cols]

# 4. Define robust function with Two-Step tmerge
run_time_varying_cox <- function(df_full, outcome_name, outcome_cols) {
  cat("\n>>> Processing outcome:", outcome_name, "\n")
  
  out_df <- df_full[, outcome_cols, drop = FALSE]
  out_df_dates <- data.frame(lapply(out_df, function(col) {
    as.Date(trimws(as.character(col)), "%Y-%m-%d")
  }), stringsAsFactors = FALSE)
  
  earliest_out_date <- do.call(pmin, c(out_df_dates, list(na.rm = TRUE)))
  
  d_sub <- data.frame(
    eid = df_full$eid,
    base_date = base_date,
    pso_date = pso_date,
    out_date = earliest_out_date,
    fu_time = df_full$fu_time,
    stringsAsFactors = FALSE
  )
  
  for (nm in cov_names) {
    d_sub[[nm]] <- df_full[[nm]]
  }
  
  # Exclude prevalent outcome cases
  has_prevalent_out <- !is.na(d_sub$out_date) & d_sub$out_date <= d_sub$base_date
  d_sub <- d_sub[!has_prevalent_out, ]
  
  # Exclude prevalent psoriasis cases at baseline
  has_prevalent_pso <- !is.na(d_sub$pso_date) & d_sub$pso_date <= d_sub$base_date
  d_sub <- d_sub[!has_prevalent_pso, ]
  
  # Calculate times
  d_sub$t_out <- as.numeric(d_sub$out_date - d_sub$base_date) / 365.25
  d_sub$t_end <- pmin(d_sub$t_out, d_sub$fu_time, na.rm = TRUE)
  d_sub <- d_sub[is.finite(d_sub$t_end) & d_sub$t_end > 0, ]
  
  d_sub$event <- ifelse(!is.na(d_sub$t_out) & d_sub$t_out <= d_sub$fu_time, 1L, 0L)
  
  cat("Valid sample size (n):", nrow(d_sub), " | Events:", sum(d_sub$event == 1, na.rm = TRUE), "\n")
  
  if (sum(d_sub$event == 1, na.rm = TRUE) < 5) {
    cat("Skipped: Events less than 5.\n")
    return(NULL)
  }
  
  # Incident psoriasis transition time
  d_sub$t_pso_incident <- ifelse(
    !is.na(d_sub$pso_date) & d_sub$pso_date > d_sub$base_date,
    as.numeric(d_sub$pso_date - d_sub$base_date) / 365.25,
    NA
  )
  d_sub$t_pso_incident[!is.na(d_sub$t_pso_incident) & d_sub$t_pso_incident >= d_sub$t_end] <- NA
  
  # Valid covariates
  valid_covs <- Filter(function(nm) {
    length(unique(na.omit(d_sub[[nm]]))) >= 2
  }, cov_names)
  cov_terms <- paste(valid_covs, collapse = " + ")
  
  # --- STEP 1: Initialize tmerge structure with outcome endpoint ---
  cat("Building tmerge step 1 (endpoint)...\n")
  base_subset <- d_sub[, c("eid", "t_end", "event", valid_covs)]
  td_data <- tmerge(base_subset, base_subset, id = eid, death = event(t_end, event))
  
  # --- STEP 2: Add time-varying exposure (tdc) ---
  cat("Building tmerge step 2 (time-varying exposure)...\n")
  pso_subset <- d_sub[!is.na(d_sub$t_pso_incident), c("eid", "t_pso_incident")]
  if (nrow(pso_subset) > 0) {
    td_data <- tmerge(td_data, pso_subset, id = eid, psoriasis_tv = tdc(t_pso_incident, init = 0))
  } else {
    td_data$psoriasis_tv <- 0
  }
  
  # --- 1. Adjusted Model ---
  cat("Fitting Adjusted Time-Varying Cox Model...\n")
  fml_adj <- as.formula(paste0("Surv(tstart, tstop, death) ~ psoriasis_tv + ", cov_terms))
  m_adj <- coxph(fml_adj, data = td_data)
  s_adj <- summary(m_adj)
  
  beta_adj <- s_adj$coef["psoriasis_tv", "coef"]
  se_adj   <- s_adj$coef["psoriasis_tv", "se(coef)"]
  pval_adj <- s_adj$coef["psoriasis_tv", "Pr(>|z|)"]
  
  ph_p_exp <- NA; ph_p_glob <- NA
  zph <- tryCatch(cox.zph(m_adj), error = function(e) NULL)
  if (!is.null(zph)) {
    if ("psoriasis_tv" %in% rownames(zph$table)) ph_p_exp <- zph$table["psoriasis_tv", "p"]
    if ("GLOBAL" %in% rownames(zph$table)) ph_p_glob <- zph$table["GLOBAL", "p"]
  }
  
  res_adj <- data.frame(
    Outcome = outcome_name,
    Model = "Adjusted",
    HR = exp(beta_adj),
    CI_low = exp(beta_adj - 1.96 * se_adj),
    CI_high = exp(beta_adj + 1.96 * se_adj),
    p = pval_adj,
    PH_p_exposure = ph_p_exp,
    PH_p_global = ph_p_glob,
    n = m_adj$n,
    events = m_adj$nevent,
    stringsAsFactors = FALSE
  )
  
  # --- 2. Unadjusted Model ---
  cat("Fitting Unadjusted Time-Varying Cox Model...\n")
  fml_unadj <- Surv(tstart, tstop, death) ~ psoriasis_tv
  m_unadj <- coxph(fml_unadj, data = td_data)
  s_unadj <- summary(m_unadj)
  
  beta_unadj <- s_unadj$coef["psoriasis_tv", "coef"]
  se_unadj   <- s_unadj$coef["psoriasis_tv", "se(coef)"]
  pval_unadj <- s_unadj$coef["psoriasis_tv", "Pr(>|z|)"]
  
  res_unadj <- data.frame(
    Outcome = outcome_name,
    Model = "Unadjusted",
    HR = exp(beta_unadj),
    CI_low = exp(beta_unadj - 1.96 * se_unadj),
    CI_high = exp(beta_unadj + 1.96 * se_unadj),
    p = pval_unadj,
    PH_p_exposure = NA,
    PH_p_global = NA,
    n = m_unadj$n,
    events = m_unadj$nevent,
    stringsAsFactors = FALSE
  )
  
  rbind(res_unadj, res_adj)
}

# 5. Batch loop
results_list <- list()
for (cat_name in names(outcome_categories)) {
  cols <- outcome_categories[[cat_name]]
  res <- run_time_varying_cox(data2, cat_name, cols)
  if (!is.null(res)) {
    results_list[[length(results_list) + 1]] <- res
  }
}

# 6. Combine and export results
if (length(results_list) > 0) {
  final_tv_summary <- do.call(rbind, results_list)
  write.csv(final_tv_summary, "cox_time_varying_respiratory_summary.csv", row.names = FALSE)
  cat("\nSUCCESS! Results successfully exported to cox_time_varying_respiratory_summary.csv\n")
} else {
  cat("\nNOTICE: results_list is empty.\n")
}



###########Association of incident psoriasis with each respiratory infection using time-varying Cox regression
library(survival)

# 1. Load datasets
data2 <- read.csv("total_data_new_407208.csv")
followup_time <- read.csv("20250613data_followuptime(50w).csv")

# 2. Merge official follow-up time (124th column)
fu_colname <- names(followup_time)[124]
data2 <- merge(
  data2,
  followup_time[, c("eid", fu_colname)],
  by = "eid",
  all.x = TRUE,
  sort = FALSE
)
names(data2)[names(data2) == fu_colname] <- "fu_time"

# 3. Define key column mappings and covariates
date_baseline_col <- "p53_i0"
date_pso_col      <- "p131742"
cov_cols          <- 25:53

# Convert baseline and psoriasis dates
base_date <- as.Date(trimws(as.character(data2[[date_baseline_col]])), "%Y-%m-%d")
pso_date  <- as.Date(trimws(as.character(data2[[date_pso_col]])), "%Y-%m-%d")

if (!is.numeric(data2$fu_time)) {
  data2$fu_time <- as.numeric(as.character(data2$fu_time))
}

for (nm in names(data2)[cov_cols]) {
  if (is.character(data2[[nm]])) data2[[nm]] <- as.factor(data2[[nm]])
}
cov_names <- names(data2)[cov_cols]

# 4. Universal function for single outcome column analysis
run_time_varying_cox_single <- function(df_full, col_index) {
  outcome_name <- names(df_full)[col_index]
  cat("\n>>> Processing column:", col_index, "(", outcome_name, ")\n")
  
  # Extract single outcome date column
  out_date_raw <- as.Date(trimws(as.character(df_full[[col_index]])), "%Y-%m-%d")
  
  d_sub <- data.frame(
    eid = df_full$eid,
    base_date = base_date,
    pso_date = pso_date,
    out_date = out_date_raw,
    fu_time = df_full$fu_time,
    stringsAsFactors = FALSE
  )
  
  for (nm in cov_names) {
    d_sub[[nm]] <- df_full[[nm]]
  }
  
  # Exclude prevalent outcome cases
  has_prevalent_out <- !is.na(d_sub$out_date) & d_sub$out_date <= d_sub$base_date
  d_sub <- d_sub[!has_prevalent_out, ]
  
  # Exclude prevalent psoriasis cases at baseline
  has_prevalent_pso <- !is.na(d_sub$pso_date) & d_sub$pso_date <= d_sub$base_date
  d_sub <- d_sub[!has_prevalent_pso, ]
  
  # Calculate times
  d_sub$t_out <- as.numeric(d_sub$out_date - d_sub$base_date) / 365.25
  d_sub$t_end <- pmin(d_sub$t_out, d_sub$fu_time, na.rm = TRUE)
  d_sub <- d_sub[is.finite(d_sub$t_end) & d_sub$t_end > 0, ]
  
  d_sub$event <- ifelse(!is.na(d_sub$t_out) & d_sub$t_out <= d_sub$fu_time, 1L, 0L)
  
  cat("Valid sample size (n):", nrow(d_sub), " | Events:", sum(d_sub$event == 1, na.rm = TRUE), "\n")
  
  # Skip if event count is less than 5
  if (sum(d_sub$event == 1, na.rm = TRUE) < 5) {
    cat("Skipped: Events less than 5.\n")
    return(NULL)
  }
  
  # Incident psoriasis transition time
  d_sub$t_pso_incident <- ifelse(
    !is.na(d_sub$pso_date) & d_sub$pso_date > d_sub$base_date,
    as.numeric(d_sub$pso_date - d_sub$base_date) / 365.25,
    NA
  )
  d_sub$t_pso_incident[!is.na(d_sub$t_pso_incident) & d_sub$t_pso_incident >= d_sub$t_end] <- NA
  
  # Valid covariates
  valid_covs <- Filter(function(nm) {
    length(unique(na.omit(d_sub[[nm]]))) >= 2
  }, cov_names)
  cov_terms <- paste(valid_covs, collapse = " + ")
  
  # --- STEP 1: Initialize tmerge structure with outcome endpoint ---
  base_subset <- d_sub[, c("eid", "t_end", "event", valid_covs)]
  td_data <- tmerge(base_subset, base_subset, id = eid, death = event(t_end, event))
  
  # --- STEP 2: Add time-varying exposure (tdc) ---
  pso_subset <- d_sub[!is.na(d_sub$t_pso_incident), c("eid", "t_pso_incident")]
  if (nrow(pso_subset) > 0) {
    td_data <- tmerge(td_data, pso_subset, id = eid, psoriasis_tv = tdc(t_pso_incident, init = 0))
  } else {
    td_data$psoriasis_tv <- 0
  }
  
  # --- 1. Adjusted Model ---
  fml_adj <- as.formula(paste0("Surv(tstart, tstop, death) ~ psoriasis_tv + ", cov_terms))
  m_adj <- tryCatch(coxph(fml_adj, data = td_data), error = function(e) NULL)
  if (is.null(m_adj)) return(NULL)
  s_adj <- summary(m_adj)
  
  if (!("psoriasis_tv" %in% rownames(s_adj$coef))) return(NULL)
  
  beta_adj <- s_adj$coef["psoriasis_tv", "coef"]
  se_adj   <- s_adj$coef["psoriasis_tv", "se(coef)"]
  pval_adj <- s_adj$coef["psoriasis_tv", "Pr(>|z|)"]
  
  ph_p_exp <- NA; ph_p_glob <- NA
  zph <- tryCatch(cox.zph(m_adj), error = function(e) NULL)
  if (!is.null(zph)) {
    if ("psoriasis_tv" %in% rownames(zph$table)) ph_p_exp <- zph$table["psoriasis_tv", "p"]
    if ("GLOBAL" %in% rownames(zph$table)) ph_p_glob <- zph$table["GLOBAL", "p"]
  }
  
  res_adj <- data.frame(
    Column_Index = col_index,
    Outcome_Name = outcome_name,
    Model = "Adjusted",
    HR = exp(beta_adj),
    CI_low = exp(beta_adj - 1.96 * se_adj),
    CI_high = exp(beta_adj + 1.96 * se_adj),
    p = pval_adj,
    PH_p_exposure = ph_p_exp,
    PH_p_global = ph_p_glob,
    n = m_adj$n,
    events = m_adj$nevent,
    stringsAsFactors = FALSE
  )
  
  # --- 2. Unadjusted Model ---
  fml_unadj <- Surv(tstart, tstop, death) ~ psoriasis_tv
  m_unadj <- tryCatch(coxph(fml_unadj, data = td_data), error = function(e) NULL)
  if (is.null(m_unadj)) return(res_adj)
  s_unadj <- summary(m_unadj)
  
  beta_unadj <- s_unadj$coef["psoriasis_tv", "coef"]
  se_unadj   <- s_unadj$coef["psoriasis_tv", "se(coef)"]
  pval_unadj <- s_unadj$coef["psoriasis_tv", "Pr(>|z|)"]
  
  res_unadj <- data.frame(
    Column_Index = col_index,
    Outcome_Name = outcome_name,
    Model = "Unadjusted",
    HR = exp(beta_unadj),
    CI_low = exp(beta_unadj - 1.96 * se_unadj),
    CI_high = exp(beta_unadj + 1.96 * se_unadj),
    p = pval_unadj,
    PH_p_exposure = NA,
    PH_p_global = NA,
    n = m_unadj$n,
    events = m_unadj$nevent,
    stringsAsFactors = FALSE
  )
  
  rbind(res_unadj, res_adj)
}

# 5. Batch loop across columns 4 to 23 individually
results_list2 <- list()
for (col_idx in 4:23) {
  res <- run_time_varying_cox_single(data2, col_idx)
  if (!is.null(res)) {
    results_list2[[length(results_list2) + 1]] <- res
  }
}

# 6. Combine and export results
if (length(results_list2) > 0) {
  final_tv_summary <- do.call(rbind, results_list2)
  write.csv(final_tv_summary, "cox_time_varying_individual_columns_summary2.csv", row.names = FALSE)
  cat("\nSUCCESS! Results successfully exported to cox_time_varying_individual_columns_summary.csv\n")
} else {
  cat("\nNOTICE: results_list is empty.\n")
}



###########Stratified and multiplicative/additive interaction analyses of incident psoriasis and risk of pneumonia using time-varying Cox regression
library(survival)

# 1. Load datasets
data2 <- read.csv("total_data_new_407208.csv")
followup_time <- read.csv("20250613data_followuptime(50w).csv")

# 2. Merge official follow-up time (124th column)
fu_colname <- names(followup_time)[124]
data2 <- merge(
  data2,
  followup_time[, c("eid", fu_colname)],
  by = "eid",
  all.x = TRUE,
  sort = FALSE
)
names(data2)[names(data2) == fu_colname] <- "fu_time"

# 3. Define key column mappings and covariates
date_baseline_col <- "p53_i0"
date_pso_col      <- "p131742"
cov_cols          <- 25:53

# Pneumonia outcome category (columns 4 to 10)
pneumonia_cols    <- 4:10

# Convert baseline and psoriasis dates
base_date <- as.Date(trimws(as.character(data2[[date_baseline_col]])), "%Y-%m-%d")
pso_date  <- as.Date(trimws(as.character(data2[[date_pso_col]])), "%Y-%m-%d")

if (!is.numeric(data2$fu_time)) {
  data2$fu_time <- as.numeric(as.character(data2$fu_time))
}

# Construct derived subgroup variables on raw data before tmerge
data2$age_group <- ifelse(data2$age >= 60, ">=60", "<60")
data2$bmi_group <- ifelse(data2$p21001_i0 >= 30, "Obese_geq30", "NonObese_lt30")

for (nm in names(data2)[cov_cols]) {
  if (is.character(data2[[nm]])) data2[[nm]] <- as.factor(data2[[nm]])
}
cov_names <- names(data2)[cov_cols]

# Ensure subgroup factors
data2$p31       <- as.factor(data2$p31)
data2$age_group <- as.factor(data2$age_group)
data2$bmi_group <- as.factor(data2$bmi_group)
data2$p20116_i0 <- as.factor(data2$p20116_i0)

cat("\n>>> Processing Pneumonia Outcome (Columns 4-10) for Stratification & Interaction <<<\n")

# Extract outcome date columns
out_df <- data2[, pneumonia_cols, drop = FALSE]
out_df_dates <- data.frame(lapply(out_df, function(col) {
  as.Date(trimws(as.character(col)), "%Y-%m-%d")
}), stringsAsFactors = FALSE)

earliest_out_date <- do.call(pmin, c(out_df_dates, list(na.rm = TRUE)))

d_sub <- data.frame(
  eid = data2$eid,
  base_date = base_date,
  pso_date = pso_date,
  out_date = earliest_out_date,
  fu_time = data2$fu_time,
  p31 = data2$p31,
  age_group = data2$age_group,
  bmi_group = data2$bmi_group,
  smoking = data2$p20116_i0,
  stringsAsFactors = FALSE
)

for (nm in cov_names) {
  d_sub[[nm]] <- data2[[nm]]
}

# Exclude prevalent outcome cases & prevalent psoriasis cases at baseline
has_prevalent_out <- !is.na(d_sub$out_date) & d_sub$out_date <= d_sub$base_date
has_prevalent_pso <- !is.na(d_sub$pso_date) & d_sub$pso_date <= d_sub$base_date
d_sub <- d_sub[!(has_prevalent_out | has_prevalent_pso), ]

# Calculate times
d_sub$t_out <- as.numeric(d_sub$out_date - d_sub$base_date) / 365.25
d_sub$t_end <- pmin(d_sub$t_out, d_sub$fu_time, na.rm = TRUE)
d_sub <- d_sub[is.finite(d_sub$t_end) & d_sub$t_end > 0, ]
d_sub$event <- ifelse(!is.na(d_sub$t_out) & d_sub$t_out <= d_sub$fu_time, 1L, 0L)

cat("Clean cohort sample size (n):", nrow(d_sub), " | Events:", sum(d_sub$event == 1, na.rm = TRUE), "\n")

# Incident psoriasis transition time
d_sub$t_pso_incident <- ifelse(
  !is.na(d_sub$pso_date) & d_sub$pso_date > d_sub$base_date,
  as.numeric(d_sub$pso_date - d_sub$base_date) / 365.25,
  NA
)
d_sub$t_pso_incident[!is.na(d_sub$t_pso_incident) & d_sub$t_pso_incident >= d_sub$t_end] <- NA

# Valid covariates (excluding active subgroup variables to avoid collinearity)
valid_covs <- Filter(function(nm) {
  length(unique(na.omit(d_sub[[nm]]))) >= 2
}, cov_names)

# --- TWO-STEP TMERGE ---
cat("Building tmerge structure...\n")
subgroup_vars <- c("p31", "age_group", "bmi_group", "smoking")
base_vars <- unique(c("eid", "t_end", "event", subgroup_vars, valid_covs))
base_subset <- d_sub[, base_vars]

td_data <- tmerge(base_subset, base_subset, id = eid, death = event(t_end, event))

pso_subset <- d_sub[!is.na(d_sub$t_pso_incident), c("eid", "t_pso_incident")]
if (nrow(pso_subset) > 0) {
  td_data <- tmerge(td_data, pso_subset, id = eid, psoriasis_tv = tdc(t_pso_incident, init = 0))
} else {
  td_data$psoriasis_tv <- 0
}

# --- FUNCTION FOR STRATIFIED & INTERACTION ANALYSIS ---
run_subgroup_analysis <- function(td_df, strat_var_name, strat_col_name, covs_vector) {
  cat("\n--- Running Subgroup Analysis for:", strat_var_name, "---\n")
  
  # Remove strat_var from covariates if present
  current_covs <- setdiff(covs_vector, strat_col_name)
  cov_terms <- paste(current_covs, collapse = " + ")
  
  strat_levels <- levels(td_df[[strat_col_name]])
  strat_results <- list()
  
  # 1. Stratum-specific models
  for (lvl in strat_levels) {
    sub_data <- td_df[td_df[[strat_col_name]] == lvl, ]
    fml_strat <- as.formula(paste0("Surv(tstart, tstop, death) ~ psoriasis_tv + ", cov_terms))
    
    res_strat <- tryCatch({
      m <- coxph(fml_strat, data = sub_data)
      s <- summary(m)
      if (!("psoriasis_tv" %in% rownames(s$coef))) stop("Exposure dropped")
      
      beta <- s$coef["psoriasis_tv", "coef"]
      se   <- s$coef["psoriasis_tv", "se(coef)"]
      
      data.frame(
        Variable = strat_var_name,
        Subgroup_Level = lvl,
        Analysis_Type = "Stratum-Specific",
        HR = exp(beta),
        CI_low = exp(beta - 1.96 * se),
        CI_high = exp(beta + 1.96 * se),
        p = s$coef["psoriasis_tv", "Pr(>|z|)"],
        Interaction_P = NA,
        n = m$n,
        events = m$nevent,
        stringsAsFactors = FALSE
      )
    }, error = function(e) { NULL })
    
    if (!is.null(res_strat)) strat_results[[length(strat_results) + 1]] <- res_strat
  }
  
  # 2. Formal Interaction Model
  fml_interact <- as.formula(paste0("Surv(tstart, tstop, death) ~ psoriasis_tv * ", strat_col_name, " + ", cov_terms))
  interact_res <- tryCatch({
    m_int <- coxph(fml_interact, data = td_df)
    s_int <- summary(m_int)
    
    # Overall likelihood ratio test or Wald test p-value for interaction
    anova_res <- anova(m_int)
    # Find interaction term row in anova or summary
    int_p_val <- NA
    match_rows <- grep(paste0("psoriasis_tv:", strat_col_name), rownames(s_int$coef))
    if (length(match_rows) > 0) {
      # Extract omnibus interaction p-value using drop1 or anova
      drop_res <- drop1(m_int, test = "Chisq")
      int_row_name <- grep(paste0("psoriasis_tv:", strat_col_name), rownames(drop_res))
      if (length(int_row_name) > 0) int_p_val <- drop_res[int_row_name, "Pr(>Chi)"]
    }
    
    # Return interaction summary row representing the global test of interaction
    data.frame(
      Variable = strat_var_name,
      Subgroup_Level = "Overall Interaction",
      Analysis_Type = "Interaction Test",
      HR = NA,
      CI_low = NA,
      CI_high = NA,
      p = NA,
      Interaction_P = int_p_val,
      n = m_int$n,
      events = m_int$nevent,
      stringsAsFactors = FALSE
    )
  }, error = function(e) { NULL })
  
  res_df <- do.call(rbind, strat_results)
  if (!is.null(interact_res)) res_df <- rbind(res_df, interact_res)
  res_df
}

# Run across all 4 subgroup variables
subgroup_specs <- list(
  list(name = "Gender (p31)", col = "p31"),
  list(name = "Age Group (<60 vs >=60)", col = "age_group"),
  list(name = "BMI Group (<30 vs >=30)", col = "bmi_group"),
  list(name = "Smoking Status (p20116_i0)", col = "smoking")
)

all_subgroup_results <- list()
for (spec in subgroup_specs) {
  res_sub <- run_subgroup_analysis(td_data, spec$name, spec$col, valid_covs)
  if (!is.null(res_sub)) {
    all_subgroup_results[[length(all_subgroup_results) + 1]] <- res_sub
  }
}

# Export combined subgroup and interaction summary
if (length(all_subgroup_results) > 0) {
  final_strat_summary <- do.call(rbind, all_subgroup_results)
  write.csv(final_strat_summary, "pneumonia_stratified_and_interaction_summary3.csv", row.names = FALSE)
  cat("\nSUCCESS! Stratified and Interaction results exported to pneumonia_stratified_and_interaction_summary.csv\n")
} else {
  cat("\nNOTICE: No results generated.\n")
}


library(epiR)
td_data$smoking_binary <- ifelse(
  as.character(td_data$smoking) %in% c("Never", "0", "Never smoked"), "Never", "Ever"
)

td_data$smoking_binary <- as.factor(td_data$smoking_binary)
td_data$smoking_binary <- relevel(td_data$smoking_binary, ref = "Never")

td_data$p31       <- as.factor(td_data$p31)       # sex
td_data$age_group <- as.factor(td_data$age_group) # "<60" vs ">=60"
td_data$bmi_group <- as.factor(td_data$bmi_group) # "NonObese_lt30" vs "Obese_geq30"

run_additive_interaction_with_ci <- function(td_df, strat_var_name, strat_col_name, covs_vector) {
  cat("\n--- Running Additive Interaction with 95% CI for:", strat_var_name, "---\n")
  
  current_covs <- setdiff(covs_vector, c(strat_col_name, "smoking"))
  cov_terms <- paste(current_covs, collapse = " + ")
  
  strat_levels <- levels(td_df[[strat_col_name]])
  if (length(strat_levels) != 2) {
    return(NULL)
  }
  
  fml_add <- as.formula(paste0("Surv(tstart, tstop, death) ~ psoriasis_tv * ", strat_col_name, " + ", cov_terms))
  m_add <- tryCatch({
    coxph(fml_add, data = td_df)
  }, error = function(e) { NULL })
  
  if (is.null(m_add)) return(NULL)
  
  coefs <- coef(m_add)
  vcv <- vcov(m_add)
  
  term_pso <- "psoriasis_tv"
  term_strat <- grep(strat_col_name, names(coefs), value = TRUE)[1]
  term_interact <- grep(paste0("psoriasis_tv:", strat_col_name), names(coefs), value = TRUE)
  
  if (length(term_interact) == 0) return(NULL)
  
  b1  <- coefs[term_pso]
  b2  <- coefs[term_strat]
  b12 <- coefs[term_interact]
  
  RR10 <- exp(b1)
  RR01 <- exp(b2)
  RR11 <- exp(b1 + b2 + b12)
  
  RERI <- RR11 - RR10 - RR01 + 1
  AP   <- RERI / RR11
  
  grad <- c(
    d_b1  = RR11 - RR10,
    d_b2  = RR11 - RR01,
    d_b12 = RR11
  )
  
  idx <- c(which(names(coefs) == term_pso),
           which(names(coefs) == term_strat),
           which(names(coefs) == term_interact))
  sub_vcv <- vcv[idx, idx]
  
  var_reri <- t(grad) %*% sub_vcv %*% grad
  se_reri  <- sqrt(var_reri)
  
  reri_lower <- RERI - 1.96 * se_reri
  reri_upper <- RERI + 1.96 * se_reri
  
  ap_lower <- AP - 1.96 * (se_reri / RR11)
  ap_upper <- AP + 1.96 * (se_reri / RR11)
  
  data.frame(
    Variable = strat_var_name,
    Subgroup_Comparison = paste0(strat_levels[2], " vs ", strat_levels[1]),
    HR_Pso_Alone = round(RR10, 3),
    HR_Sub_Alone = round(RR01, 3),
    HR_Both      = round(RR11, 3),
    RERI         = round(RERI, 3),
    RERI_95CI    = paste0("[", round(reri_lower, 3), ", ", round(reri_upper, 3), "]"),
    AP           = round(AP, 3),
    AP_95CI      = paste0("[", round(ap_lower, 3), ", ", round(ap_upper, 3), "]"),
    n            = m_add$n,
    events       = m_add$nevent,
    stringsAsFactors = FALSE
  )
}

additive_specs <- list(
  list(name = "Gender (p31)", col = "p31"),
  list(name = "Age Group (<60 vs >=60)", col = "age_group"),
  list(name = "BMI Group (<30 vs >=30)", col = "bmi_group"),
  list(name = "Smoking Status (Binary)", col = "smoking_binary")
)

all_additive_results <- list()
for (spec in additive_specs) {
  res_add <- run_additive_interaction_with_ci(td_data, spec$name, spec$col, valid_covs)
  if (!is.null(res_add)) {
    all_additive_results[[length(all_additive_results) + 1]] <- res_add
  }
}

if (length(all_additive_results) > 0) {
  final_additive_summary <- do.call(rbind, all_additive_results)
  write.csv(final_additive_summary, "pneumonia_additive_interaction_with_CIs.csv", row.names = FALSE)
  cat("\nSUCCESS! Additive interaction results with 95% CIs exported successfully!\n")
} else {
  cat("\nNOTICE: No results generated.\n")
}

==============================================================================
# Script 3 (Expanded): Prospective Cox Proportional Hazards Regression
# Description: Evaluates prospective associations between baseline psoriasis and all respiratory outcomes (composite categories and individual columns 55:74) using external follow-up time data and proper prevalent case exclusion.

==============================================================================
library(survival)

# 1. Load dataset
data2 <- read.csv("total_data_new_407208.csv")

# 2. Define key column mappings
date_baseline_col <- "p53_i0"
date_pso_col      <- "p131742"
admin_cens_date   <- as.Date("2024-03-31", "%Y-%m-%d") # Update to your actual censoring date
cov_cols          <- 25:53

# Define outcome categories and their respective column ranges
outcome_categories <- list(
  "Pneumonia"       = 4:10,
  "Upper_RI"        = 11:13,
  "Lower_RI"        = 14:20,
  "Influenza"       = 21:23
)

# Convert baseline and psoriasis dates
base_date <- as.Date(trimws(as.character(data2[[date_baseline_col]])), "%Y-%m-%d")
pso_date  <- as.Date(trimws(as.character(data2[[date_pso_col]])), "%Y-%m-%d")

# Convert character covariates to factors
for (nm in names(data2)[cov_cols]) {
  if (is.character(data2[[nm]])) data2[[nm]] <- as.factor(data2[[nm]])
}
cov_names <- names(data2)[cov_cols]

# 3. Universal function for time-varying Cox regression per outcome category
run_time_varying_cox <- function(df_full, outcome_name, outcome_cols) {
  
  # Extract earliest diagnosis date for the outcome category across columns
  out_dates <- apply(df_full[, outcome_cols, drop = FALSE], 1, function(row) {
    dates <- as.Date(trimws(as.character(row)), "%Y-%m-%d")
    dates <- dates[!is.na(dates)]
    if (length(dates) > 0) min(dates) else as.Date(NA)
  })
  
  # Construct base working dataframe
  d_sub <- data.frame(
    eid = df_full$eid,
    base_date = base_date,
    pso_date = pso_date,
    out_date = out_dates,
    stringsAsFactors = FALSE
  )
  
  for (nm in cov_names) {
    d_sub[[nm]] <- df_full[[nm]]
  }
  
  # Exclude prevalent cases (outcome occurred at or before baseline)
  has_prevalent <- !is.na(d_sub$out_date) & d_sub$out_date <= d_sub$base_date
  d_sub <- d_sub[!has_prevalent, ]
  
  # Determine follow-up end date (earliest of outcome date, death/admin censoring date)
  d_sub$end_date <- pmin(d_sub$out_date, admin_cens_date, na.rm = TRUE)
  
  # Calculate time differences in days from baseline
  d_sub$t_end <- as.numeric(d_sub$end_date - d_sub$base_date)
  d_sub$t_pso <- as.numeric(d_sub$pso_date - d_sub$base_date)
  d_sub$t_out <- as.numeric(d_sub$out_date - d_sub$base_date)
  
  # Filter valid positive follow-up duration
  d_sub <- d_sub[is.finite(d_sub$t_end) & d_sub$t_end > 0, ]
  
  # Define event indicator (1 if outcome occurred at end_date, 0 otherwise)
  d_sub$event <- ifelse(!is.na(d_sub$t_out) & d_sub$t_out == d_sub$end_date, 1L, 0L)
  
  # Check for minimum event count
  if (sum(d_sub$event == 1, na.rm = TRUE) < 5) return(NULL)
  
  # Filter valid covariates with >= 2 unique levels
  valid_covs <- Filter(function(nm) {
    length(unique(na.omit(d_sub[[nm]]))) >= 2
  }, cov_names)
  cov_terms <- paste(valid_covs, collapse = " + ")
  
  # Use tmerge to create time-varying exposure data structure
  td_data <- tryCatch({
    tmerge(
      data1 = d_sub[, c("eid", "t_end", "event", valid_covs)],
      data2 = d_sub[, c("eid", "t_pso")],
      id = eid,
      death = event(t_end, event),
      psoriasis_tv = tdc(t_pso)
    )
  }, error = function(e) { NULL })
  
  if (is.null(td_data)) return(NULL)
  
  # --- 1. Adjusted Time-Varying Cox Model ---
  fml_adj <- as.formula(paste0("Surv(tstart, tstop, death) ~ psoriasis_tv + ", cov_terms))
  res_adj <- tryCatch({
    m_adj <- coxph(fml_adj, data = td_data)
    s_adj <- summary(m_adj)
    if (!("psoriasis_tv" %in% rownames(s_adj$coef))) stop("Exposure dropped")
    
    beta <- s_adj$coef["psoriasis_tv", "coef"]
    se   <- s_adj$coef["psoriasis_tv", "se(coef)"]
    pval <- s_adj$coef["psoriasis_tv", "Pr(>|z|)"]
    
    ph_p_exp <- NA; ph_p_glob <- NA
    tryCatch({
      zph <- cox.zph(m_adj)
      if ("psoriasis_tv" %in% rownames(zph$table)) ph_p_exp <- zph$table["psoriasis_tv", "p"]
      if ("GLOBAL" %in% rownames(zph$table)) ph_p_glob <- zph$table["GLOBAL", "p"]
    }, error = function(e) {})
    
    data.frame(
      Outcome = outcome_name,
      Model = "Adjusted",
      HR = exp(beta),
      CI_low = exp(beta - 1.96 * se),
      CI_high = exp(beta + 1.96 * se),
      p = pval,
      PH_p_exposure = ph_p_exp,
      PH_p_global = ph_p_glob,
      n = m_adj$n,
      events = m_adj$nevent,
      stringsAsFactors = FALSE
    )
  }, error = function(e) { NULL })
  
  # --- 2. Unadjusted Time-Varying Cox Model ---
  fml_unadj <- Surv(tstart, tstop, death) ~ psoriasis_tv
  res_unadj <- tryCatch({
    m_unadj <- coxph(fml_unadj, data = td_data)
    s_unadj <- summary(m_unadj)
    
    beta <- s_unadj$coef["psoriasis_tv", "coef"]
    se   <- s_unadj$coef["psoriasis_tv", "se(coef)"]
    pval <- s_unadj$coef["psoriasis_tv", "Pr(>|z|)"]
    
    ph_p_exp0 <- NA; ph_p_glob0 <- NA
    tryCatch({
      zph0 <- cox.zph(m_unadj)
      if ("psoriasis_tv" %in% rownames(zph0$table)) ph_p_exp0 <- zph0$table["psoriasis_tv", "p"]
      if ("GLOBAL" %in% rownames(zph0$table)) ph_p_glob0 <- zph0$table["GLOBAL", "p"]
    }, error = function(e) {})
    
    data.frame(
      Outcome = outcome_name,
      Model = "Unadjusted",
      HR = exp(beta),
      CI_low = exp(beta - 1.96 * se),
      CI_high = exp(beta + 1.96 * se),
      p = pval,
      PH_p_exposure = ph_p_exp0,
      PH_p_global = ph_p_glob0,
      n = m_unadj$n,
      events = m_unadj$nevent,
      stringsAsFactors = FALSE
    )
  }, error = function(e) { NULL })
  
  rbind(res_unadj, res_adj)
}

# 4. Batch loop across all outcome categories
results_list <- list()
for (cat_name in names(outcome_categories)) {
  cols <- outcome_categories[[cat_name]]
  res <- run_time_varying_cox(data2, cat_name, cols)
  if (!is.null(res)) {
    results_list[[length(results_list) + 1]] <- res
  }
}

# 5. Combine and export results
final_tv_summary <- do.call(rbind, results_list)
write.csv(final_tv_summary, "cox_time_varying_respiratory_summary.csv", row.names = FALSE)
# 7. Combine and export final results
final_cox_summary <- do.call(rbind, results_list)
write.csv(final_cox_summary, "cox_all_respiratory_outcomes_summary.csv", row.names = FALSE)

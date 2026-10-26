# ==============================================================================
# Script 2 (Expanded): Cross-Sectional Analysis for All Individual & Composite Outcomes
# Description: Runs cross-sectional logistic regressions for all individual and composite outcomes under two datasets:
#              1. Main Analysis (Keeping incident followup cases as unexposed)
#              2. Sensitivity Analysis (Excluding incident followup cases)
# ==============================================================================
# 1. Load dataset
data2 <- read.csv("total_data_new_407208.csv")

# 2. Define basic variables and indices
cov_cols  <- 25:53
x_col     <- 54
x_name    <- names(data2)[x_col]

# 3. Binarize exposure function
make_bin <- function(v) {
  if (is.factor(v)) v <- as.character(v)
  if (is.character(v)) {
    out <- ifelse(v == "baseline", 1L, 0L)
    out[is.na(v)] <- 0L
    return(out)
  } else {
    out <- as.integer(v)
    out[is.na(out)] <- 0L
    return(out)
  }
}

# 4. Universal outcome construction and model fitting function
run_all_outcomes_for_df <- function(df_input, dataset_label) {
  
  df <- df_input
  df[[x_name]] <- make_bin(df[[x_col]])
  
  # Convert covariates to factors
  cov_names <- names(df)[cov_cols]
  for (nm in cov_names) {
    if (is.character(df[[nm]])) df[[nm]] <- as.factor(df[[nm]])
  }
  
  # Define outcome sets
  composite_outcomes <- list(
    "Pneumonia_Composite"      = 55:61,
    "Lower_RI_Composite"       = 62:64,
    "Upper_RI_Composite"       = 65:71,
    "Influenza_Composite"      = 72:74
  )
  
  for (comp_name in names(composite_outcomes)) {
    cols <- composite_outcomes[[comp_name]]
    df[[comp_name]] <- apply(df[, cols], 1, function(row) {
      row_char <- as.character(row)
      if (any(row_char == "baseline", na.rm = TRUE)) 1L else 0L
    })
  }
  
  individual_outcomes <- 55:74
  for (col_idx in individual_outcomes) {
    col_orig_name <- names(df)[col_idx]
    new_col_name <- paste0("indiv_", col_idx, "_", col_orig_name)
    df[[new_col_name]] <- make_bin(df[[col_idx]])
  }
  
  all_target_outcomes <- c(
    names(composite_outcomes),
    paste0("indiv_", 55:74, "_", names(df)[55:74])
  )
  
  # Regression function for individual outcomes
  run_logistic_models <- function(data_subset, outcome_col_name, cov_list) {
    y_vals <- data_subset[[outcome_col_name]]
    if (length(unique(na.omit(y_vals))) < 2 || sum(y_vals == 1, na.rm = TRUE) < 5) {
      return(NULL)
    }
    
    valid_covs <- Filter(function(nm) {
      length(unique(na.omit(data_subset[[nm]]))) >= 2
    }, cov_list)
    
    cov_terms <- paste(valid_covs, collapse = " + ")
    
    # Adjusted Model
    f_adj <- as.formula(paste0(outcome_col_name, " ~ ", x_name, " + ", cov_terms))
    res_adj_df <- tryCatch({
      m_adj <- glm(f_adj, data = data_subset, family = binomial())
      cf_adj <- summary(m_adj)$coefficients
      if (!(x_name %in% rownames(cf_adj))) stop("Exposure dropped")
      
      beta <- cf_adj[x_name, "Estimate"]
      se   <- cf_adj[x_name, "Std. Error"]
      pval <- cf_adj[x_name, "Pr(>|z|)"]
      
      data.frame(
        Analysis_Set = dataset_label,
        Outcome = outcome_col_name,
        Model = "Adjusted",
        OR = exp(beta),
        CI_low = exp(beta - 1.96 * se),
        CI_high = exp(beta + 1.96 * se),
        p = pval,
        n = nobs(m_adj),
        events = sum(y_vals == 1, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }, error = function(e) { NULL })
    
    # Unadjusted Model
    f_unadj <- as.formula(paste0(outcome_col_name, " ~ ", x_name))
    res_unadj_df <- tryCatch({
      m_unadj <- glm(f_unadj, data = data_subset, family = binomial())
      cf_unadj <- summary(m_unadj)$coefficients
      
      beta <- cf_unadj[x_name, "Estimate"]
      se   <- cf_unadj[x_name, "Std. Result"]
      pval <- cf_unadj[x_name, "Pr(>|z|)"]
      
      data.frame(
        Analysis_Set = dataset_label,
        Outcome = outcome_col_name,
        Model = "Unadjusted",
        OR = exp(beta),
        CI_low = exp(beta - 1.96 * se),
        CI_high = exp(beta + 1.96 * se),
        p = pval,
        n = nobs(m_unadj),
        events = sum(y_vals == 1, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    }, error = function(e) { NULL })
    
    rbind(res_unadj_df, res_adj_df)
  }
  
  # Loop through all outcomes
  results_list <- list()
  for (out_name in all_target_outcomes) {
    res_res <- run_logistic_models(df, out_name, cov_names)
    if (!is.null(res_res)) {
      results_list[[length(results_list) + 1]] <- res_res
    }
  }
  
  do.call(rbind, results_list)
}

# 5. Execute Analysis 1: Main Analysis (Keep all / Do not exclude followup)
df_main <- data2
res_main <- run_all_outcomes_for_df(df_main, dataset_label = "Main_KeepFollowup")

# 6. Execute Analysis 2: Sensitivity Analysis (Exclude followup cases)
df_sens <- data2
exclude_idx <- !is.na(df_sens[[x_col]]) & trimws(as.character(df_sens[[x_col]])) == "followup"
df_sens <- df_sens[!exclude_idx, ]
res_sens <- run_all_outcomes_for_df(df_sens, dataset_label = "Sensitivity_ExcludeFollowup")

# 7. Combine both sets of results and export to a single CSV
final_combined_table <- rbind(res_main, res_sens)
output_filename <- "logistic_main_vs_sensitivity_summary.csv"
write.csv(final_combined_table, output_filename, row.names = FALSE)

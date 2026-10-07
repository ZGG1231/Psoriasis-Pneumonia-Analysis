==============================================================================
# Script 3: Prospective Longitudinal Analysis (Cox Proportional Hazards)
# Description: Derives an at-risk prospective cohort by excluding baseline prevalent cases, calculates incident events, and estimates Hazard Ratios (HR) with Proportional Hazards (PH) validation.
==============================================================================
library(survival)

data2 <- read.csv("total_data_new_407208.csv")

# Step 1: Integrate updated follow-up time data
followup_time <- read.csv("20250613data_followuptime(50w).csv")
fu_colname <- names(followup_time)[124]

data2 <- merge(
  data2,
  followup_time[, c("eid", fu_colname)],
  by = "eid",
  all.x = TRUE,
  sort = FALSE
)
names(data2)[names(data2) == fu_colname] <- "fu_time"

df_raw <- data2 
pneu_names <- names(df_raw)[55:61]
x_name <- names(df_raw)[54]

# Step 2: Exclude prevalent cases to establish the prospective 'at-risk' cohort
is_baseline_case <- apply(df_raw[, pneu_names], 1, function(row) {
  any(as.character(row) == "baseline", na.rm = TRUE)
})
df_cox <- df_raw[!is_baseline_case, ]

# Step 3: Define incident pneumonia events during the follow-up period
df_cox$incident_pneumonia <- apply(df_cox[, pneu_names], 1, function(row) {
  if (any(as.character(row) == "followup", na.rm = TRUE)) 1L else 0L
})

# Format exposure, follow-up time, and covariates
df_cox[[x_name]] <- ifelse(as.character(df_cox[[x_name]]) == "baseline", 1L, 0L)
df_cox[[x_name]][is.na(df_cox[[x_name]])] <- 0L

if (!is.numeric(df_cox$fu_time)) {
  df_cox$fu_time <- as.numeric(as.character(df_cox$fu_time))
}

cov_names <- names(df_cox)[25:53]
for (nm in cov_names) {
  if (is.character(df_cox[[nm]])) df_cox[[nm]] <- as.factor(df_cox[[nm]])
}

# Define clinical subgroup variables
df_cox$age_group <- ifelse(df_cox$age >= 60, ">=60", "<60")
df_cox$bmi_group <- ifelse(df_cox$p21001_i0 >= 30, "Obese_geq30", "NonObese_lt30")

# Step 4: Universal Cox Regression function with error handling and PH testing
run_cox_model <- function(data_subset, strata_label, exclude_cov = NULL) {
  
  # Filter out singular covariates for the current subset
  valid_covs <- c()
  for (nm in cov_names) {
    if (!is.null(exclude_cov) && nm %in% exclude_cov) next
    if (length(unique(na.omit(data_subset[[nm]]))) > 1) valid_covs <- c(valid_covs, nm)
  }
  
  current_cov_terms <- paste(valid_covs, collapse = " + ")
  f <- as.formula(paste0("Surv(fu_time, incident_pneumonia) ~ ", x_name, " + ", current_cov_terms))
  
  # Ensure sufficient event counts for stable model convergence
  if (sum(data_subset$incident_pneumonia == 1, na.rm = TRUE) < 5) {
    return(data.frame(Analysis = strata_label, HR = NA, CI_low = NA, CI_high = NA, p = NA, 
                      PH_p_exposure = NA, PH_p_global = NA, n = NA, events = NA, stringsAsFactors = FALSE))
  } 
  
  tryCatch({
    m <- coxph(f, data = data_subset)
    coefs <- summary(m)$coefficients
    
    # Check for collinearity-induced dropouts
    if (!(x_name %in% rownames(coefs))) stop("Exposure dropped")
    
    beta  <- coefs[x_name, "coef"]
    se    <- coefs[x_name, "se(coef)"]
    
    # Assess Proportional Hazards (PH) assumption using Schoenfeld residuals
    ph_p_exp <- NA; ph_p_glob <- NA
    tryCatch({
      zph <- cox.zph(m)
      if (x_name %in% rownames(zph$table)) ph_p_exp <- zph$table[x_name, "p"]
      if ("GLOBAL" %in% rownames(zph$table)) ph_p_glob <- zph$table["GLOBAL", "p"]
    }, error = function(e) { })
    
    data.frame(
      Analysis = strata_label,
      HR = exp(beta),
      CI_low = exp(beta - 1.96 * se),
      CI_high = exp(beta + 1.96 * se),
      p = coefs[x_name, "Pr(>|z|)"],
      PH_p_exposure = ph_p_exp,
      PH_p_global = ph_p_glob,
      n = m$n,
      events = m$nevent,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    message("x ", strata_label, " Failed: ", e$message)
    data.frame(Analysis = strata_label, HR = NA, CI_low = NA, CI_high = NA, p = NA, 
               PH_p_exposure = NA, PH_p_global = NA, n = NA, events = NA, stringsAsFactors = FALSE)
  })
}

# Step 5: Execute main adjusted model and stratifications
res_list_cox <- list()
res_list_cox[[1]] <- run_cox_model(df_cox, "Main_All", exclude_cov = NULL)

# Iterative stratification
for (val in na.omit(unique(df_cox$p31))) res_list_cox[[length(res_list_cox) + 1]] <- run_cox_model(df_cox[df_cox$p31 == val, ], paste0("Gender_", val), "p31")
for (val in na.omit(unique(df_cox$age_group))) res_list_cox[[length(res_list_cox) + 1]] <- run_cox_model(df_cox[df_cox$age_group == val, ], paste0("Age_", val), "age")
for (val in na.omit(unique(df_cox$bmi_group))) res_list_cox[[length(res_list_cox) + 1]] <- run_cox_model(df_cox[df_cox$bmi_group == val, ], paste0("BMI_", val), "p21001_i0")
for (val in na.omit(unique(df_cox$p20116_i0))) res_list_cox[[length(res_list_cox) + 1]] <- run_cox_model(df_cox[df_cox$p20116_i0 == val, ], paste0("Smoking_", val), "p20116_i0")

# Compile and export results
final_res_cox <- do.call(rbind, res_list_cox)
write.csv(final_res_cox, "Cox_Incident_Pneumonia_Stratified.csv", row.names = FALSE)

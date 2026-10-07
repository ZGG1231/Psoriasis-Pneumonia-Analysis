==============================================================================
# Script 2: Cross-Sectional Analysis (Logistic Regression)
# Description: Evaluates cross-sectional associations between psoriasis and composite pneumonia using both unadjusted and adjusted models, including predefined subgroup stratifications.
==============================================================================
data2 <- read.csv("total_data_new_407208.csv")

# Define variable indices
cov_cols  <- 25:53
x_col     <- 54
pneu_cols <- 55:61  

df <- data2

# Step 1: Binarize baseline exposures and outcomes (Baseline = 1, Other = 0)
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

x_name <- names(df)[x_col]
df[[x_name]] <- make_bin(df[[x_col]])

# Step 2: Construct composite outcome (Any Pneumonia)
pneu_names <- names(df)[pneu_cols]
df$pneumonia_any <- apply(df[, pneu_names], 1, function(row) {
  row <- as.character(row)
  if (any(row == "baseline", na.rm = TRUE)) 1L else 0L
})
y_name <- "pneumonia_any"

# Convert character covariates to factors
cov_names <- names(df)[cov_cols]
for (nm in cov_names) {
  if (is.character(df[[nm]])) df[[nm]] <- as.factor(df[[nm]])
}

# Step 3: Unadjusted Logistic Regression Model
f0 <- as.formula(paste0(y_name, " ~ ", x_name))
m0 <- glm(f0, data = df, family = binomial())
coefs0 <- summary(m0)$coefficients

res0 <- data.frame(
  outcome = "Any_pneumonia (Unadjusted)",
  OR = exp(coefs0[x_name, "Estimate"]),
  CI_low = exp(coefs0[x_name, "Estimate"] - 1.96 * coefs0[x_name, "Std. Error"]),
  CI_high = exp(coefs0[x_name, "Estimate"] + 1.96 * coefs0[x_name, "Std. Error"]),
  p = coefs0[x_name, "Pr(>|z|)"],
  n = nobs(m0),
  events = sum(df[[y_name]] == 1),
  stringsAsFactors = FALSE
)

# Step 4: Adjusted Models and Stratification Setup
# Define clinical subgroup variables
df$age_group <- ifelse(df$age >= 60, ">=60", "<60")
df$bmi_group <- ifelse(df$p21001_i0 >= 30, "Obese_geq30", "NonObese_lt30")

# Universal function for adjusted logistic regression with dynamic covariate filtering
run_pneu_model <- function(data_subset, strata_label, exclude_cov = NULL) {
  
  # Remove covariates with zero variance within the specific stratum
  valid_covs <- c()
  for (nm in cov_names) {
    if (!is.null(exclude_cov) && nm %in% exclude_cov) next
    if (length(unique(na.omit(data_subset[[nm]]))) > 1) {
      valid_covs <- c(valid_covs, nm)
    }
  }
  
  current_cov_terms <- paste(valid_covs, collapse = " + ")
  f <- as.formula(paste0("pneumonia_any ~ ", x_name, " + ", current_cov_terms))
  
  # Skip model fitting if the outcome lacks variation in the current subset
  if (length(unique(na.omit(data_subset$pneumonia_any))) < 2) {
    return(data.frame(Analysis = strata_label, OR = NA, CI_low = NA, CI_high = NA, p = NA, n = NA, events = NA, stringsAsFactors = FALSE))
  }
  
  tryCatch({
    m <- glm(f, data = data_subset, family = binomial())
    coefs <- summary(m)$coefficients
    if (!(x_name %in% rownames(coefs))) stop("Exposure dropped")
    
    beta  <- coefs[x_name, "Estimate"]
    se    <- coefs[x_name, "Std. Error"]
    
    data.frame(
      Analysis = strata_label,
      OR = exp(beta),
      CI_low = exp(beta - 1.96 * se),
      CI_high = exp(beta + 1.96 * se),
      p = coefs[x_name, "Pr(>|z|)"],
      n = nobs(m),
      events = sum(data_subset$pneumonia_any == 1, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    data.frame(Analysis = strata_label, OR = NA, CI_low = NA, CI_high = NA, p = NA, n = NA, events = NA, stringsAsFactors = FALSE)
  })
}

# Step 5: Execute main adjusted model and stratifications
res_list <- list()
res_list[[1]] <- run_pneu_model(df, "Main_All", exclude_cov = NULL)

# Iterative stratification by Gender, Age, BMI, and Smoking Status
for (val in na.omit(unique(df$p31))) res_list[[length(res_list) + 1]] <- run_pneu_model(df[df$p31 == val, ], paste0("Gender_", val), "p31")
for (val in na.omit(unique(df$age_group))) res_list[[length(res_list) + 1]] <- run_pneu_model(df[df$age_group == val, ], paste0("Age_", val), "age")
for (val in na.omit(unique(df$bmi_group))) res_list[[length(res_list) + 1]] <- run_pneu_model(df[df$bmi_group == val, ], paste0("BMI_", val), "p21001_i0")
for (val in na.omit(unique(df$p20116_i0))) res_list[[length(res_list) + 1]] <- run_pneu_model(df[df$p20116_i0 == val, ], paste0("Smoking_", val), "p20116_i0")

# Compile and export results
final_res <- do.call(rbind, res_list)
write.csv(final_res, "logistic_pneumonia_any_Main_and_Stratified.csv", row.names = FALSE)

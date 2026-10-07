==============================================================================
# Script 1: Data Cleaning and Preprocessing
# Description: Merges baseline phenotype data with follow-up datasets, handles 
#              missing values, and validates date formats for longitudinal analysis.
==============================================================================
library(dplyr)

# Step 1: Load raw data and follow-up time dataset
data <- read.csv("psoriasis_respiratory_infections.csv")
followup_time <- read.csv("20250612data_followuptime.csv")

# Remove duplicated participant IDs in the follow-up dataset
followup_time2 <- followup_time[!duplicated(followup_time$eid), c("eid", "age")]

# Match and merge age data to the primary dataset
idx <- match(data$eid, followup_time2$eid)
age_vec <- followup_time2$age[idx]

pos <- 24
data2 <- cbind(
  data[, 1:pos, drop = FALSE],
  age = age_vec,
  if (ncol(data) > pos) data[, (pos + 1):ncol(data), drop = FALSE] else NULL
)

# Step 2: Define column indices for baseline, diseases, and covariates
baseline_col   <- 2
disease_cols   <- 3:24
filter_cols    <- 25:53

# Identify and remove rows with generalized missing values in covariates
missing_matrix <- sapply(data2[, filter_cols], function(x) {
  x_char <- trimws(as.character(x)) 
  is.na(x) | x_char == "NA" | x_char == "N/A" | x_char == "null" | x_char == ""
})
na_counts <- colSums(missing_matrix, na.rm = TRUE)

rm_row <- apply(missing_matrix, 1, any, na.rm = TRUE)
data2 <- data2[!rm_row, ]

# Step 3: Standardize date columns to Date objects
base_date <- as.Date(as.character(data2[[baseline_col]]), format = "%Y-%m-%d")
for (j in disease_cols) {
  data2[[j]] <- as.Date(as.character(data2[[j]]), format = "%Y-%m-%d")
}

# Step 4: Categorize disease status based on baseline date
disease_names <- names(data2)[disease_cols]
group_cols <- paste0(disease_names, "_group")

for (k in seq_along(disease_cols)) {
  j  <- disease_cols[k]
  nm <- group_cols[k]
  d  <- data2[[j]]
  data2[[nm]] <- ifelse(
    is.na(d) | is.na(base_date), NA_character_,
    ifelse(d <= base_date, "baseline", "followup")
  )
}

# Generate summary counts for prevalent (baseline) and incident (followup) cases
counts <- do.call(rbind, lapply(group_cols, function(g) {
  v <- data2[[g]]
  c(
    baseline = sum(v == "baseline", na.rm = TRUE),
    followup = sum(v == "followup", na.rm = TRUE),
    missing  = sum(is.na(v))
  )
}))

# Step 5: Strict validation of date formats using regular expressions
cols <- 3:24
x_chr <- as.data.frame(lapply(data2[, cols], function(col) trimws(as.character(col))))
is_blank <- is.na(x_chr) | x_chr == ""
is_date <- as.data.frame(lapply(x_chr, function(v) grepl("^\\d{4}-\\d{2}-\\d{2}$", v)))

# Remove rows containing non-empty, non-date string artifacts
bad_cell <- (!is_blank) & (!as.matrix(is_date))
rm_row <- apply(bad_cell, 1, any)

data2 <- data2[!rm_row, ]

# Export cleaned dataset
write.csv(data2, "total_data_new_407208.csv", row.names = FALSE)

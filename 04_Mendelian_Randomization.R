==============================================================================
# Script 4: Mendelian Randomization and Steiger Test (04_Mendelian_Randomization.R)
# Description: Instrumental Variable (IV) selection, EAF alignment, local clumping, TwoSampleMR analysis, sensitivity tests, and Steiger directionality.
==============================================================================
library(data.table)
library(TwoSampleMR)
library(writexl)
# if(!require("plinkbinr")) remotes::install_github("explodecomputer/plinkbinr")
library(plinkbinr)

# Part 1: IV Selection & EAF Alignment
# Load and format outcome (pneumonia) GWAS summary statistics
gwas_pneu <- fread("pneumoniaGWAS.h.tsv")
setnames(gwas_pneu, 
         old = c("SNP", "CHR", "BP", "A1", "A2", "BETA", "sebeta", "P", "FRQ"), 
         new = c("SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P", "EAF"), 
         skip_absent = TRUE)

mr_cols <- c("SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P", "EAF")
pneu_mr_ready <- gwas_pneu[, ..mr_cols]
pneu_mr_ready[, BP := as.numeric(BP)]

# Exclude the Major Histocompatibility Complex (MHC) region due to complex LD patterns
pneu_mr_ready <- pneu_mr_ready[!(pneu_mr_ready$CHR %in% c(6, "6", "chr6") & pneu_mr_ready$BP >= 28510120 & pneu_mr_ready$BP <= 33480577), ]
write.table(pneu_mr_ready, file = "pneu_mr_ready.txt", sep = "\t", row.names = FALSE, quote = FALSE) 

# Load and format exposure (psoriasis) GWAS summary statistics
gwas_pso <- fread("psoriasisGWAS.h.tsv")
setnames(gwas_pso, 
         old = c("rsid", "chromosome", "base_pair_location", "effect_allele", "other_allele", "beta", "standard_error", "p_value"), 
         new = c("SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P"), 
         skip_absent = TRUE)

gwas_pso[, BP := as.numeric(BP)]  
gwas_pso[, P := as.numeric(P)]
gwas_pso[, BETA := as.numeric(BETA)]
gwas_pso[, SE := as.numeric(SE)]

# Exclude MHC region from exposure data
gwas_pso <- gwas_pso[!(CHR %in% c(6, "6", "chr6") & BP >= 28510120 & BP <= 33480577)]

# Select genome-wide significant SNPs and filter out weak instruments (F-statistic <= 10)
pso_ivs <- gwas_pso[P < 5e-8]
pso_ivs[, F_stat := (BETA / SE)^2]
pso_ivs_strong <- pso_ivs[F_stat > 10]

pso_cols <- c("SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P", "F_stat")
pso_mr_ready <- pso_ivs_strong[, ..pso_cols]

# Extract Effect Allele Frequencies (EAF) from outcome data as reference
setDT(pneu_mr_ready)
setDT(pso_mr_ready)
ref_freq <- pneu_mr_ready[, .(SNP, ref_A1 = A1, ref_A2 = A2, ref_EAF = EAF)]

pso_merged <- merge(pso_mr_ready, ref_freq, by = "SNP", all.x = TRUE)

# Standardize allele capitalization to prevent mismatch
cols_to_upper <- c("A1", "A2", "ref_A1", "ref_A2")
pso_merged[, (cols_to_upper) := lapply(.SD, toupper), .SDcols = cols_to_upper]
pso_merged[, EAF := NA_real_]

# Align EAF based on matching alleles
pso_merged[A1 == ref_A1, EAF := ref_EAF]
pso_merged[A1 == ref_A2, EAF := 1 - ref_EAF]

final_cols <- c("SNP", "CHR", "BP", "A1", "A2", "BETA", "SE", "P", "EAF", "F_stat")
pso_exposure_final <- pso_merged[!is.na(EAF), ..final_cols]
write.table(pso_exposure_final, file = "pso_mr_ready.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Part 2: Formatting & Local Clumping
pso_exposure_final <- as.data.frame(pso_exposure_final)
pneu_mr_ready <- as.data.frame(pneu_mr_ready)

# Format for TwoSampleMR package
exp_dat <- format_data(
   pso_exposure_final, type = "exposure", snp_col = "SNP", beta_col = "BETA",
   se_col = "SE", effect_allele_col = "A1", other_allele_col = "A2",
   pval_col = "P", eaf_col = "EAF"
)

# Set up local 1000 Genomes reference panel for LD clumping
options(timeout = 3600)
ref_dir <- "1kg_reference_local"
tgz_file <- file.path(ref_dir, "1kg.v3.tgz")
if(file.exists(tgz_file)) {
   file.remove(tgz_file)
}
# download.file("http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz", destfile = tgz_file, mode = "wb")
# untar(tgz_file, exdir = ref_dir)

# Perform local LD clumping (r2 < 0.001, 10000 kb window)
bfile_path <- file.path(ref_dir, "EUR")
exp_dat_clumped <- clump_data(
   exp_dat, clump_kb = 10000, clump_r2 = 0.001, bfile = bfile_path,
   plink_bin = plinkbinr::get_plink_exe()
)

# Part 3: Outcome Extraction & Harmonization
# Extract clumped IVs from the outcome dataset
out_dat <- format_data(
   pneu_mr_ready, type = "outcome", snps = exp_dat_clumped$SNP, snp_col = "SNP",
   beta_col = "BETA", se_col = "SE", effect_allele_col = "A1",
   other_allele_col = "A2", pval_col = "P", eaf_col = "EAF"
)

# Harmonize data (action = 2: assume forward strand and align alleles)
dat <- harmonise_data(exposure_dat = exp_dat_clumped, outcome_dat = out_dat, action = 2)

# Part 4: Main MR & Sensitivity Analyses
dat_clean <- dat 

# Main Mendelian Randomization analyses
res <- mr(dat_clean)
res_format <- generate_odds_ratios(res)
print(res_format[, c("method", "nsnp", "or", "or_lci95", "or_uci95", "pval")])

# Cochran's Q test for heterogeneity
het <- mr_heterogeneity(dat_clean)
print(het[, c("method", "Q", "Q_df", "Q_pval")])

# MR-Egger intercept test for directional horizontal pleiotropy
pleio <- mr_pleiotropy_test(dat_clean)
print(pleio[, c("egger_intercept", "se", "pval")])

# Leave-one-out sensitivity analysis
loo_clean <- mr_leaveoneout(dat_clean)
p_loo <- mr_leaveoneout_plot(loo_clean)
print(p_loo)

# Part 5: Single SNP Analysis & Steiger Directionality Test
# Calculate Wald ratio for each single SNP
res_single <- mr_singlesnp(dat_clean)
res_single_snps <- res_single[!grepl("All", res_single$SNP), ]
res_single_snps$FDR <- p.adjust(res_single_snps$p, method = "fdr")
res_single_snps <- res_single_snps[order(res_single_snps$FDR), ]
write.table(res_single_snps, file = "114_IVs_mr_wald.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Define parameters for Steiger directionality test
ncase_pso <- 36466; ncontrol_pso <- 458078; prev_pso <- 0.02
ncase_pneu <- 78791; ncontrol_pneu <- 421557; prev_pneu <- 0.05

dat_clean$r.exposure <- get_r_from_lor(lor = dat_clean$beta.exposure, af = dat_clean$eaf.exposure, ncase = ncase_pso, ncontrol = ncontrol_pso, prevalence = prev_pso)
dat_clean$r.outcome <- get_r_from_lor(lor = dat_clean$beta.outcome, af = dat_clean$eaf.outcome, ncase = ncase_pneu, ncontrol = ncontrol_pneu, prevalence = prev_pneu)

# Perform Steiger test to verify causal direction (Exposure -> Outcome)
steiger_res <- directionality_test(dat_clean)
steiger_res$steiger_pval_sci <- sprintf("%.2e", steiger_res$steiger_pval)
print(steiger_res[, c("snp_r2.exposure", "snp_r2.outcome", "correct_causal_direction", "steiger_pval_sci")])

# Merge Steiger results with Single SNP Wald ratios
cols_to_extract <- c("SNP", "r.exposure", "r.outcome", "is_direction_correct")
steiger_info <- dat_clean[, cols_to_extract]
res_single_snps <- merge(res_single_snps, steiger_info, by = "SNP", all.x = TRUE)
res_single_snps <- res_single_snps[order(res_single_snps$p), ]
write_xlsx(res_single_snps, "114IVs_MR&steiger_all_results.xlsx")

# Filter for significant SNPs with correct causal direction
dat_clean_steiger <- dat_clean[dat_clean$is_direction_correct == TRUE, ]
sig_snps <- res_single_snps[res_single_snps$p < 0.05, "SNP"]
final_target_snps <- intersect(sig_snps, dat_clean_steiger$SNP)

final_targets_info <- res_single_snps[res_single_snps$SNP %in% final_target_snps, ]
final_targets_info <- final_targets_info[order(final_targets_info$p), ]
write.table(final_targets_info, file = "11_IVs_wald_P_0.05_and_steiger.txt", sep = "\t", row.names = FALSE, quote = FALSE)

==============================================================================
# Script 5: Bayesian Colocalization Analysis (05_Colocalization.R)
# Description: Performs coloc.abf to identify shared causal variants between exposure and outcome traits within a +/- 500kb window of lead SNPs.
==============================================================================
library(coloc)
library(data.table)

# Ensure full summary statistics are in data.table format for fast subsetting
setDT(gwas_pso)
setDT(pneu_mr_ready)

all_coloc_results <- data.frame()

# Iterate colocalization analysis over significant SNPs from MR
for (i in 1:nrow(res_single_snps)) {
   my_target <- res_single_snps$SNP[i]
   
   # Locate the target SNP in the outcome dataset
   snp_info <- pneu_mr_ready[pneu_mr_ready$SNP == my_target, ]
   if(nrow(snp_info) == 0) {
      next
   }
   
   # Define a 1Mb window (+/- 500kb) centered on the lead SNP
   chr_val <- snp_info$CHR[1]
   bp_val <- snp_info$BP[1]
   w_start <- bp_val - 500000
   w_end <- bp_val + 500000

   # Extract regional summary statistics for both traits
   reg_pso <- gwas_pso[CHR == chr_val & BP >= w_start & BP <= w_end]
   reg_pneu <- pneu_mr_ready[CHR == chr_val & BP >= w_start & BP <= w_end]
   
   # Sort by P-value and remove duplicate SNPs to retain the most significant associations
   reg_pso <- reg_pso[order(P)]
   reg_pso <- reg_pso[!duplicated(reg_pso$SNP)]
   
   reg_pneu <- reg_pneu[order(P)]
   reg_pneu <- reg_pneu[!duplicated(reg_pneu$SNP)]
   
   # Harmonize SNPs between the two datasets
   shared_snps <- intersect(reg_pso$SNP, reg_pneu$SNP)
   
   # Skip regions with insufficient overlapping SNPs to ensure robust Bayesian inference
   if(length(shared_snps) < 10) {
            next
   }
   
   reg_pso <- reg_pso[SNP %in% shared_snps]
   reg_pneu <- reg_pneu[SNP %in% shared_snps]
   setorder(reg_pso, SNP)
   setorder(reg_pneu, SNP)

   # Prepare dataset list for exposure (Psoriasis)
   dataset_pso <- list(
      pvalues = reg_pso$P,
      N = 494544,                                      
      s = 36466 / 494544,                              
      type = "cc",                                     
      MAF = reg_pneu$EAF,                              
      snp = reg_pso$SNP
   )
   
   # Prepare dataset list for outcome (Pneumonia)
   dataset_pneu <- list(
      pvalues = reg_pneu$P,
      N = 500348,                                      
      s = 78791 / 500348,                              
      type = "cc",
      MAF = reg_pneu$EAF,
      snp = reg_pneu$SNP
   )
   
   # Execute Bayesian colocalization with error handling for convergence issues
   coloc_run <- tryCatch({
      coloc.abf(dataset1 = dataset_pso, dataset2 = dataset_pneu)
   }, error = function(e) {
      cat("Failed:", conditionMessage(e), "\n")
      return(NULL)
   })
   
   if (is.null(coloc_run)) next 
   
   # Extract and aggregate posterior probabilities (PP.H0 to PP.H4)
   res_row <- as.data.frame(t(coloc_run$summary))
   
   res_row$Target_SNP <- my_target
   res_row$CHR <- chr_val
   res_row$BP <- bp_val
   
   all_coloc_results <- rbind(all_coloc_results, res_row)
}

# Format final results table and sort by PP.H4 (probability of a shared causal variant)
all_coloc_results <- all_coloc_results[, c("Target_SNP", "CHR", "BP", "nsnps", "PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf")]
all_coloc_results <- all_coloc_results[order(-all_coloc_results$PP.H4.abf), ]

# Export colocalization results
write.csv(all_coloc_results, file = "Coloc_Batch_Results_114IVs.csv", row.names = FALSE)

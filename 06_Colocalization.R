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

###Forest Plot of Mendelian Randomization & Coloc Results
data <- read.delim(file.choose(), header = TRUE, sep = "\t")
library(ggplot2)
library(dplyr)
plot_data <- data %>%
   mutate(
      # Calculate OR and 95% confidence interval
      OR = exp(b),
      LCI = exp(b - 1.96 * se),
      UCI = exp(b + 1.96 * se),
      # Format P values and FDR for display on the plot (optional)
      p_label = sprintf("%.4f", p),
      FDR_label = sprintf("%.4f", FDR)
   )
# To make the forest plot more aesthetically pleasing, reorder SNPs by OR value
plot_data$SNP <- factor(plot_data$SNP, levels = plot_data$SNP[order(plot_data$OR)])

p <- ggplot(plot_data, aes(y = SNP, x = OR)) +
   # Draw confidence interval error bars
   geom_errorbarh(aes(xmin = LCI, xmax = UCI), height = 0.2, color = "#2c3e50", size = 0.8) +
   # Draw effect size points (squares)
   geom_point(shape = 15, size = 3.5, color = "#e74c3c") +
   # Add null line (vertical line at OR = 1)
   geom_vline(xintercept = 1, linetype = "dashed", color = "darkgray", size = 1) +
   # Use a clean academic theme
   theme_classic() +
   # Adjust fonts and axis labels
   theme(
      axis.text.y = element_text(size = 11, face = "bold", color = "black"),
      axis.text.x = element_text(size = 11, color = "black"),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
   ) +
   # Set title and axis labels
   labs(
      title = "Forest Plot of Individual SNPs",
      x = "Odds Ratio (95% CI)",
      y = "SNP"
   )
print(p)

### MR scatter plot
library(TwoSampleMR)
mr_results <- mr(dat_clean_steiger)
p1 <- mr_scatter_plot(mr_results, dat_clean_steiger)
p1[[1]]

# Use ggsave to save as a high-resolution image
library(ggplot2)
ggsave("MR_scatter_plot.pdf", plot = p1[[1]], width = 7, height = 5, dpi = 300)


### Colocalization scatter comparison plot (LocusCompare Plot)
devtools::install_github("boxiangliu/locuscomparer")
library(locuscomparer)
colnames(gwas_pso)[colnames(gwas_pso) == "SNP"] <- "rsid"
colnames(gwas_pso)[colnames(gwas_pso) == "P"] <- "pval"
colnames(gwas_pneu)[colnames(gwas_pneu) == "SNP"] <- "rsid"
colnames(gwas_pneu)[colnames(gwas_pneu) == "P"] <- "pval"
gwas_pso_clean2 <- gwas_pso %>%
   group_by(rsid) %>%
   slice_min(order_by = pval, n = 1, with_ties = FALSE) %>%
   ungroup()
# 2. Deduplicate pneumonia data: group by rsid, select the row with the minimum pval in each group
gwas_pneu_clean2 <- gwas_pneu %>%
   group_by(rsid) %>%
   slice_min(order_by = pval, n = 1, with_ties = FALSE) %>%
   ungroup()
# Take the intersection of the two datasets by rsid and merge them
merged_data <- inner_join(gwas_pso_clean2, gwas_pneu_clean2, by = "rsid", suffix = c(".pso", ".pneu"))
# Convert P values to -log10 format, which is the standard form for plotting
merged_data <- merged_data %>%
   mutate(
      logP_pso = -log10(pval.pso),
      logP_pneu = -log10(pval.pneu)
   )
# Specify the core SNP to highlight
target_snp <- "rs35741374"
# Plot using ggplot2
p_custom <- ggplot(merged_data, aes(x = logP_pso, y = logP_pneu)) +
   geom_point(color = "darkgray", size = 2, alpha = 0.6) +
   
   geom_point(data = filter(merged_data, rsid == target_snp), 
              color = "red", size = 4, shape = 18) +
   
   geom_text(data = filter(merged_data, rsid == target_snp),
             aes(label = rsid), vjust = -1.5, color = "darkred", fontface = "bold") +
   
   geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "blue", alpha = 0.5) +
   
   theme_bw() +
   labs(
      title = paste("Colocalization at", target_snp, "Region"),
      x = "-log10(P) for Psoriasis",
      y = "-log10(P) for Pneumonia"
   ) +
   theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.title = element_text(face = "bold")
   )

print(p_custom)

### locuszoom
library(ggplot2)
library(dplyr)
library(patchwork)
target_snp <- "rs35741374"
center_info <- gwas_pso[gwas_pso$rsid == target_snp, ]

chr_target <- center_info$CHR[1]
bp_target <- center_info$BP[1]

bp_min <- bp_target - 500000
bp_max <- bp_target + 500000

region_pso <- gwas_pso %>% 
   filter(CHR == chr_target & BP >= bp_min & BP <= bp_max)

region_pneu <- gwas_pneu %>% 
   filter(CHR == chr_target & BP >= bp_min & BP <= bp_max)

target_snp <- "rs35741374"

pso_plot_data <- region_pso %>%
   mutate(BP_Mb = BP / 1000000, logP = -log10(pval))
pneu_plot_data <- region_pneu %>%
   mutate(BP_Mb = BP / 1000000, logP = -log10(pval))

target_pso <- pso_plot_data %>% filter(rsid == target_snp)
target_pneu <- pneu_plot_data %>% filter(rsid == target_snp)

x_limits <- c(min(pso_plot_data$BP_Mb), max(pso_plot_data$BP_Mb))

p1 <- ggplot(pso_plot_data, aes(x = BP_Mb, y = logP)) +
   geom_point(color = "#7f8c8d", size = 2, alpha = 0.6) +
   geom_point(data = target_pso, fill = "#c0392b", color = "black", shape = 23, size = 4) +
   geom_vline(xintercept = target_pso$BP_Mb, linetype = "dotted", color = "#c0392b", size = 0.8) +
   theme_classic() +
   scale_x_continuous(limits = x_limits) +
   labs(
      title = paste("Colocalization at", target_snp, "Region"),
      y = bquote("-log"[10]~"(P) Psoriasis")
   ) +
   theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      axis.title.x = element_blank(), # Hide the X-axis title of the top panel to make the plot more compact
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_text(face = "bold", size = 11, color = "black")
   )

p2 <- ggplot(pneu_plot_data, aes(x = BP_Mb, y = logP)) +
   geom_point(color = "#7f8c8d", size = 2, alpha = 0.6) +
   geom_point(data = target_pneu, fill = "#c0392b", color = "black", shape = 23, size = 4) +
   geom_vline(xintercept = target_pneu$BP_Mb, linetype = "dotted", color = "#c0392b", size = 0.8) +
   theme_classic() +
   scale_x_continuous(limits = x_limits) +
   labs(
      x = paste("Position on Chromosome", pso_plot_data$CHR[1], "(Mb)"),
      y = bquote("-log"[10]~"(P) Pneumonia")
   ) +
   theme(
      axis.title.x = element_text(face = "bold", size = 12),
      axis.title.y = element_text(face = "bold", size = 11, color = "black")
   )

p_combined <- p1 / p2

print(p_combined)

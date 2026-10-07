#!/bin/bash
# ==============================================================================
# Script 4: Linkage Disequilibrium Score Regression (04_LDSC_Genetic_Correlation.sh)
# Description: Formats GWAS summary statistics (munging) and estimates the 
#              genome-wide genetic correlation (Rg) between psoriasis and pneumonia.
# ==============================================================================

# Define paths to the LDSC software suite and the 1000 Genomes European reference panel
# (Users replicating this study should modify these paths to their local installation)
LDSC_DIR="./ldsc"
EUR_LD_DIR="./ldsc/eur_w_ld_chr"

# =======================================================
# Part 1: Format Summary Statistics (Munging)
# =======================================================
echo "Step 1: Munging Psoriasis summary statistics..."
python ${LDSC_DIR}/munge_sumstats.py \
    --sumstats Psoriasis_cleaned_for_LDSC_hg19.txt \
    --out psoriasis_final \
    --snp SNP \
    --a1 A1 \
    --a2 A2 \
    --p P \
    --signed-sumstats BETA,0 \
    --N-cas 36466 \
    --N-con 458078 \
    --chunksize 500000 \
    --merge-alleles ${EUR_LD_DIR}/w_hm3.snplist

# Note: The pneumonia summary statistics should be munged using the same 
# procedure to generate 'pneumonia_final.sumstats.gz' before proceeding to Part 2.
# python ${LDSC_DIR}/munge_sumstats.py \
#     --sumstats Pneumonia_cleaned_for_LDSC_hg19.txt \
#     --out pneumonia_final \
#     ... (parameters omitted for brevity)

# =======================================================
# Part 2: Genetic Correlation Estimation (Rg)
# =======================================================
echo "Step 2: Calculating genome-wide genetic correlation (Rg)..."
python ${LDSC_DIR}/ldsc.py \
    --rg psoriasis_final.sumstats.gz,pneumonia_final.sumstats.gz \
    --ref-ld-chr ${EUR_LD_DIR}/ \
    --w-ld-chr ${EUR_LD_DIR}/ \
    --out pso_pneu_rg

echo "LDSC analysis complete. Results saved in pso_pneu_rg.log"

# Psoriasis-Pneumonia-Analysis
Analytical pipeline for psoriasis and pneumonia risk, integrating UK Biobank observational cohort analysis with genetic causal inference (MR; Coloc) using FinnGen and GWAS Catalog data.

## Data Availability
To protect participant privacy and comply with regulatory requirements, raw individual-level data is not hosted in this repository. 

*   **UK Biobank Data:** Access to the individual-level data used in scripts `01-03` requires an approved application via the [UK Biobank Access Management System (AMS)](https://ams.ukbiobank.ac.uk/). 
*   **GWAS Summary Statistics:** The summary-level data used in scripts `04-05` for genetic causal inference are publicly available and can be accessed as follows:
*   **Psoriasis GWAS:** Data was obtained from the GWAS Catalog (Study Accession: [GCST90472771](https://www.ebi.ac.uk/gwas/studies/GCST90472771)), as published in *Nature Communications* ([DOI: 10.1038/s41467-025-56719-8](https://www.nature.com/articles/s41467-025-56719-8)).
*   **Pneumonia GWAS:** Data was obtained from the [FinnGen Consortium](https://risteys.finngen.fi/endpoints/J10_PNEUMONIA) public data releases.

## Repository Structure & Analytic Pipeline
The analytical pipeline is divided into five sequential R scripts:

### 1. Observational Epidemiology (UK Biobank)
*   **`01_Data_Cleaning.R`**
    *   Merges baseline phenotype data with longitudinal follow-up datasets.
    *   Performs strict handling of generalized missing values and regex-based date validation.
    *   Defines clinical covariates and subgroups (Gender, Age </>= 60, WHO Obesity categories, smoking status).
*   **`02_Cross_Sectional_Logistic.R`**
    *   Implements multivariate Logistic Regression models to evaluate cross-sectional associations between psoriasis exposure and prevalent pneumonia at baseline.
    *   Executes automated dynamic covariate filtering and stratification analyses (by Gender, Age, BMI, and Smoking status).
*   **`03_Prospective_Cox.R`**
    *   Derives the prospective 'at-risk' cohort by rigorously excluding prevalent baseline cases.
    *   Defines incident pneumonia events and calculates person-years of follow-up.
    *   Fits Cox Proportional Hazards models with automated Proportional Hazards (PH) assumption validation using Schoenfeld residuals.

### 2. Genetic Causal Inference (Summary-Level Data)
*   **`04_Mendelian_Randomization.R`**
    *   **Instrument Selection:** Identifies valid instrumental variables (IVs) with strict criteria (P < 5e10-8, F-statistic > 10) and excludes the highly pleiotropic Major Histocompatibility Complex (MHC) region.
    *   **Clumping & Harmonization:** Performs local LD clumping utilizing the 1000 Genomes EUR reference panel via `plinkbinr` and harmonizes effect allele frequencies (EAF).
    *   **MR Analyses & Sensitivity:** Executes TwoSampleMR methods (IVW, MR-Egger, Weighted Median), alongside Cochran's Q test, MR-Egger intercept pleiotropy test, leave-one-out validation, and the Steiger directionality test.
*   **`05_Colocalization.R`**
    *   Employs the `coloc` package to evaluate the probability of shared causal variants between psoriasis and pneumonia.
    *   Utilizes high-performance `data.table` slicing to extract regional summary statistics within a 1Mb window (±500kb) surrounding the identified lead SNPs.

## Software Requirements and Dependencies
The analyses were conducted in the **R statistical computing environment** (version 4.4.2). To run the scripts, the following R packages are required:
*   **Data Wrangling:** `dplyr`, `data.table`, `writexl`
*   **Survival Analysis:** `survival`
*   **Genetic Analysis:** `TwoSampleMR`, `coloc`, `plinkbinr`

## Reproducibility
Researchers with approved UK Biobank access and the downloaded GWAS summary statistics can replicate the findings of this study by placing the respective data files in the working directory and executing the scripts sequentially from `01` to `05`.

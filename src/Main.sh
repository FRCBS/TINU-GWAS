#!/bin/bash

## -----------------------------------------------------------
## GWAS of TIN patients vs. FIN HSCT donors 
## -----------------------------------------------------------


# get a list of common SNPs between imputed IC1 and IC3 FIN data
module load r-env
Rscript ./src/extractCommonSNPs.R

# extract common SNP and merge into a single file
./src/mergeIC.sh

# pick only donors and TIN patients
module load r-env
Rscript ./src/getSamples.R

# create a Plink dosage genotype file of samples to be included
./src/extractSamples.sh

# association test between IC1 and IC3
module load plink/1.90b4.1
plink --bfile ./data/IC1-IC3_imputed.merged --allow-no-sex --assoc --pheno ./data/IC1-IC3_case-control.pheno --out ./data/IC1-IC3_case-control

# pick significant SNPs
module load r-env
Rscript ./src/IC1-IC3-case-control_removeSNPs.R

# remove cohort-associationg SNPs, output raw dosage and bed formats
module load plink/1.90b4.1
plink --bfile ./data/ICFIN_imp --exclude ./data/IC1-IC3_case-control.snps --geno 0.1 --recode --alleleACGT --out ./data/ICFIN_imp_filtered
plink --bfile ./data/ICFIN_imp --exclude ./data/IC1-IC3_case-control.snps --geno 0.1 --make-bed --out ./data/ICFIN_imp_filtered

# LD prune filtered IC3 FIN data, run PCA
plink --bfile ./data/ICFIN_imp_filtered --indep-pairwise 50 5 0.8 --out ./data/ICFIN_imp_filtered_LDpruned
plink --bfile ./data/ICFIN_imp_filtered --exclude ./data/ICFIN_imp_filtered_LDpruned.prune.out --make-bed --out ./data/ICFIN_imp_filtered_LDpruned
plink --bfile ./data/ICFIN_imp_filtered_LDpruned --pca 30 header --out ./data/pca

# run SPAtest for assoction
module load r-env
Rscript ./src/runSPAtest.R

# HLA imputation, analyses, plots
Rscript ./src/analysis.R

# Replication of LINGO2 SNPs in FinnGen eye diseases
Rscript ./src/SNPs_FG_replication.R



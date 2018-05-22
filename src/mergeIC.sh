#!/bin/bash

module load plink/1.90b4.1

plink --bfile ./data/IC1_chr1_22_all.impute2_filtered0.5 --extract ./data/common.imputed.snps --make-bed --out ./data/IC1_imp_commonSNPs 

plink --file ./data/IC3F_chr1_22_all.impute2_filtered0.5 --extract ./data/common.imputed.snps --make-bed --out ./data/IC3_imp_commonSNPs 

plink --bfile ./data/IC1_imp_commonSNPs --bmerge ./data/IC3_imp_commonSNPs --make-bed --out ./data/IC1-IC3_imputed.merged



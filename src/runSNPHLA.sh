#!/bin/bash

## extract HLA region and donors+TIN from IC FIN data

module load plink/1.90b4.1

# IC1
cat ./data/IC3F.fam | grep 'DT\|TIN' > ./data/indiv # extract donors and TIN subjects
plink -./data/bfile IC3F --keep ./data/indiv --chr 6 --from-mb 29 --to-mb 33.5 --make-bed --out ./data/IC3F_TIN_HLA
rm -f ./data/indiv

# IC3
cat ./data/IC1_cleaned.fam | grep 'DT\|TIN' > ./data/indiv
# manual step: DT5462 & DT5302 replicates removed with nano
plink -bfile ./data/IC1_cleaned --keep ./data/indiv --chr 6 --from-mb 29 --to-mb 33.5 --make-bed --out ./data/IC1_HLA
rm -f ./data/indiv

# merge IC1 & IC3 FIN HLA
plink --bfile ./data/IC1_HLA --bmerge ./data/IC3F_TIN_HLA --make-bed --allow-no-sex --geno 0.01 --out ./data/IC1_IC3_FIN_HLA


## impute HLA types to IC1&IC3 data
module load biokit
./src/SNP2HLA_package_v1.0.3/SNP2HLA/SNP2HLA.csh \
	./data/IC1_IC3_FIN_HLA \
	./src/SNP2HLA_package_v1.0.3/SNP2HLA/T1DGC/T1DGC_REF \
	./data/IC1_IC3_FIN_HLA_SNP2HLA \
	./src/plink 6000 1000 # Plink v1.07


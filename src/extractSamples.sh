#!/bin/bash

module load plink/1.90b4.1

# plink binary output
plink --bfile ./data/IC1-IC3_imputed.merged --keep ./data/samples.list --make-bed --out ./data/ICFIN_imp 


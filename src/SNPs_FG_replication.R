
library(data.table)
library(tidyverse)
library(ggpubr)

## ------------------------------------------------------------
## Data
## ------------------------------------------------------------

# FinnGen association data
fg.pheno <- fread('./data/FG/finngen_R4_endpoint_definition.tsv', data.table=F) %>% 
  filter(., !grepl('diabet|cancer|tumor|neoplasm', LONGNAME, ignore.case=T))

rs188424115 <- fread('./data/FG/9_28362193_A_T_phenotype_associations.tsv', data.table=F) %>% 
  filter(!grepl('diabet|cancer|tumor|neoplasm', phenostring, ignore.case=T), pval<1e-2)
rs192601611 <- fread('./data/FG/9_28512770_G_A_phenotype_associations.tsv', data.table=F) %>% 
  filter(!grepl('diabet|cancer|tumor|neoplasm', phenostring, ignore.case=T), pval<1e-2)

# Eye diseases in LINGO2 replicated variants
snps.fg.replication <- rbind(
  data.frame(SNP='rs192601611', rbind(
    rs192601611[which(grepl('eye', rs192601611$category, ignore.case=T)), ])),
  data.frame(SNP='rs188424115', rbind(
    rs188424115[which(grepl('eye', rs188424115$category, ignore.case=T)), ]))
)

# SNP metadata from PheWeb
rs192601611.meta <- data.frame(N_cases=c(304, 1201, 1238, 110, 939, 4198),
                               N_controls=c(164762, 164973, 164973, 165083, 164973, 164973),
                               MAF_cases=c(4.2, 1, 1.1, 5, 1.1, 1.6),
                               MAF_controls=c(2, 2, 2, 2, 2, 2))

rs188424115.meta <- data.frame(N_cases=c(366, 1201),
                               N_controls=c(165083, 164973),
                               MAF_cases=c(0.15, 0.83),
                               MAF_controls=c(1.5, 1.5))

snps.fg.replication <- cbind(snps.fg.replication, rbind(rs192601611.meta, rs188424115.meta))


## ------------------------------------------------------------
## Test enrichment against finngen phenotypes
## ------------------------------------------------------------

# count the number of occurrences of terms in FinnGen phenotypes list
fg.pheno.eye     <- fg.pheno[which(grepl('eye', fg.pheno$NAME, ignore.case=T)), ] %>% nrow

# count the number of terms in the SNPs and FinnGen pheno list
fg.rs188424115.phenogroups <- rs188424115$phenostring %>% length
fg.rs192601611.phenogroups <- rs192601611$phenostring %>% length
fg.phenogroups             <- fg.pheno$NAME %>% length

# run hypergeometric test for enrichment of terms and calculate fold enrichment
fg.rs188424115.uveitis.hg <- phyper(snps.fg.replication %>% filter(SNP=='rs188424115') %>% nrow, 
                                    fg.pheno.eye, fg.phenogroups-fg.pheno.eye, fg.rs188424115.phenogroups, F)

fg.rs192601611.uveitis.hg <- phyper(snps.fg.replication %>% filter(SNP=='rs192601611') %>% nrow, 
                                    fg.pheno.eye, fg.phenogroups-fg.pheno.eye, fg.rs192601611.phenogroups, F)

fg.rs192601611.uveitis.hg.2 <- phyper(snps.fg.replication %>% filter(SNP=='rs192601611', beta>0) %>% nrow, 
                                    fg.pheno.eye, fg.phenogroups-fg.pheno.eye, filter(rs192601611, beta>0) %>% nrow, F)

## Write output tables
write.table(
  data.frame(
    SNP=c('rs188424115', 'rs192601611','rs192601611'),
    Beta=c('any', 'any', '>0'),
    Fold_enrichment=c(
      ((snps.fg.replication %>% filter(SNP=='rs188424115') %>% nrow) / fg.rs188424115.phenogroups) /
        (fg.pheno.eye / fg.phenogroups),
      ((snps.fg.replication %>% filter(SNP=='rs192601611') %>% nrow) / fg.rs192601611.phenogroups) /
        (fg.pheno.eye / fg.phenogroups),
      ((snps.fg.replication %>% filter(SNP=='rs192601611', beta>0) %>% nrow) / (filter(rs192601611, beta>0) %>% nrow)) /
        (fg.pheno.eye / fg.phenogroups)
    ), 
    Enrichment_pval=c(fg.rs188424115.uveitis.hg, fg.rs192601611.uveitis.hg, fg.rs192601611.uveitis.hg.2)
  ), './results/FinnGen_SNPs_Eye_enrichment.tsv', sep='\t', row.names=F)

write.table(snps.fg.replication, './results/FinnGen_SNPs_Eye.tsv', sep='\t', row.names=F)



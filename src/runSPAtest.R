## ---------------------------------------------------------------
## TIN(U) patients vs. FIN HSCT donors population 
## SPAtest for association
## ---------------------------------------------------------------


library(SPAtest)
library(data.table)
library(stringr)


## Data

# genotypes
gen <- fread('./data/IC1_imp_commonSNPs_samples.raw', data.table=F)
gen.meta <- gen[ , 1:6]
gen <- gen[, 7:ncol(gen)]

# metadata
met <- read.delim('./data/TIN.metadata', stringsAsFactors=F)
sam <- gen.meta[, 1]
sam.vec <- sam.vec.tin <- sam.vec.tinu <- sam.vec.chronic <- rep(NA, length(sam))
mapfile <- fread('./data/IC3_imp_commonSNPs.bim', header=F, data.table=F)


## Phenotypes

# TIN
sam.vec.tin[ grepl('DT', sam) ] <- 0
sam.vec.tin[ !grepl('DT', sam) ] <- 1

# TINU
sam.vec.tinu[ gsub('SPR.' ,'', sam, fixed=F) %in% met[met[, 'Uveitis']=='y', 1]  ] <- 1
sam.vec.tinu[ grepl('DT', sam) ] <- 0
is.na(sam.vec.tinu)

# Chronic Uveitis
sam.vec.chronic[ gsub('SPR.' ,'', sam, fixed=F) %in% met[met[, 'ChronicUv']=='y', 1]  ] <- 1
sam.vec.chronic[ grepl('DT', sam) ] <- 0
is.na(sam.vec.chronic)


## Association & write output

# TIN
spa.out.tin <- ScoreTest_SPA(genos=t(gen), pheno=sam.vec.tin, beta.out=T, beta.Cutoff=1*10^-6)
spa.out.tin.signif <- which(spa.out.tin$p.value < 5e-8 & spa.out.tin$Is.converge==T)
spa.out.tin.signif.pos <- mapfile[which(mapfile[, 2] %in% gsub('.{2}$', '', colnames(gen)[spa.out.tin.signif])), c(1, 4)] 
write.table(data.frame(spa.out.tin$p.value[spa.out.tin.signif], spa.out.tin$beta[spa.out.tin.signif], spa.out.tin.signif.pos), 
            './data/TIN_SPA_signif', quote=F, sep='\t', col.names=c('PVALUE', 'BETA', 'CHR', 'POS'))
write(na.omit(spa.out.tin$p.value), './data/TIN_SPA_allp')

# TINU
spa.out.tinu <- ScoreTest_SPA(genos=t(gen[!is.na(sam.vec.tinu), ]), pheno=sam.vec.tinu[!is.na(sam.vec.tinu)], beta.out=T, beta.Cutoff=1*10^-6)
spa.out.tinu.signif <- which(spa.out.tinu$p.value < 5e-8 & spa.out.tinu$Is.converge==T)
spa.out.tinu.signif.pos <- mapfile[which(mapfile[, 2] %in% gsub('.{2}$', '', colnames(gen)[spa.out.tinu.signif])), c(1, 4)] 
write.table(data.frame(spa.out.tinu$p.value[spa.out.tinu.signif], spa.out.tinu$beta[spa.out.tinu.signif], spa.out.tinu.signif.pos), 
            './data/TINU_SPA_signif', quote=F, sep='\t', col.names=c('PVALUE', 'BETA', 'CHR', 'POS'))
write(na.omit(spa.out.tinu$p.value), './data/TINU_SPA_allp')

# Chronic Uveitis
spa.out.chronic <- ScoreTest_SPA(genos=t(gen[!is.na(sam.vec.chronic), ]), pheno=sam.vec.chronic[!is.na(sam.vec.chronic)], beta.out=T, beta.Cutoff=1*10^-6)
spa.out.chronic.signif <- which(spa.out.chronic$p.value < 5e-8 & spa.out.chronic$Is.converge==T)
spa.out.chronic.signif.pos <- mapfile[which(mapfile[, 2] %in% gsub('.{2}$', '', colnames(gen)[spa.out.chronic.signif])), c(1, 4)] 
write.table(data.frame(spa.out.chronic$p.value[spa.out.chronic.signif], spa.out.chronic$beta[spa.out.chronic.signif], spa.out.chronic.signif.pos), 
            './data/Chronic_SPA_signif', quote=F, sep='\t', col.names=c('PVALUE', 'BETA', 'CHR', 'POS'))
write(na.omit(spa.out.chronic$p.value), './data/Chronic_SPA_allp')



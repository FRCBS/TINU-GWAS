
## ---------------------------------------------------------------
## Analysis of TIN(U) association results
##
## ---------------------------------------------------------------

# Libs
library(SPAtest)
library(data.table)
library(biomaRt)
library(tidyverse)
library(ggpubr)
library(HIBAG)
library(RColorBrewer)




## --------------------------------------------------------
## Data
## --------------------------------------------------------

# TIN
SPA.TIN <- fread('./data/TIN_SPA', data.table=F, header=F, skip=1)
colnames(SPA.TIN) <- c('PVALUE', 'BETA', 'CHR', 'POS', 'SNP')
SPA.TIN$SNP <- str_split_fixed(SPA.TIN$SNP, ':', 5)[, 1]
SPA.TIN$SNP[!grepl('rs', SPA.TIN$SNP)] <- NA
SPA.TIN.lingo.ind <- which(SPA.TIN$CHR==9 & SPA.TIN$POS>28e6 & SPA.TIN$POS<29.5e6)
SPA.TIN.mhc.ind <- which(SPA.TIN$CHR==6 & SPA.TIN$POS>29.5e6 & SPA.TIN$POS<33.5e6)

# TIN + uveitis
SPA.TINU <- fread('./data/TINU_SPA', data.table=F, header=F, skip=1)
colnames(SPA.TINU) <- c('PVALUE', 'BETA', 'CHR', 'POS', 'SNP')
SPA.TINU$SNP <- str_split_fixed(SPA.TINU$SNP, ':', 5)[, 1]
SPA.TINU$SNP[!grepl('rs', SPA.TINU$SNP)] <- NA
SPA.TINU.lingo.ind <- which(SPA.TINU$CHR==9 & SPA.TINU$POS>27.8e6 & SPA.TINU$POS<29.3e6)
SPA.TINU.mhc.ind <- which(SPA.TINU$CHR==6 & SPA.TINU$POS>29.5e6 & SPA.TINU$POS<33.5e6)

# TIN + chronic uveitis
SPA.TINUC <- fread('./data/TINC_SPA', data.table=F, header=F, skip=1)
colnames(SPA.TINUC) <- c('PVALUE', 'BETA', 'CHR', 'POS', 'SNP')
SPA.TINUC$SNP <- str_split_fixed(SPA.TINUC$SNP, ':', 5)[, 1]
SPA.TINUC$SNP[!grepl('rs', SPA.TINUC$SNP)] <- NA
SPA.TINUC.lingo.ind <- which(SPA.TINUC$CHR==9 & SPA.TINUC$POS>28e6 & SPA.TINUC$POS<29.5e6)
SPA.TINUC.mhc.ind <- which(SPA.TINUC$CHR==6 & SPA.TINUC$POS>29.5e6 & SPA.TINUC$POS<33.5e6)



## --------------------------------------------------------
## blood cell expression profile of LINGO2
## --------------------------------------------------------

meltList <- function(l) {
  tmp <- sapply(l, length)
  nam <- rep(names(tmp), tmp)
  tmp <- unlist(l)
  data.frame(Var=nam, Value=tmp)
}

bce <- fread('./data/HaemAtlasMKEBNormalizedIntensities.csv', data.table=F)[-1, ] %>% data.frame
bce.arr <- fread('./data/A-MEXP-930.adf.txt', skip=20, data.table=F)
bce.sam <- fread('./data/E-TABM-633.sdrf.txt', data.table=F)
colnames(bce)[2:ncol(bce)] <- gsub('X', '', colnames(bce)[2:ncol(bce)])
bce.sam <- bce.sam[match(colnames(bce)[2:ncol(bce)], bce.sam[, 'Hybridization Name']), ]

bce.lingo.ind <- which(bce[, 1]==bce.arr[grepl('LINGO2', bce.arr[, 3]), 1])

bce.mean <- mean(bce[, 2:ncol(bce)] %>% unlist %>% as.numeric)

bce.lingo.cell <- tapply(2:ncol(bce), bce.sam[, 7], function(x) bce[bce.lingo.ind, x] %>% as.numeric) %>% meltList

pdf('./results/ImmuneCells_LINGO2_expression.pdf', height=5)
ggplot(bce.lingo.cell, aes(Var %>% factor, Value)) +
  geom_boxplot() +
  xlab('Immune cell type') +
  ylab('Normalized expression') +
  geom_hline(yintercept=bce.mean, linetype='dashed') +
  ggtitle('LINGO2 mRNA') +
  theme_light() + 
  theme(axis.text.x=element_text(angle=45,hjust=1))
dev.off()


## --------------------------------------------------------
## Manhattans
## --------------------------------------------------------

                         
source('./scripts/GGManhattan.R')
# SNP CHR BP P

# Whole genome

tmp <- SPA.TIN
colnames(tmp) <- c('P', 'BETA', 'CHR', 'BP', 'SNP')
tmp <- tmp[, c('SNP', 'CHR', 'BP', 'P')]

jpeg('./results/TIN_genomic.jpg', height=3.2, width=9, units='in', type='cairo', res=800)
gg.manhattan(tmp %>% na.omit, threshold=filter(tmp, P<5e-8) %>% arrange(P) %>% .[1:3,] %>% .$P %>% max, 
             hlight=NA, sig=5e-8, 
             sugg=max(tmp$P[which((tmp$P %>% p.adjust(., method='BH'))<0.05)]), pointsize=1,
             col=brewer.pal(4, 'Blues')[3:4], ylims=c(0, 9), title='TIN with or without uveitis')
dev.off()


tmp <- SPA.TINU
colnames(tmp) <- c('P', 'BETA', 'CHR', 'BP', 'SNP')
tmp <- tmp[, c('SNP', 'CHR', 'BP', 'P')]

jpeg('./results/TINU_genomic.jpg', height=3.2, width=9, units='in', type='cairo', res=800)
gg.manhattan(tmp %>% na.omit, threshold=filter(tmp, P<5e-8) %>% arrange(P) %>% .[1:3,] %>% .$P %>% max, 
             hlight=NA, sig=5e-8, 
             sugg=max(tmp$P[which((tmp$P %>% p.adjust(., method='BH'))<0.05)]), pointsize=1,
             col=brewer.pal(4, 'Blues')[3:4], ylims=c(0, 9), title='TIN with uveitis')
dev.off()


tmp <- SPA.TINUC
colnames(tmp) <- c('P', 'BETA', 'CHR', 'BP', 'SNP')
tmp <- tmp[, c('SNP', 'CHR', 'BP', 'P')]

jpeg('./results/TINUC_genomic.jpg', height=3.2, width=9, units='in', type='cairo', res=800)
gg.manhattan(tmp %>% na.omit, threshold=filter(tmp, P<5e-8) %>% arrange(P) %>% .[1:3,] %>% .$P %>% max, 
             hlight=NA, sig=5e-8, 
             sugg=max(tmp$P[which((tmp$P %>% p.adjust(., method='BH'))<0.05)]), pointsize=1,
             col=brewer.pal(4, 'Blues')[3:4], ylims=c(0, 9), title='TIN with chronic uveitis')
dev.off()



## --------------------------------------------------------
## Blood cis-eQTL
## --------------------------------------------------------

lingo.cis <- fread('./data/LINGO2_cis', data.table=F)
lingo.cis <- arrange(lingo.cis, V4)
lingo.cis <- filter(lingo.cis, V5>0, V4>27.8e6)

ensembl <- useMart(biomart="ENSEMBL_MART_FUNCGEN", host="grch37.ensembl.org", path="/biomart/martservice", 
                     dataset="hsapiens_regulatory_feature")

lingo2.ensembl <- getBM(attributes=c('activity', 'feature_type_description', 'chromosome_start','chromosome_end', 
                                     'feature_type_name', 'regulatory_stable_id') , 
                           filters='chromosomal_region', values='9:27800000:29800000', mart=ensembl)
lingo2.ensembl <- lingo2.ensembl %>% filter(., activity=='ACTIVE') 
lingo2.ensembl <- filter(lingo2.ensembl, chromosome_start>27.8e6, chromosome_start<29.3e6)
lingo2.ensembl$feature_type_name <- gsub(' ', '_', lingo2.ensembl$feature_type_name)
lingo2.ensembl$feature_type_name %>% table

atrack1 <- AnnotationTrack(start=lingo2.ensembl$chromosome_start, end=lingo2.ensembl$chromosome_end, 
                           feature=lingo2.ensembl$feature_type_name, id=lingo2.ensembl$regulatory_stable_id,
                           chromosome='chr9', genome='hg19', name="REG", stacking='dense', 
                           Promoter_Flanking_Region='grey30', TF_binding_site='grey30',
                           CTCF_Binding_Site='grey90', Enhancer='grey90', Open_chromatin='grey90',
                           background.title="grey80", fontface.axis=1, fontface.title=2)
 
dtrack1 <- DataTrack(data=lingo.cis$V1 %>% -log10(.), start=lingo.cis$V4, end=lingo.cis$V4, 
                     chromosome='chr9', genome='hg19', name="eQTL", type='h', 
                     background.title="grey80", fontface.axis=1, fontface.title=2)

dtrack2 <- DataTrack(data=SPA.TINU[SPA.TINU.lingo.ind, ]$PVALUE %>% -log10(.), 
                     start=SPA.TINU[SPA.TINU.lingo.ind, ]$POS, end=SPA.TINU[SPA.TINU.lingo.ind, ]$POS, 
                     chromosome='chr9', genome='hg19', name="TINU GWAS", type='h', aggregation='max', window=100,
                     baseline=-log10(9.151029e-06), col.baseline='skyblue1', lty.baseline='dashed',
                     background.title="grey80", fontface.axis=1, fontface.title=2)

biomTrack <- BiomartGeneRegionTrack(genome="hg19", symbol='LINGO2', name="ENSEMBL Gene",
                                    transcriptAnnotation="symbol", collapse=F,
                                    pseudogene='red', utr5='orange', non_coding='cyan', min.height=3,
                                    background.title="grey80", fontface.axis=1, fontface.title=2)
biomTrack %>% feature %>% unique

cairo_pdf('./results/eQTL_tracks.pdf', height=4, width=8)
plotTracks(list(IdeogramTrack(genome='hg19', chromosome='chr9'),
                GenomeAxisTrack(), dtrack1, dtrack2,
                biomTrack, atrack1),
           col=NULL, col.line=NULL, sizes=c(1, 1.5, 2.5, 2.5, 4.5, 1), from=27.83e6, to=29.3e6)
dev.off()


## --------------------------------------------------------
## LINGO2 region haplotype
## --------------------------------------------------------

# genotype data
ped <- fread('./data/IC3_TINU_LINGO.ped', data.table=F)
ped.gen <- ped[, 7:ncol(ped)]
ped.gen <- t(sapply(seq(1, ncol(ped.gen), by=2), function(x) apply(cbind(ped.gen[, x], ped.gen[, x+1]), 1, paste, collapse=' ')))
ped.met <- ped[, 1:6]
map <- fread('./data/IC3_TINU_LINGO.map', data.table=F)
map[, 2] <- sapply(map[, 2], function(x) strsplit(x, ':', fixed=T)[[1]][1])
rownames(ped.gen) <- map[, 2]
colnames(ped.gen) <- gsub('SPR.', '', ped.met[, 1], fixed=T)

# check sample order
(ped.met[, 1]==tin.pheno[, 1]) %>% all

# plotting

# genotype into tsygocity numerical coding
geno2num <- function(prow) { 
  tt <- names(sort(table(prow)))
  prow[prow==tt[1]] <- 'minor'
  prow[prow==tt[2]] <- 'heterozygote'
  prow[prow==tt[3]] <- 'major'
  prow
}

ped.gen.num <- t(apply(ped.gen, 1, geno2num))
colnames(ped.gen.num) <- colnames(ped.gen)#sapply(colnames(ped.gen), str_pad, side='right', width=6)
ped.gen.num <- data.frame(t(ped.gen.num) %>% melt, ChronicUv=tin.pheno$ChronicUv, stringsAsFactors=F)
ped.gen.num$Var1 <- gsub('TINU', '', ped.gen.num$Var1)
ped.gen.num$ChronicUv[ped.gen.num$ChronicUv=='y'] <- 'with chronic uveitis'
ped.gen.num$ChronicUv[ped.gen.num$ChronicUv=='n'] <- 'without chronic uveitis'

cairo_pdf(height=5.5, width=5, file='./results/LINGO2_SNPs_haplo.pdf')
ggplot(ped.gen.num[nrow(ped.gen.num):1, ], 
       aes(Var1, Var2, fill=value %>% factor(., levels=c('minor', 'heterozygote', 'major')))) +
  geom_tile() +
  facet_wrap(~ChronicUv, scales='free_x') +
  xlab('TIN(U) patients') +
  ylab('LINGO2 SNPs') +
  theme_light() +
  theme(axis.text.x=element_text(angle=90, hjust=1, size=7), panel.grid=element_blank()) +
  theme(strip.background=element_rect(fill="white"), legend.position="right") +
  theme(strip.text=element_text(colour='grey20', size=11)) +
  theme(legend.position = "bottom") +
  guides(fill=guide_legend(title="Genotype"))
dev.off()




## --------------------------------------------------------
## HLA imputation with HIBAG
## --------------------------------------------------------

# genotype data
jia.mhc <- hlaBED2Geno(bed.fn="./data/VPU_ILLUMINA_AUG_2018_MHC.bed", 
                       fam.fn="./data/VPU_ILLUMINA_AUG_2018_MHC.fam", 
                       bim.fn="./data/VPU_ILLUMINA_AUG_2018_MHC.bim", assembly='hg19')
tin.mhc <- hlaBED2Geno(bed.fn="./data/IC1_IC3_FIN_HLA.bed", #ICFIN_imp_filtered_HLA.bed", 
                       fam.fn="./data/IC1_IC3_FIN_HLA.fam", #ICFIN_imp_filtered_HLA.fam", 
                       bim.fn="./data/IC1_IC3_FIN_HLA.bim", #ICFIN_imp_filtered_HLA.bim", 
                       assembly='hg19')
# default european reference
model.list <- get(load('./data/European-HLA4-hg19.RData'))
model.a    <- hlaModelFromObj(model.list[['A']])
model.b    <- hlaModelFromObj(model.list[['B']])
model.c    <- hlaModelFromObj(model.list[['C']])
model.drb1 <- hlaModelFromObj(model.list[['DRB1']])
model.dqa1 <- hlaModelFromObj(model.list[['DQA1']])
model.dqb1 <- hlaModelFromObj(model.list[['DQB1']])
model.dpb1 <- hlaModelFromObj(model.list[['DPB1']])

# HLA imputation
tin.mhc.a    <- predict(model.a,    tin.mhc, type="response+prob", match.type="Position")
tin.mhc.b    <- predict(model.b,    tin.mhc, type="response+prob", match.type="Position")
tin.mhc.c    <- predict(model.c,    tin.mhc, type="response+prob", match.type="Position")
tin.mhc.drb1 <- predict(model.drb1, tin.mhc, type="response+prob", match.type="Position")
tin.mhc.dqa1 <- predict(model.dqa1, tin.mhc, type="response+prob", match.type="Position")
tin.mhc.dqb1 <- predict(model.dqb1, tin.mhc, type="response+prob", match.type="Position")
tin.mhc.dpb1 <- predict(model.dpb1, tin.mhc, type="response+prob", match.type="Position")

# TIN phenotypes
tin.pheno    <- fread('./data/TIN-koodit', data.table=F)[, c(1,4,5)]
tin.pheno$V1 <- paste0('SPR.', tin.pheno$V1)
tin.pheno    <- tin.pheno[match(tin.mhc$sample.id[grepl('TINU', tin.mhc$sample.id)], tin.pheno$V1), ]
tin.upos     <- tin.uneg <- tin.all <- tin.cpos <- tin.cneg <- rep(NA, length(tin.mhc$sample.id))

tin.all[tin.mhc$sample.id %in% tin.pheno$V1] <- 1
tin.all[!(tin.mhc$sample.id %in% tin.pheno$V1)] <- 0

tin.uneg[tin.mhc$sample.id %in% tin.pheno$V1[tin.pheno[,2]=='n']] <- 1
tin.uneg[!(tin.mhc$sample.id %in% tin.pheno$V1)] <- 0

tin.upos[tin.mhc$sample.id %in% tin.pheno$V1[tin.pheno[,2]=='y']] <- 1
tin.upos[!(tin.mhc$sample.id %in% tin.pheno$V1)] <- 0

tin.cpos[tin.mhc$sample.id %in% tin.pheno$V1[tin.pheno[,3]=='y']] <- 1
tin.cpos[!(tin.mhc$sample.id %in% tin.pheno$V1)] <- 0

tin.cneg[tin.mhc$sample.id %in% tin.pheno$V1[tin.pheno[,3]=='n']] <- 1
tin.cneg[!(tin.mhc$sample.id %in% tin.pheno$V1)] <- 0


# write out imputed HLA types of cases

# TIN all cases
tin.all.types <- cbind(tin.mhc.a$value[tin.all==1, 1:4], tin.mhc.b$value[tin.all==1, 2:4], tin.mhc.c$value[tin.all==1, 2:4],
                       tin.mhc.drb1$value[tin.all==1, 2:4], tin.mhc.dqa1$value[tin.all==1, 2:4], 
                       tin.mhc.dqb1$value[tin.all==1, 2:4], tin.mhc.dpb1$value[tin.all==1, 2:4])
colnames(tin.all.types) <- c('Sample ID', paste(rep(c('A', 'B', 'C', 'DRB1', 'DQA1', 'DQB1', 'DPB1'), each=3), 
                                                colnames(tin.all.types)[-1], sep='_'))
tin.all.types[, colnames(tin.all.types) %>% grepl('prob', .)] <- tin.all.types[, colnames(tin.all.types) %>% 
                                                                                 grepl('prob', .)] %>% round(., 3)

# TIN uveitis cases
tin.uve.types <- cbind(tin.mhc.a$value[tin.upos==1, 1:4], tin.mhc.b$value[tin.upos==1, 2:4], tin.mhc.c$value[tin.upos==1, 2:4],
                       tin.mhc.drb1$value[tin.upos==1, 2:4], tin.mhc.dqa1$value[tin.upos==1, 2:4], 
                       tin.mhc.dqb1$value[tin.upos==1, 2:4], tin.mhc.dpb1$value[tin.upos==1, 2:4]) %>% na.omit
colnames(tin.uve.types) <- c('Sample ID', paste(rep(c('A', 'B', 'C', 'DRB1', 'DQA1', 'DQB1', 'DPB1'), each=3), 
                                                colnames(tin.uve.types)[-1], sep='_'))
tin.uve.types[, colnames(tin.uve.types) %>% grepl('prob', .)] <- 
  tin.uve.types[, colnames(tin.uve.types) %>% grepl('prob', .)] %>% round(., 3)

# TIN chronic uveitis cases
tin.cuve.types <- cbind(tin.mhc.a$value[tin.cpos==1, 1:4], tin.mhc.b$value[tin.cpos==1, 2:4], tin.mhc.c$value[tin.cpos==1, 2:4],
                       tin.mhc.drb1$value[tin.cpos==1, 2:4], tin.mhc.dqa1$value[tin.cpos==1, 2:4], 
                       tin.mhc.dqb1$value[tin.cpos==1, 2:4], tin.mhc.dpb1$value[tin.cpos==1, 2:4]) %>% na.omit
colnames(tin.cuve.types) <- c('Sample ID', paste(rep(c('A', 'B', 'C', 'DRB1', 'DQA1', 'DQB1', 'DPB1'), each=3), 
                                                colnames(tin.cuve.types)[-1], sep='_'))
tin.cuve.types[, colnames(tin.cuve.types) %>% grepl('prob', .)] <- 
  tin.cuve.types[, colnames(tin.cuve.types) %>% grepl('prob', .)] %>% round(., 3)


# Output tables
write.table(tin.all.types, './results/tin.all.types',     sep='\t', row.names=F)
write.table(tin.uve.types, './results/tin.uveitis.types', sep='\t', row.names=F)
write.table(tin.cuve.types, './results/tin.chronic.uveitis.types', sep='\t', row.names=F)


## --------------------------------------------------------
## Fisher association tests
## --------------------------------------------------------

hlaFisher <- function(dat, phe, allele) {
  tmp <- cbind(phe, dat$value[, 2:4])
  tmp[tmp$prob<0.5, 'prob'] <- NA
  tmp <- tmp[complete.cases(tmp), ]
  colnames(tmp)[2:3] <- c('allele', 'allele')
  tmp <- rbind(tmp[, 1:2], tmp[, c(1,3)])
  c(fisher.test(tmp[, 1] %>% factor, (tmp$allele==allele) %>% as.numeric %>% as.factor), 
    Cases=sum(tmp[, 1]), Controls=sum(tmp[, 1]==0),
    Cases.prct=(sum(tmp[tmp$allele==allele, 1])/sum(tmp[, 1]))*100,
    Controls.prct=(sum((tmp[, 1]==0)[tmp$allele==allele])/sum(tmp[, 1]==0))*100)
}

runHLAfisher <- function(d, phe, locus) {
  alleles <- d$value[, 2:3] %>% unlist %>% na.omit %>% table
  alleles <- alleles[alleles>2]
  alleles <- names(alleles)
  tt <- lapply(alleles, function(x) {
    out <- hlaFisher(d, phe, x)
    c(paste0(locus, '*', x), out[[1]], out[[3]], out[[2]][1], out[[2]][2], out[[8]], out[[9]], out[[10]], out[[11]])
  }) %>% do.call(rbind, .)
  colnames(tt) <- c('Allele', 'p-value', 'OR', '95% CI lower', '95% CI upper', 'Cases', 'Controls', 
                    'Allele % cases', 'Allele % controls') 
  data.frame(tt)
} 

# HLA Fisher's test for TIN(U)
write.table(rbind(runHLAfisher(tin.mhc.a, tin.all, 'A'), 
                  runHLAfisher(tin.mhc.b, tin.all, 'B'),
                  runHLAfisher(tin.mhc.c, tin.all, 'C'), 
                  runHLAfisher(tin.mhc.drb1, tin.all, 'DRB1'), 
                  runHLAfisher(tin.mhc.dqa1, tin.all, 'DQA1'), 
                  runHLAfisher(tin.mhc.dqb1, tin.all, 'DQB1'), 
                  runHLAfisher(tin.mhc.dpb1, tin.all, 'DPB1')), 
            './results/TIN_HLA_Fisher.tsv', row.names=F, sep='\t') 

write.table(rbind(runHLAfisher(tin.mhc.a, tin.upos, 'A'), 
                  runHLAfisher(tin.mhc.b, tin.upos, 'B'),
                  runHLAfisher(tin.mhc.c, tin.upos, 'C'), 
                  runHLAfisher(tin.mhc.drb1, tin.upos, 'DRB1'), 
                  runHLAfisher(tin.mhc.dqa1, tin.upos, 'DQA1'), 
                  runHLAfisher(tin.mhc.dqb1, tin.upos, 'DQB1'), 
                  runHLAfisher(tin.mhc.dpb1, tin.upos, 'DPB1')), 
            './results/TINU_HLA_Fisher.tsv', row.names=F, sep='\t') 

write.table(rbind(runHLAfisher(tin.mhc.a, tin.cpos, 'A'), 
                  runHLAfisher(tin.mhc.b, tin.cpos, 'B'),
                  runHLAfisher(tin.mhc.c, tin.cpos, 'C'), 
                  runHLAfisher(tin.mhc.drb1, tin.cpos, 'DRB1'), 
                  runHLAfisher(tin.mhc.dqa1, tin.cpos, 'DQA1'), 
                  runHLAfisher(tin.mhc.dqb1, tin.cpos, 'DQB1'), 
                  runHLAfisher(tin.mhc.dpb1, tin.cpos, 'DPB1')), 
            './results/TINUC_HLA_Fisher.tsv', row.names=F, sep='\t') 

write.table(rbind(runHLAfisher(tin.mhc.a, tin.uneg, 'A'), 
                  runHLAfisher(tin.mhc.b, tin.uneg, 'B'),
                  runHLAfisher(tin.mhc.c, tin.uneg, 'C'), 
                  runHLAfisher(tin.mhc.drb1, tin.uneg, 'DRB1'), 
                  runHLAfisher(tin.mhc.dqa1, tin.uneg, 'DQA1'), 
                  runHLAfisher(tin.mhc.dqb1, tin.uneg, 'DQB1'), 
                  runHLAfisher(tin.mhc.dpb1, tin.uneg, 'DPB1')), 
            './results/TINUnegative_HLA_Fisher.tsv', row.names=F, sep='\t')  

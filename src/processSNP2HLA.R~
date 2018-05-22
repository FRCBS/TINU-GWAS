
## Process results of HLA imputation;
## SNP2HLA T1DGC imputed IC1-IC3 HLA types


## functions

cbindHLA <- function(x, y) { # collect hla types into a matrix
  all.names <- c("A", "B", "C", "DPA1", "DPB1", "DQA1", "DQB1", "DRB1")
  max.l <- max(length(x), length(y))
  dat <- cbind(rep(NA, length(all.names)), rep(NA, length(all.names)))
  rownames(dat) <- all.names
  dat[match(names(x), all.names), 1] <- x 
  dat[match(names(y), all.names), 2] <- y 
  dat
}

freqTable <- function(x) table(x)/sum(table(x))


## data

pha <- fread('./data/IC1_IC3_FIN_HLA_SNP2HLA.bgl.phased', skip=1, header=T, data.table=F) # read SNP2HLA output
pha <- pha[grepl('HLA_', pha[, 'id']), ] # select imputed HLA types
rownames(pha) <- 1:nrow(pha)
pha.hla <- str_split_fixed(pha[, 'id'], '_', 3)
pha <- pha[, 3:ncol(pha)]
pha.samples <- gsub('.1', '', colnames(pha), fixed=T)
pha4 <- pha[nchar(pha.hla[, 3])==4, ] # 4-digit type

# parse 4-digit resolution output into type-by-sample matrix
pha.types.4 <- tapply(1:ncol(pha4), pha.samples, function(x) {
  tmp1 <- pha.hla[as.numeric(rownames(pha4[, x][pha4[, x][, 1]=='P', ])), ]
  tmp1 <- tapply(tmp1[, 3], tmp1[, 2], paste, collapse='/')
  tmp2 <- pha.hla[as.numeric(rownames(pha4[, x][pha4[, x][, 2]=='P', ])), ]
  tmp2 <- tapply(tmp2[, 3], tmp2[, 2], paste, collapse='/')
  
  tmp3 <- cbindHLA(tmp1, tmp2)
  colnames(tmp3) <-  colnames(pha4[, x])
  tmp3
})
pha.types.4 <- do.call(cbind, pha.types.4)

tin.cli <- read.delim('./data/TIN-codes', stringsAsFactors=F)
tin.uveitis <- paste('SPR.', rownames(tin.cli[tin.cli[, 'Uveitis']=='y', ]), sep='')
tin.uveitis <- c(tin.uveitis, paste(tin.uveitis, '.1', sep=''))
tin.uveitis.no <- paste('SPR.', rownames(tin.cli[tin.cli[, 'Uveitis']=='n', ]), sep='')
tin.uveitis.no <- c(tin.uveitis.no, paste(tin.uveitis.no, '.1', sep=''))

# ouput HLA type table

hla.types.out <- pha.types.4
colnames(hla.types.out) <- gsub('SPR.DT|SPR.SPR.DT', 'Control-', colnames(hla.types.out))
colnames(hla.types.out)[colnames(hla.types.out) %in% tin.uveitis] <- 
  gsub('SPR.', '', colnames(hla.types.out)[colnames(hla.types.out) %in% tin.uveitis])
colnames(hla.types.out)[colnames(hla.types.out) %in% tin.uveitis.no] <- 
  gsub('SPR.TINU', 'TIN', colnames(hla.types.out)[colnames(hla.types.out) %in% tin.uveitis.no])
write.table(hla.types.out %>% t, './data/HLA_types_table.txt', quote=F, sep='\t')


## check type differences between TINU and TIN

# A
tmp1 <- table(pha.types.4[1, colnames(pha.types.4) %in% tin.uveitis])
tmp2 <- table(pha.types.4[1, colnames(pha.types.4) %in% tin.uveitis.no])
tmp3 <- cbind(CASE=tmp1, CONTROL=tmp2[match(names(tmp1), names(tmp2))]) %>% na.omit
fisher.test(cbind(colSums(tmp3)-tmp3[3, ], tmp3[3, ])) 
fisher.test(cbind(colSums(tmp3)-tmp3[2, ], tmp3[2, ])) 
# B
tmp1 <- table(pha.types.4[2, colnames(pha.types.4) %in% tin.uveitis])
tmp2 <- table(pha.types.4[2, colnames(pha.types.4) %in% tin.uveitis.no])
tmp3 <- cbind(CASE=tmp1, CONTROL=tmp2[match(names(tmp1), names(tmp2))]) %>% na.omit
fisher.test(cbind(colSums(tmp3)-tmp3[3, ], tmp3[3, ])) 
fisher.test(cbind(colSums(tmp3)-tmp3[2, ], tmp3[2, ])) 
# C
tmp1 <- table(pha.types.4[3, colnames(pha.types.4) %in% tin.uveitis])
tmp2 <- table(pha.types.4[3, colnames(pha.types.4) %in% tin.uveitis.no])
tmp3 <- cbind(CASE=tmp1, CONTROL=tmp2[match(names(tmp1), names(tmp2))]) %>% na.omit
fisher.test(cbind(colSums(tmp3)-tmp3[3, ], tmp3[3, ])) 
fisher.test(cbind(colSums(tmp3)-tmp3[2, ], tmp3[2, ])) 
fisher.test(cbind(colSums(tmp3)-tmp3[4, ], tmp3[3, ])) 
fisher.test(cbind(colSums(tmp3)-tmp3[5, ], tmp3[2, ])) 

# DPB1
tmp1 <- table(pha.types.4[5, colnames(pha.types.4) %in% tin.uveitis])
tmp2 <- table(pha.types.4[5, colnames(pha.types.4) %in% tin.uveitis.no])
tmp3 <- cbind(CASE=tmp1, CONTROL=tmp2[match(names(tmp1), names(tmp2))]) %>% na.omit
fisher.test(cbind(colSums(tmp3)-tmp3[3, ], tmp3[3, ])) 
fisher.test(cbind(colSums(tmp3)-tmp3[2, ], tmp3[2, ])) 
# DQA1
tmp1 <- table(pha.types.4[6, colnames(pha.types.4) %in% tin.uveitis])
tmp2 <- table(pha.types.4[6, colnames(pha.types.4) %in% tin.uveitis.no])
tmp3 <- cbind(CASE=tmp1, CONTROL=tmp2[match(names(tmp1), names(tmp2))]) %>% na.omit
fisher.test(cbind(colSums(tmp3)-tmp3[3, ], tmp3[3, ])) 
fisher.test(cbind(colSums(tmp3)-tmp3[2, ], tmp3[2, ])) 
# DQB1
tmp1 <- table(pha.types.4[7, colnames(pha.types.4) %in% tin.uveitis])
tmp2 <- table(pha.types.4[7, colnames(pha.types.4) %in% tin.uveitis.no])
tmp3 <- cbind(CASE=tmp1, CONTROL=tmp2[match(names(tmp1), names(tmp2))]) %>% na.omit
fisher.test(cbind(colSums(tmp3)-tmp3[3, ], tmp3[3, ])) 
fisher.test(cbind(colSums(tmp3)-tmp3[2, ], tmp3[2, ])) 
# DRB1
tmp1 <- table(pha.types.4[8, colnames(pha.types.4) %in% tin.uveitis])
tmp2 <- table(pha.types.4[8, colnames(pha.types.4) %in% tin.uveitis.no])
tmp3 <- cbind(CASE=tmp1, CONTROL=tmp2[match(names(tmp1), names(tmp2))]) %>% na.omit
fisher.test(cbind(colSums(tmp3)-tmp3[3, ], tmp3[3, ])) 
fisher.test(cbind(colSums(tmp3)-tmp3[2, ], tmp3[2, ])) 


## check type differences between TIN(U) and controls

# 4-digit resolution
# DPB1
tmp1 <- table(pha.types.4[5, grepl('TIN', colnames(pha.types.4))])
tmp2 <- table(pha.types.4[5, grepl('TIN', colnames(pha.types.4))==F])
tmp3 <- cbind(CASE=tmp1, CONTROL=tmp2[match(names(tmp1), names(tmp2))])
dpb1.4 <- fisher.test(cbind(colSums(tmp3)-tmp3[3, ], tmp3[3, ])) # 0301
tmp1 <- freqTable(pha.types.4[5, grepl('TIN', colnames(pha.types.4))])
tmp2 <- freqTable(pha.types.4[5, grepl('TIN', colnames(pha.types.4))==F])
tmp3 <- cbind(CASE=tmp1, CONTROL=tmp2[match(names(tmp1), names(tmp2))])
dpb1.4.fq <- tmp3[3, ]
# DQA1
tmp1 <- table(pha.types.4[6, grepl('TIN', colnames(pha.types.4))])
tmp2 <- table(pha.types.4[6, grepl('TIN', colnames(pha.types.4))==F])
tmp3 <- cbind(CASE=tmp1, CONTROL=tmp2[match(names(tmp1), names(tmp2))])
dqa1.4 <- fisher.test(cbind(colSums(tmp3)-tmp3[5, ], tmp3[5, ])) # 0401
tmp1 <- freqTable(pha.types.4[6, grepl('TIN', colnames(pha.types.4))])
tmp2 <- freqTable(pha.types.4[6, grepl('TIN', colnames(pha.types.4))==F])
tmp3 <- cbind(CASE=tmp1, CONTROL=tmp2[match(names(tmp1), names(tmp2))])
dqa1.4.fq <- tmp3[5, ]
# DQB1
tmp1 <- table(pha.types.4[7, grepl('TIN', colnames(pha.types.4))])
tmp2 <- table(pha.types.4[7, grepl('TIN', colnames(pha.types.4))==F])
tmp3 <- cbind(CASE=tmp1, CONTROL=tmp2[match(names(tmp1), names(tmp2))])
dqb1.4 <- fisher.test(cbind(colSums(tmp3)-tmp3[5, ], tmp3[5, ])) # 0402
tmp1 <- freqTable(pha.types.4[7, grepl('TIN', colnames(pha.types.4))])
tmp2 <- freqTable(pha.types.4[7, grepl('TIN', colnames(pha.types.4))==F])
tmp3 <- cbind(CASE=tmp1, CONTROL=tmp2[match(names(tmp1), names(tmp2))])
dqb1.4.fq <- tmp3[5, ]
# DRB1
tmp1 <- table(pha.types.4[8, grepl('TIN', colnames(pha.types.4))])
tmp2 <- table(pha.types.4[8, grepl('TIN', colnames(pha.types.4))==F])
tmp3 <- cbind(CASE=tmp1, CONTROL=tmp2[match(names(tmp1), names(tmp2))])
drb1.4 <- fisher.test(cbind(colSums(tmp3)-tmp3[3, ], tmp3[3, ])) # 0802
tmp1 <- freqTable(pha.types.4[8, grepl('TIN', colnames(pha.types.4))])
tmp2 <- freqTable(pha.types.4[8, grepl('TIN', colnames(pha.types.4))==F])
tmp3 <- cbind(CASE=tmp1, CONTROL=tmp2[match(names(tmp1), names(tmp2))])
drb1.4.fq <- tmp3[3, ]


# 4-digit resolution for TINU only
# DPB1
tmp1 <- table(pha.types.4[5, match(tin.uveitis, colnames(pha.types.2))])
tmp2 <- table(pha.types.4[5, grepl('DT', colnames(pha.types.4))])
tmp3 <- cbind(CASE=tmp1, CONTROL=tmp2[match(names(tmp1), names(tmp2))])
dpb1.4.u <- fisher.test(cbind(colSums(tmp3)-tmp3[3, ], tmp3[3, ])) # 0301
tmp1 <- freqTable(pha.types.4[5, match(tin.uveitis, colnames(pha.types.2))])
tmp2 <- freqTable(pha.types.4[5, grepl('DT', colnames(pha.types.4))])
tmp3 <- cbind(CASE=tmp1, CONTROL=tmp2[match(names(tmp1), names(tmp2))])
dpb1.4.fq.u <- tmp3[3, ]
# DQA1
tmp1 <- table(pha.types.4[6, match(tin.uveitis, colnames(pha.types.2))])
tmp2 <- table(pha.types.4[6, grepl('DT', colnames(pha.types.4))])
tmp3 <- cbind(CASE=tmp1, CONTROL=tmp2[match(names(tmp1), names(tmp2))])
dqa1.4.u <- fisher.test(cbind(colSums(tmp3)-tmp3[5, ], tmp3[5, ])) # 0401
tmp1 <- freqTable(pha.types.4[6, match(tin.uveitis, colnames(pha.types.2))])
tmp2 <- freqTable(pha.types.4[6, grepl('DT', colnames(pha.types.4))])
tmp3 <- cbind(CASE=tmp1, CONTROL=tmp2[match(names(tmp1), names(tmp2))])
dqa1.4.fq.u <- tmp3[5, ]
# DQB1
tmp1 <- table(pha.types.4[7, match(tin.uveitis, colnames(pha.types.2))])
tmp2 <- table(pha.types.4[7, grepl('DT', colnames(pha.types.4))])
tmp3 <- cbind(CASE=tmp1, CONTROL=tmp2[match(names(tmp1), names(tmp2))])
dqb1.4.u <- fisher.test(cbind(colSums(tmp3)-tmp3[4, ], tmp3[4, ])) # 0402
tmp1 <- freqTable(pha.types.4[7, match(tin.uveitis, colnames(pha.types.2))])
tmp2 <- freqTable(pha.types.4[7, grepl('DT', colnames(pha.types.4))])
tmp3 <- cbind(CASE=tmp1, CONTROL=tmp2[match(names(tmp1), names(tmp2))])
dqb1.4.fq.u <- tmp3[4, ]
# DRB1
tmp1 <- table(pha.types.4[8, match(tin.uveitis, colnames(pha.types.2))])
tmp2 <- table(pha.types.4[8, grepl('DT', colnames(pha.types.4))])
tmp3 <- cbind(CASE=tmp1, CONTROL=tmp2[match(names(tmp1), names(tmp2))])
drb1.4.u <- fisher.test(cbind(colSums(tmp3)-tmp3[3, ], tmp3[3, ])) # 0802
tmp1 <- freqTable(pha.types.4[8, match(tin.uveitis, colnames(pha.types.2))])
tmp2 <- freqTable(pha.types.4[8, grepl('DT', colnames(pha.types.4))])
tmp3 <- cbind(CASE=tmp1, CONTROL=tmp2[match(names(tmp1), names(tmp2))])
drb1.4.fq.u <- tmp3[3, ]

# 4-digit resolution for TIN only
# DPB1
tmp1 <- table(pha.types.4[5, match(tin.uveitis.no, colnames(pha.types.2))])
tmp2 <- table(pha.types.4[5, grepl('DT', colnames(pha.types.4))])
tmp3 <- cbind(CASE=tmp1, CONTROL=tmp2[match(names(tmp1), names(tmp2))])
dpb1.4.n <- fisher.test(cbind(colSums(tmp3)-tmp3[2, ], tmp3[2, ])) # 0301
tmp1 <- freqTable(pha.types.4[5, match(tin.uveitis.no, colnames(pha.types.2))])
tmp2 <- freqTable(pha.types.4[5, grepl('DT', colnames(pha.types.4))])
tmp3 <- cbind(CASE=tmp1, CONTROL=tmp2[match(names(tmp1), names(tmp2))])
dpb1.4.fq.n <- tmp3[2, ]
# DQA1
tmp1 <- table(pha.types.4[6, match(tin.uveitis.no, colnames(pha.types.2))])
tmp2 <- table(pha.types.4[6, grepl('DT', colnames(pha.types.4))])
tmp3 <- cbind(CASE=tmp1, CONTROL=tmp2[match(names(tmp1), names(tmp2))])
dqa1.4.n <- fisher.test(cbind(colSums(tmp3)-tmp3[4, ], tmp3[4, ])) # 0401
tmp1 <- freqTable(pha.types.4[6, match(tin.uveitis.no, colnames(pha.types.2))])
tmp2 <- freqTable(pha.types.4[6, grepl('DT', colnames(pha.types.4))])
tmp3 <- cbind(CASE=tmp1, CONTROL=tmp2[match(names(tmp1), names(tmp2))])
dqa1.4.fq.n <- tmp3[4, ]
# DQB1
tmp1 <- table(pha.types.4[7, match(tin.uveitis.no, colnames(pha.types.2))])
tmp2 <- table(pha.types.4[7, grepl('DT', colnames(pha.types.4))])
tmp3 <- cbind(CASE=tmp1, CONTROL=tmp2[match(names(tmp1), names(tmp2))])
dqb1.4.n <- fisher.test(cbind(colSums(tmp3)-tmp3[4, ], tmp3[4, ])) # 0402
tmp1 <- freqTable(pha.types.4[7, match(tin.uveitis.no, colnames(pha.types.2))])
tmp2 <- freqTable(pha.types.4[7, grepl('DT', colnames(pha.types.4))])
tmp3 <- cbind(CASE=tmp1, CONTROL=tmp2[match(names(tmp1), names(tmp2))])
dqb1.4.fq.n <- tmp3[4, ]
# DRB1
tmp1 <- table(pha.types.4[8, match(tin.uveitis.no, colnames(pha.types.2))])
tmp2 <- table(pha.types.4[8, grepl('DT', colnames(pha.types.4))])
tmp3 <- cbind(CASE=tmp1, CONTROL=tmp2[match(names(tmp1), names(tmp2))])
drb1.4.n <- fisher.test(cbind(colSums(tmp3)-tmp3[3, ], tmp3[3, ])) # 0802
tmp1 <- freqTable(pha.types.4[8, match(tin.uveitis.no, colnames(pha.types.2))])
tmp2 <- freqTable(pha.types.4[8, grepl('DT', colnames(pha.types.4))])
tmp3 <- cbind(CASE=tmp1, CONTROL=tmp2[match(names(tmp1), names(tmp2))])
drb1.4.fq.n <- tmp3[3, ]


# 4-digit haplotypes
tmp1 <- table(apply(pha.types.4[5:8, grepl('TIN', colnames(pha.types.4))], 2, paste, collapse='-'))
tmp2 <- table(apply(pha.types.4[5:8, grepl('TIN', colnames(pha.types.4))==F], 2, paste, collapse='-'))
tmp3 <- cbind(CASE=tmp1, CONTROL=tmp2[match(names(tmp1), names(tmp2))])
tmp3 <- tmp3[complete.cases(tmp3), ]
hapl.4 <- fisher.test(cbind(colSums(tmp3)-tmp3[11, ], tmp3[11, ])) # 0301-0401-0402-0802
tmp1 <- freqTable(apply(pha.types.4[5:8, grepl('TIN', colnames(pha.types.4))], 2, paste, collapse='-'))
tmp2 <- freqTable(apply(pha.types.4[5:8, grepl('TIN', colnames(pha.types.4))==F], 2, paste, collapse='-'))
tmp3 <- cbind(CASE=tmp1, CONTROL=tmp2[match(names(tmp1), names(tmp2))])
tmp3 <- tmp3[complete.cases(tmp3), ]
hapl.4.fq <- tmp3[11, ]

# 4-digit haplotypes for TINU only
tmp1 <- table(apply(pha.types.4[5:8, match(tin.uveitis, colnames(pha.types.2))], 2, paste, collapse='-'))
tmp2 <- table(apply(pha.types.4[5:8, grepl('DT', colnames(pha.types.4))], 2, paste, collapse='-'))
tmp3 <- cbind(CASE=tmp1, CONTROL=tmp2[match(names(tmp1), names(tmp2))])
tmp3 <- tmp3[complete.cases(tmp3), ]
hapl.4.u <- fisher.test(cbind(colSums(tmp3)-tmp3[10, ], tmp3[10, ])) # 0301-0401-0402-0802
tmp1 <- freqTable(apply(pha.types.4[5:8, match(tin.uveitis, colnames(pha.types.2))], 2, paste, collapse='-'))
tmp2 <- freqTable(apply(pha.types.4[5:8, grepl('DT', colnames(pha.types.4))], 2, paste, collapse='-'))
tmp3 <- cbind(CASE=tmp1, CONTROL=tmp2[match(names(tmp1), names(tmp2))])
tmp3 <- tmp3[complete.cases(tmp3), ]
hapl.4.fq.u <- tmp3[10, ]

# 4-digit haplotypes for TIN only
tmp1 <- table(apply(pha.types.4[5:8, match(tin.uveitis.no, colnames(pha.types.2))], 2, paste, collapse='-'))
tmp2 <- table(apply(pha.types.4[5:8, grepl('DT', colnames(pha.types.4))], 2, paste, collapse='-'))
tmp3 <- cbind(CASE=tmp1, CONTROL=tmp2[match(names(tmp1), names(tmp2))])
tmp3 <- tmp3[complete.cases(tmp3), ]
hapl.4.n <- fisher.test(cbind(colSums(tmp3)-tmp3[4, ], tmp3[4, ])) # 0301-0401-0402-0802
tmp1 <- freqTable(apply(pha.types.4[5:8, match(tin.uveitis.no, colnames(pha.types.2))], 2, paste, collapse='-'))
tmp2 <- freqTable(apply(pha.types.4[5:8, grepl('DT', colnames(pha.types.4))], 2, paste, collapse='-'))
tmp3 <- cbind(CASE=tmp1, CONTROL=tmp2[match(names(tmp1), names(tmp2))])
tmp3 <- tmp3[complete.cases(tmp3), ]
hapl.4.fq.n <- tmp3[4, ]


# collect results with freqs & write 
hla.tin.fisher.fq <- rbind(
  DPB1_0301=c(dpb1.4.fq, dpb1.4$estimate, dpb1.4$conf.int, dpb1.4$p.value),
  DQA1_0401=c(dqa1.4.fq, dqa1.4$estimate, dqa1.4$conf.int, dqa1.4$p.value),
  DQB1_0402=c(dqb1.4.fq, dqb1.4$estimate, dqb1.4$conf.int, dqb1.4$p.value),
  DRB1_0802=c(drb1.4.fq, drb1.4$estimate, drb1.4$conf.int, drb1.4$p.value),
  HAPLOTYPE_4=c(hapl.4.fq, hapl.4$estimate, hapl.4$conf.int, hapl.4$p.value)
)
colnames(hla.tin.fisher.fq) <- c('Freq. Cases', 'Freq. Controls', 'OR', 'CI 2.5%', 'CI 97.5%', 'Fisher test p-value')
write.table(hla.tin.fisher.fq, './data/HLA_4digFQ_TIN_Fisher', quote=F, sep='\t')



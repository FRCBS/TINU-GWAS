library(data.table)

tmp <- fread('./data/IC1-IC3_imputed.merged.fam', data.table=F)

tmp <- tmp[grepl('SPR.D|TIN', tmp[, 1]), ]

write.table(tmp[, 1:2], './data/samples.list', quote=F, sep='\t', row.names=F, col.names=F)


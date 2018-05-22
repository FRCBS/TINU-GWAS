
library(data.table)

ic3.map <- fread('./data/IC3F_chr1_22_all.impute2_filtered0.5.map', data.table=F)
ic1.map <- fread('./data/IC1_chr1_22_all.impute2_filtered0.5.map', data.table=F)
ic1.ic3.map <- intersect(ic1.map[, 2], ic3.map[, 2])

write.table(ic1.ic3.map, './data/common.imputed.snps', quote=F, sep='\t', col.names=F, row.names=F)




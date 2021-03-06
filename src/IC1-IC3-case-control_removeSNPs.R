# Cohort association test results

library(data.table)

cc.assoc <- fread('./data/IC1-IC3_case-control.assoc', data.table=F)
cc.assoc.signif <- na.omit(cc.assoc[cc.assoc[, 'P'] < 0.001, 'SNP'])
cc.assoc.signif <- colnames(gen)[which(gsub('_[0-9]$', '', colnames(gen)) %in% cc.assoc.signif)]
write.table(cc.assoc.signif, './data/IC1-IC3_case-control.snps', quote=F, row.names=F, col.names=F)


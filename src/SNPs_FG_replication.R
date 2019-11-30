
library(data.table)
library(tidyverse)
library(ggpubr)

## Chronic iridocyclitis, ICD-10 code H2O.1

# iridochronic replication results rs192601611
r4.dat <- fread('./data/FinnGen_R4_iridochronic_SNPs.tsv', data.table=F)
r4.dat <- mutate(r4.dat, stdplus=q50+stderr, stdminus=q50-stderr)

# 90% CI for sampling distribution of allele freqs
cairo_pdf('./results/R4_iridochronic_rs192601611.pdf', height=4, width=4)
ggplot(r4.dat, aes(age, q50, ymin=stdminus, ymax=stdplus)) +
  geom_pointrange() +
  ylab('minor allele frequency') + xlab('age group') +
  facet_wrap(~group) +
  theme(panel.background=element_blank(), axis.line=element_line(colour="black", size=.2)) +
  theme(strip.text=element_text(colour='grey20', size=11), 
        axis.text.x=element_text(size=10)) 
dev.off()


library(data.table)
library(biomaRt)
library(tidyverse)
library(ggpubr)


## Plotting 

# plots of p-values by position

# data; TIN+TINU associations
spa.tin.p.range      <- fread('./results/tin.range.best.snps', data.table=F)
spa.tin.p.range$Pval[is.na(spa.tin.p.range$Pval)] <- 1
spa.tin.p.range$CHR  <- paste('Chr', spa.tin.p.range$CHR, sep='')
spa.tin.p.range$POS  <- spa.tin.p.range$POS/1e6
# data; TINU associations
spa.tinu.p.range     <- fread('./results/tinu.range.best.snps', data.table=F)
spa.tinu.p.range$Pval[is.na(spa.tinu.p.range$Pval)] <- 1
spa.tinu.p.range     <- spa.tinu.p.range %>% filter(CHR=='9')
spa.tinu.p.range$CHR <- paste('Chr', spa.tinu.p.range$CHR, sep='')
spa.tinu.p.range$POS <- spa.tinu.p.range$POS/1e6

p1 <- ggplot(spa.tin.p.range[6000:30000, ], aes(x=POS, y=-log10(Pval))) +
  geom_point(size=0.01, alpha=0.8, color='blue3') +
  geom_hline(yintercept=7.30103, color='red', alpha=0.55, size=0.4) +
  annotate("rect", xmin=33.01, xmax=33.2, ymin=0, ymax=8.1, alpha=.1) +
  annotate("rect", xmin=32.55, xmax=32.7, ymin=0, ymax=8.1, alpha=.1) +
  annotate("text", x=33.01+2e-1, y=8.3, label='HLA-DPA1/DPB1', colour='grey20', size=3) +
  annotate("text", x=32.55-1e-1, y=8.3, label='HLA-DRB1', colour='grey20', size=3) +
  geom_point(data=spa.tin.p.range[6000:30000, ] %>% filter(Pval<5e-8), aes(x=POS, y=-log10(Pval) %>% data.frame), 
             color='red', shape=5, fill='red', alpha=0.8) +
  facet_grid(. ~ CHR) +
  labs(y=expression("-log"[10]*"(p-value)")) +
  xlab('Position (Mb)')

p2 <- ggplot(spa.tinu.p.range[4000:15500, ], aes(x=POS, y=-log10(Pval))) +
  geom_point(size=0.01, alpha=0.8, color='blue3') +
  geom_hline(yintercept=7.30103, color='red', alpha=0.55, size=0.4) +
  annotate("rect", xmin=27.948, xmax=28.670, ymin=0, ymax=8.1, alpha=.1) +
  annotate("text", x=28.30, y=8.3, label='LINGO2', colour='grey20', size=3.5) +
  geom_point(data=spa.tinu.p.range[4000:15500, ] %>% filter(Pval<5e-8), aes(x=POS, y=-log10(Pval) %>% data.frame), 
             color='red', shape=5, fill='red', alpha=0.8) +
  facet_grid(. ~ CHR) +
  labs(y=expression("-log"[10]*"(p-value)")) +
  xlab('Position (Mb)')

pdf(width=9, height=4, file='./results/signif_pvals.pdf')
ggpubr::ggarrange(p1, p2, ncol=2, nrow=1, labels=c('a', 'b'), align='h', 
                  font.label=list(size=16))
dev.off()


# TINU and GTEx whole blood eqtl
lingo2.gtex           <- fread('./data/Lingo2_GTEx Portal.csv', data.table=F)
colnames(lingo2.gtex) <- colnames(lingo2.gtex) %>% map_chr(., gsub, pattern=' |-', replacement='_')
lingo2.gtex           <- lingo2.gtex %>% mutate(., POS=str_split_fixed(Variant_Id, '_', 3)[, 2] %>% as.numeric)
lingo2.gtex           <- lingo2.gtex %>% filter(Tissue=='Whole Blood')
lingo2.gtex$POS       <- lingo2.gtex$POS/1e6 
lingo2.gtex$P_Value   <- as.numeric(lingo2.gtex$P_Value)
lingo2.tinu           <- full_join(lingo2.gtex, spa.tinu.p.range[8000:9600, ], by='POS')
lingo2.tinu           <- lingo2.tinu %>% gather(., Source, Value, Pval, P_Value)
lingo2.tinu$Value     <- lingo2.tinu$Value %>% -log10(.) 
lingo2.tinu$Source    <- c('TINU association', 'Whole Blood EQTL')[lingo2.tinu$Source %>% factor]

pdf(width=5, height=4, file='./results/eqtl_pvals.pdf')
ggplot(lingo2.tinu, aes(x=POS, y=Value)) +
  geom_line() +
  geom_point() +
  coord_cartesian(xlim=c(28.5, 29.5), ylim=c(4.8, 8)) +
  facet_grid(Source ~., switch='both') +
  xlab('Position (Mb)') + 
  labs(y=expression("-log"[10]*"(p-value)")) 
dev.off()



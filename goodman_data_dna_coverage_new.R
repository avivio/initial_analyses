load.packages()

fitesq.location <- 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\fitseq_data_with_meta_day_12_no_umi_05_08.csv'

fitseq_data_residuals <- load.fitseq.data(fitesq.location)

goodman.data.dna <- fitseq_data_residuals %>% filter(day == 12, lineage == 'A') %>% select(Name, new.name, Promoter.Display, RBS.Display, Count.DNA, Prot, frequency)
                                    
goodman.data.dna.no.wt <-goodman.data.dna %>% filter(RBS.Display != 'WT') 
goodman.data.dna.no.wt <- collect(goodman.data.dna.no.wt)
goodman.data.dna.no.wt$RBS.Display <- factor(goodman.data.dna.no.wt$RBS.Display,
                                             levels = c("Strong", "Mid", "Weak"))
goodman.data.dna.no.wt$Promoter.Display <- factor(goodman.data.dna.no.wt$Promoter.Display,
                                             levels = c("High", "Low"))

png('C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\results\\for_thesis\\goodman_dna_vs_promoter_rbs.png',
    type="cairo",    units="in", width=10, height=6, pointsize=12, res=500)    
ggplot(goodman.data.dna.no.wt, aes(x =interaction(RBS.Display,Promoter.Display), y = log2(Count.DNA+1) )) +
  geom_boxplot(fill = hue_pal()(2)[2]) +
  ylab('Log 2 DNA coverage') +
  xlab('Promoter and RBS type') +
  scale_x_discrete(labels=c("High Strong", "High Mid", "High Weak","Low Strong", "Low Mid", "Low Weak")) +
  theme_aviv
dev.off()

goodman.data.dna.no.wt.1 <- goodman.data.dna.no.wt %>% filter(Promoter.Display == 'Low',RBS.Display == 'Mid')
goodman.data.dna.no.wt.2 <- goodman.data.dna.no.wt %>% filter(Promoter.Display == 'Low',RBS.Display == 'Weak')

wilcox.test(goodman.data.dna.no.wt.1$Count.DNA,goodman.data.dna.no.wt.2$Count.DNA,alternative = 'less')

# protein level by promoter rbs
png('C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\results\\for_thesis\\goodman_prot_vs_promoter_rbs.png',
    type="cairo",    units="in", width=10, height=6, pointsize=12, res=500) 
ggplot(goodman.data.dna.no.wt, aes(x =interaction(RBS.Display,Promoter.Display), y = log2(Prot + 1) )) +
  geom_boxplot(fill = hue_pal()(2)[2]) +
  ylab('Log 2 expression level') +
  xlab('Promoter and RBS type') +
  scale_x_discrete(labels=c("High Strong", "High Mid", "High Weak","Low Strong", "Low Mid", "Low Weak")) +
  theme_aviv
dev.off()
# differences in dna coverage, fill
png('C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\results\\for_thesis\\dna_coverage_rbs_fill.png',
    type="cairo",    units="in", width=8, height=8, pointsize=12, res=500)    
ggplot(goodman.data.dna.no.wt, aes(log10(Count.DNA), fill = RBS.Display)) +
  stat_bin(geom='area',position = 'identity',alpha=.5) +
  theme_aviv + 
  ylab('# Designs') +
  xlab('Log 10 DNA coverage') +
  guides(fill = guide_legend("RBS"))
dev.off()


# differences in dna coverage, fill
png('C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\results\\for_thesis\\dna_coverage_rbs_facet.png',
    type="cairo",    units="in", width=8, height=8, pointsize=12, res=500)    
ggplot(goodman.data.dna.no.wt, aes(log10(Count.DNA))) +
  stat_bin(geom='area',position = 'identity', fill = hue_pal()(2)[2]) +
  theme_aviv + 
  ylab('# Designs') +
  xlab('Log 10 DNA coverage') +
  facet_grid(RBS.Display~.)
dev.off()
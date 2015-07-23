load.packages()

fitseq_data_residuals <- load.db()

goodman.data.dna <- fitseq_data_residuals %>% filter(day == 12, lineage == 'A') %>% select(Name, new_name, Promoter_Display, RBS_Display, Count_DNA)
                                    
goodman.data.dna.no.wt <-goodman.data.dna %>% filter(RBS_Display != 'WT') 
goodman.data.dna.no.wt <- collect(goodman.data.dna.no.wt)
goodman.data.dna.no.wt$RBS_Display <- factor(goodman.data.dna.no.wt$RBS_Display,
                                             levels = c("Strong", "Mid", "Weak"))
goodman.data.dna.no.wt$Promoter_Display <- factor(goodman.data.dna.no.wt$Promoter_Display,
                                             levels = c("High", "Low"))


ggplot(goodman.data.dna.no.wt, aes(x =interaction(RBS_Display,Promoter_Display), y = log2(Count_DNA+1) )) +
  geom_boxplot() +
  ggtitle('DNA Coverage in by promoter and RBS type') +
  ylab('Log 2 DNA coverage') +
  xlab('Promoter and RBS type') +
  scale_x_discrete(labels=c("High Strong", "High Mid", "High Weak","Low Strong", "Low Mid", "Low Weak")) +
  theme_aviv


goodman.data.dna.no.wt.1 <- goodman.data.dna.no.wt %>% filter(Promoter_Display == 'Low',RBS_Display == 'Mid')
goodman.data.dna.no.wt.2 <- goodman.data.dna.no.wt %>% filter(Promoter_Display == 'Low',RBS_Display == 'Weak')

wilcox.test(goodman.data.dna.no.wt.1$Count_DNA,goodman.data.dna.no.wt.2$Count_DNA,alternative = 'less')

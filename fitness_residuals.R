require(dplyr)
require(tidyr)
require(ggplot2)
theme_aviv <-     theme_minimal() + 
  theme( axis.line = element_line(colour = "black"),
         axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold"),
         strip.text.x = element_text(size=16,face="bold"), strip.text.y = element_text(size=16,face="bold"), 
         plot.title = element_text(size = 25,face = "bold"))

fitseq.data.location  <- 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\clean_data\\goodman_salis_tuller_fitseq_2_mismatch_with_0.csv'
fitseq.data = read.csv(fitseq.data.location)
fitseq.data = transform(fitseq.data, Prot = as.numeric(as.character(Prot))
                        ,Bin.Pct.1 = as.numeric(as.character((Bin.Pct.1))),
                        Count.RNA = as.numeric(as.character((Count.RNA))))
fitseq.data.tidy <- fitseq.data %>% select(-Count.A.DNA,-Count.A.RNA,-Count.B.DNA,-Count.B.RNA,
                                           -Bin.1,-Bin.2,-Bin.3,-Bin.4,-Bin.5,-Bin.6,-Bin.7,
                                           -Bin.8,-Bin.9,-Bin.10,-Bin.11,-Bin.12,-RNA.A,-RNA.B,
                                           -Bin.Pct.1,-Bin.Pct.2,-Bin.Pct.3,-Bin.Pct.4,-Bin.Pct.5,
                                           -Bin.Pct.6,-Bin.Pct.7,-Bin.Pct.8,-Bin.Pct.9,-Bin.Pct.10,
                                           -Bin.Pct.11,-Bin.Pct.12,-Insuff.Prot,-Insuff.DNA,
                                           -Insuff.RNA,-Fltr.BelowRange,-Fltr.AboveRange ,
                                           -Fltr.SetGood,-full.seq)
fitseq.data.tidy  <- fitseq.data.tidy %>%   gather(sample, frequency, A_24:F_0)
fitseq.data.tidy  <- fitseq.data.tidy %>% separate(sample, c("lineage", "day"), '_',convert = T)
fitseq.data.tidy$RBS.Display <- factor(fitseq.data.tidy$RBS.Display,
                                       levels = c("Strong", "Mid", "Weak", "WT"))




fitseq.data.tidy <- fitseq.data.tidy %>% 
  group_by(day,lineage) %>% 
  mutate(sum.sample= sum(frequency),norm.freq = frequency/sum.sample,
         sum.anc= sum(anc_1),norm.anc = anc_1/sum.anc,freq.norm.anc = norm.freq/norm.anc,
         sum.sample.1= sum(frequency+1),norm.freq.1 = (frequency+1)/sum.sample.1,
         sum.anc.1= sum(anc_1+1),norm.anc.1 = (anc_1+1)/sum.anc.1,freq.norm.anc.1 = norm.freq.1/norm.anc.1)

day.number = 12
lineage.letter = 'C'
fitseq.current.sample <-  fitseq.data.tidy %>% filter(day == day.number, lineage == lineage.letter)




protein.fitness.lm <-  lm(freq.norm.anc.1 ~ Prot, data=fitseq.current.sample)
protein.fitness.coeffs <- coefficients(protein.fitness.lm)

fitseq.data.prediction <- fitseq.current.sample %>% 
  mutate(prediction =protein.fitness.coeffs[1] + protein.fitness.coeffs[2]*Prot,
         fit.resid = freq.norm.anc.1 -  prediction)

fitseq.data.prediction <- fitseq.data.prediction %>% filter(fit.resid < 2) 

fitseq.data.prediction.no.wt <- fitseq.data.prediction %>% 
  filter(RBS.Display!='WT')

ggplot(fitseq.data.prediction, aes(fit.resid)) + 
  geom_histogram() + 
  theme_aviv +
  xlab('Fitness residuals') +
  ggtitle(paste('Fitness residual histogram \nfor lineage', lineage.letter,'day', day.number))


ggplot(fitseq.data.prediction.no.wt, aes(fit.resid)) + 
  geom_histogram() + 
  theme_aviv +
  xlab('Fitness residuals') +
  ggtitle(paste('Fitness residual histogram \nfor lineage', lineage.letter,'day', day.number)) +
  facet_grid(Promoter.Display~RBS.Display)


ggplot(fitseq.data.prediction, aes(y = log2(fit.resid),x = log2(salis.init))) + 
  geom_point() + 
  theme_aviv +
  ylab('Log 2 fitness residuals') +
  xlab('GC content') +
  ggtitle(paste('Log 2 Fitness residuals vs GC content \nfor lineage', lineage.letter,'day', day.number)) 



ggplot(fitseq.data.prediction.no.wt, aes(y = log2(fit.resid),x = TASEP.avgRiboNum)) + 
  geom_point() + 
  theme_aviv +
  ylab('Log 2 fitness residuals') +
  xlab('GC content') +
  ggtitle(paste('Log 2 Fitness residuals vs GC content \nfor lineage', lineage.letter,'day', day.number)) +
  facet_grid(Promoter.Display~RBS.Display)

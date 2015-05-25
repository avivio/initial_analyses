#clean data

require(tidyr)
require(ggplot2)
require(dplyr)
require(scales)
require(plotrix)

fitseq.data.location  <- 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\clean_data\\goodman_salis_tuller_fitseq_2_mismatch_with_0.csv'
summary.data.location  <- 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\clean_data\\final_summary_14-05-15-1551.csv'
fitseq.data = read.csv(fitseq.data.location)
fitseq.data.tidy <- fitseq.data %>% select(-Count.A.DNA,-Count.A.RNA,-Count.B.DNA,-Count.B.RNA,
                                                -Bin.1,-Bin.2,-Bin.3,-Bin.4,-Bin.5,-Bin.6,-Bin.7,
                                                -Bin.8,-Bin.9,-Bin.10,-Bin.11,-Bin.12,-RNA.A,-RNA.B,
                                                -Bin.Pct.1,-Bin.Pct.2,-Bin.Pct.3,-Bin.Pct.4,-Bin.Pct.5,
                                                -Bin.Pct.6,-Bin.Pct.7,-Bin.Pct.8,-Bin.Pct.9,-Bin.Pct.10,
                                                -Bin.Pct.11,-Bin.Pct.12,-Insuff.Prot,-Insuff.DNA,
                                                -Insuff.RNA,-Fltr.BelowRange,-Fltr.AboveRange ,
                                                -Fltr.SetGood)
fitseq.data.tidy  <- fitseq.data.tidy %>%   gather(sample, frequency, A_24:F_0)
fitseq.data.tidy  <- fitseq.data.tidy %>% separate(sample, c("lineage", "day"), '_',convert = T)

#mean of reads over total reads per group over time

mean.freq.summary  <- 
  fitseq.data.tidy %>%
  group_by(lineage,day) %>% 
  mutate(sum.freq = sum(frequency),norm.freq = frequency/sum.freq) %>%
  group_by(day,lineage,Promoter.Display,RBS.Display) %>%
  summarise(mean.freq = mean(norm.freq),se = std.error(norm.freq),sd = sd(norm.freq)) %>%
  filter(lineage != 'anc',RBS.Display != 'WT')




ggplot(
  mean.freq.summary,
  aes(x=day,y=mean.freq,
      color=interaction(Promoter.Display,RBS.Display),shape=interaction(Promoter.Display,RBS.Display))) +
  geom_line() + 
  facet_wrap(~lineage) +
  geom_point(size = 5) +
  geom_errorbar(aes(ymax = mean.freq + se, ymin=mean.freq - se)) +
  theme_minimal() + 
  theme( axis.line = element_line(colour = "black"),
                     axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold"),
                     legend.text=element_text(size=14), legend.title=element_text(size=16,face="bold"),
                     plot.title = element_text(size = 25,face = "bold")) +
  guides(colour = guide_legend("Promoters and RBSs"),shape = guide_legend("Promoters and RBSs")) +
  ylab('Mean of frequency normalized to total frequency\n') + ggtitle('Mean Frequency Over Days Per Lineage error bars standard error') +
  coord_cartesian(ylim=c(0,max(mean.freq.summary$se+mean.freq.summary$mean.freq)))
# coord_cartesian(ylim=c(0,3e-04))

#median of reads over total reads per group over time


median.freq.summary  <- 
  fitseq.data.tidy %>%
  group_by(lineage,day) %>% 
  
  mutate(sum.freq = sum(frequency),norm.freq = frequency/sum.freq) %>%
  group_by(day,lineage,Promoter.Display,RBS.Display) %>%
  summarise(median.freq = median(norm.freq),sd = sd(norm.freq)) %>%
  filter(lineage != 'anc',RBS.Display != 'WT')

ggplot(
  median.freq.summary,
  aes(x=day,y=median.freq,
      color=interaction(Promoter.Display,RBS.Display),shape=interaction(Promoter.Display,RBS.Display))) +
  geom_line() + 
  geom_errorbar(aes(ymax = median.freq + sd, ymin=median.freq - sd)) +
  facet_wrap(~lineage) +
  geom_point(size = 5) +
  theme_minimal() + 
  theme( axis.line = element_line(colour = "black"),
         axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold"),
         legend.text=element_text(size=14), legend.title=element_text(size=16,face="bold"),
         plot.title = element_text(size = 25,face = "bold")) +
  guides(colour = guide_legend("Promoters and RBSs"),shape = guide_legend("Promoters and RBSs")) +
  ylab('Median of frequency normalized to total frequency\n') + ggtitle('Median Frequency Over Days Per Lineage error bars standard deviation') +
#   coord_cartesian(ylim=c(0,max(median.freq.summary$sd+median.freq.summary$median.freq)))
coord_cartesian(ylim=c(0,2e-04))



#mean of reads over total reads per group normalized to ancestor for each read over time

mean.freq.norm.anc.per.design.summary  <- 
  fitseq.data.tidy %>%
  group_by(lineage,day) %>% 
  
  mutate(sum.freq = sum(frequency),norm.freq.total = frequency/sum.freq,norm.anc = anc_1/3641463,
         freq.norm.anc = norm.freq.total/norm.anc) %>% 
  filter(!is.na(freq.norm.anc),freq.norm.anc!=Inf) %>%
  group_by(day,lineage,Promoter.Display,RBS.Display) %>%
  summarise(mean.freq = mean(freq.norm.anc),se = std.error(freq.norm.anc), sd = sd(freq.norm.anc),
            median.freq = median(freq.norm.anc)) %>%
  filter(lineage != 'anc',RBS.Display != 'WT')




ggplot(
  mean.freq.norm.anc.per.design.summary,
  aes(x=day,y=mean.freq,
      color=interaction(Promoter.Display,RBS.Display),shape=interaction(Promoter.Display,RBS.Display))) +
  geom_line() + 
  facet_wrap(~lineage) +
  geom_point(size = 5) +
  geom_errorbar(aes(ymax = mean.freq + sd, ymin=mean.freq - sd)) +
  theme_minimal() + 
  theme( axis.line = element_line(colour = "black"),
         axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold"),
         legend.text=element_text(size=14), legend.title=element_text(size=16,face="bold"),
         plot.title = element_text(size = 25,face = "bold")) +
  guides(colour = guide_legend("Promoters and RBSs"),shape = guide_legend("Promoters and RBSs")) +
  ylab('Mean change in normalized frequency compared to ancestor\n') +
  ggtitle('Mean change in normalized frequency compared to ancestor over time')  +
#    coord_cartesian(ylim=c(0,max(mean.freq.norm.anc.per.design.summary$sd+mean.freq.norm.anc.per.design.summary$mean.freq)))
    coord_cartesian(ylim=c(0,3))


#median of reads over total reads per group normalized to ancestor for each read over time

ggplot(
  mean.freq.norm.anc.per.design.summary,
  aes(x=day,y=median.freq,
      color=interaction(Promoter.Display,RBS.Display),shape=interaction(Promoter.Display,RBS.Display))) +
  geom_line() + 
  facet_wrap(~lineage) +
  geom_point(size = 5) +
  geom_errorbar(aes(ymax = median.freq + sd, ymin=median.freq - sd)) +
  theme_minimal() + 
  theme( axis.line = element_line(colour = "black"),
         axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold"),
         legend.text=element_text(size=14), legend.title=element_text(size=16,face="bold"),
         plot.title = element_text(size = 25,face = "bold")) +
  guides(colour = guide_legend("Promoters and RBSs"),shape = guide_legend("Promoters and RBSs")) +
  ylab('Median change in normalized frequency compared to ancestor\n') +
  ggtitle('Median  change in normalized frequency compared to ancestor over time')  +
#      coord_cartesian(ylim=c(0,max(mean.freq.norm.anc.per.design.summary$sd+mean.freq.norm.anc.per.design.summary$median.freq)))
  coord_cartesian(ylim=c(0,2))





#mean of reads over total reads per group normalized to ancestor for each read over time

mean.freq.over.anc  <- 
  fitseq.data.tidy %>%
  group_by(lineage,day) %>% 
  
  mutate(sum.freq = sum(frequency),norm.freq.total = frequency/sum.freq,norm.anc = anc_1/3641463,
         freq.norm.anc = norm.freq.total/norm.anc) %>% 
  filter(!is.na(freq.norm.anc),freq.norm.anc!=Inf) %>%
  group_by(day,lineage,Promoter.Display,RBS.Display) %>%
  summarise(mean.norm.freq.total = mean(norm.freq.total),
            mean.norm.anc = mean(norm.anc),mean.freq.over.anc = mean.norm.freq.total/mean.norm.anc,
            median.norm.freq.total = median(norm.freq.total),
            median.norm.anc = median(norm.anc),median.freq.over.anc = median.norm.freq.total/median.norm.anc) %>%
  filter(lineage != 'anc',RBS.Display != 'WT')




ggplot(
  mean.freq.over.anc,
  aes(x=day,y=mean.freq.over.anc,
      color=interaction(Promoter.Display,RBS.Display),shape=interaction(Promoter.Display,RBS.Display))) +
  geom_line() + 
  facet_wrap(~lineage) +
  geom_point(size = 5) +
#   geom_errorbar(aes(ymax = mean.freq + se, ymin=mean.freq - se)) +
  theme_minimal() + 
  theme( axis.line = element_line(colour = "black"),
         axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold"),
         legend.text=element_text(size=14), legend.title=element_text(size=16,face="bold"),
         plot.title = element_text(size = 25,face = "bold")) +
  guides(colour = guide_legend("Promoters and RBSs"),shape = guide_legend("Promoters and RBSs")) +
  ylab('Mean frequency over mean ancestor frequency') + ggtitle('Mean frequency over mean ancestor frequency over time')

#median

ggplot(
  mean.freq.over.anc,
  aes(x=day,y=median.freq.over.anc,
      color=interaction(Promoter.Display,RBS.Display),shape=interaction(Promoter.Display,RBS.Display))) +
  geom_line() + 
  facet_wrap(~lineage) +
  geom_point(size = 5) +
  #   geom_errorbar(aes(ymax = mean.freq + se, ymin=mean.freq - se)) +
  theme_minimal() + 
  theme( axis.line = element_line(colour = "black"),
         axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold"),
         legend.text=element_text(size=14), legend.title=element_text(size=16,face="bold"),
         plot.title = element_text(size = 25,face = "bold")) +
  guides(colour = guide_legend("Promoters and RBSs"),shape = guide_legend("Promoters and RBSs")) +
  ylab('Median frequency over median ancestor frequency') + ggtitle('Median frequency over median ancestor frequency over time')

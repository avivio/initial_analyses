#clean data

require(tidyr)
require(ggplot2)
require(dplyr)


result.dir = 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\results\\no_umi_rerun\\frequency_histograms_days\\'

fitseq.data.location  <- 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\fitseq_data_with_meta_no_umi_generations.csv'
summary.data.location  <- 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\clean_data\\final_summary_14-05-15-1551.csv'
fitseq.data = read.csv(fitseq.data.location)
fitseq.data.tidy <- fitseq.data %>% select(-Count.A.DNA,-Count.A.RNA,-Count.B.DNA,-Count.B.RNA,
                                                -Bin.1,-Bin.2,-Bin.3,-Bin.4,-Bin.5,-Bin.6,-Bin.7,
                                                -Bin.8,-Bin.9,-Bin.10,-Bin.11,-Bin.12,-RNA.A,-RNA.B,
                                                -Bin.Pct.1,-Bin.Pct.2,-Bin.Pct.3,-Bin.Pct.4,-Bin.Pct.5,
                                                -Bin.Pct.6,-Bin.Pct.7,-Bin.Pct.8,-Bin.Pct.9,-Bin.Pct.10,
                                                -Bin.Pct.11,-Bin.Pct.12,-Insuff.Prot,-Insuff.DNA,
                                                -Insuff.RNA,-Fltr.BelowRange,-Fltr.AboveRange ,
                                                -Fltr.SetGood,-full.seq)
fitseq.data.tidy  <- fitseq.data.tidy %>%   gather(sample, frequency, A_168:F_0)
fitseq.data.tidy  <- fitseq.data.tidy %>% separate(sample, c("lineage", "generation"), '_',convert = T)
fitseq.data.tidy$RBS.Display <- factor(fitseq.data.tidy$RBS.Display,
                              levels = c("Strong", "Mid", "Weak", "WT"))
fitseq.data.tidy$lineage <- factor(fitseq.data.tidy$lineage,
                                       levels = c("Ancestor", "A", "B", "C", "D", "E", "F"))


#historgrams of normalized frequency + 1 log transformed for different day, rbs-promoter groups, and lineages,
fitseq.data.tidy  <- fitseq.data.tidy %>% 
  group_by(generation,lineage) %>%
  mutate(sum.sample= sum(frequency),norm.freq = frequency/sum.sample,
         sum.anc= sum(anc_1),norm.anc = anc_1/sum.anc,freq.norm.anc = norm.freq/norm.anc,
         sum.sample.1= sum(frequency+1),norm.freq.1 = (frequency+1)/sum.sample.1,
         sum.anc.1= sum(anc_1+1),norm.anc.1 = (anc_1+1)/sum.anc.1,freq.norm.anc.1 = norm.freq.1/norm.anc.1)



line.day.histogram   <- function(lineage_letter,days,fitseq.data.tidy){
  days_text  <- paste(days,collapse = ', ')
  png(paste0(result.dir,'days_fill_normalized_frequency_histograms_lineage_',lineage_letter,'_hist.png'),
      type="cairo",    units="in", width=10, height=6, pointsize=12, res=500)    
 p =  ggplot(filter(fitseq.data.tidy,day %in% days,lineage == lineage_letter,RBS.Display != 'WT') ,
         aes(log2(freq.norm.anc.1),fill = factor(day)))  +
   stat_bin(geom='area',position = 'identity',alpha=.5) +
   theme_minimal() + 
    facet_grid(RBS.Display~Promoter.Display)  +
   theme_aviv +
   guides(color = guide_legend("Day"),fill = guide_legend("Day")) +
   xlab('\nLog 2 of normalized frequency + 1') +
   ylab('Count\n') +
   expand_limits(y=900)
 
  print(p)
  dev.off()
  
}


line.rbs.histogram   <- function(lineage_letter,generations,fitseq.data.tidy){
  generations_text  <- paste(generations,collapse = ', ')
  levels(fitseq.data.tidy$Promoter.Display) <- c("High promoter", "Low promoter")
  png(paste0(result.dir,'rbs_fill_frequency_histograms_lineage_',lineage_letter,'.png'),
      type="cairo",    units="in", width=9, height=11, pointsize=12, res=500)    
  p =  ggplot(filter(fitseq.data.tidy,generation %in% generations,lineage == lineage_letter,RBS.Display != 'WT') ,
              aes(log2(norm.freq.1),fill = factor(RBS.Display)))  +
    stat_bin(geom='area',position = 'identity',alpha=.5) +
    theme_minimal() + 
    facet_grid(generation~Promoter.Display)  +
    theme_aviv +
    guides(color = guide_legend("Generation"),fill = guide_legend("RBS")) +
    xlab('\nLog 2 frequency') +
    ylab('# designs\n') +
    expand_limits(y=1000) +
    coord_cartesian(xlim = c(-24,-10))
  
  print(p)
  dev.off()
  
}

for (lineage in c('A','B','C','D','E','F')){
  line.day.histogram(lineage,c(4,12,28),fitseq.data.tidy)
  
}

for (lineage in c('A','B','C','D','E','F')){
  line.rbs.histogram(lineage,c(28,56,84,112,140,168,196),fitseq.data.tidy)
  
}



ggplot(filter(fitseq.data.tidy,day %in% c(4,24,28),lineage == 'A',RBS.Display != 'WT') ,
       aes(log10(norm.freq),fill = factor(day),color = factor(day)))  +
  geom_density(alpha=0.5) +
  theme_minimal() + 
  facet_grid(RBS.Display~Promoter.Display)



ggplot(filter(fitseq.data.tidy,day %in% c(4,12,28),lineage == 'A',RBS.Display != 'WT') ,
       aes(log10(norm.freq),fill = factor(day),color = factor(day)))  +
  theme_minimal() + 
  facet_grid(RBS.Display~Promoter.Display)




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

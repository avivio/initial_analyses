require(dplyr)
require(tidyr)
require(ggplot2)


result.dir = 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\results\\'


pre.umi.data.location  <- 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\clean_data\\final_match_count_14-05-15-1551_2_mismatches.csv'
pre.umi.data = read.csv(pre.umi.data.location)
pre.umi.data.tidy  <- pre.umi.data %>%   gather(sample, frequency, A_24:F_4)
pre.umi.data.tidy  <- pre.umi.data.tidy %>% separate(sample, c("lineage", "day"), '_',convert = T)
pre.umi.data.tidy[pre.umi.data.tidy[,'lineage'] == 'anc','day'] <- 0

counts <- pre.umi.data.tidy %>% filter(frequency >10000) %>% group_by(day,lineage) %>% summarise(count = n())
all.day.lineage <- pre.umi.data.tidy %>% distinct(day,lineage) %>% select(day,lineage)
over.10k <- left_join(all.day.lineage,counts,,by = c('day','lineage'))
over.10k[is.na(over.10k)] <- 0                    

png(paste0(result.dir,'pre_umi_matches_above_10k.png'),units="in",  width=20, height=12, res=70)


ggplot(over.10k,aes(x = day, y = count, color = lineage, shape = lineage)) + 
  geom_point(size = 4) + geom_line(size = 1.2) + scale_shape_manual(values=c(0,15,16,1,17,2,25)) +
  theme_minimal() + 
  theme( axis.line = element_line(colour = "black"),
         axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold"),
         legend.text	 = element_text(size=14), legend.title = element_text(size=16,face="bold"), 
         plot.title = element_text(size = 25,face = "bold"),legend.position="none") +
  ggtitle('Pre umi match counts above 10,000 over time') +
  ylab('# Pre umi match counts above 10,000') +
  xlab('Day')


dev.off()
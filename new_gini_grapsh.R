require(dplyr)
require(tidyr)
require(ggplot2)
require(ineq)
require(scales)

theme_aviv <-     theme_minimal() + 
  theme( axis.line = element_line(colour = "black"),
         axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold"),
         strip.text.x = element_text(size=16,face="bold"), strip.text.y = element_text(size=16,face="bold"), 
         plot.title = element_text(size = 25,face = "bold"))

post.umi.data.location  <- 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\clean_data\\final_result_14-05-15-1551_2_mismatches.csv'
post.umi.data = read.csv(post.umi.data.location)
colnames(post.umi.data)[40] <- 'C_0'
post.umi.data  <- post.umi.data %>%   gather(sample, frequency, A_24:F_4)
post.umi.data  <- post.umi.data %>% separate(sample, c("lineage", "day"), '_',convert = T)
post.umi.data <- post.umi.data %>%
  group_by(day,lineage) %>% 
  mutate(sum.sample= sum(frequency),norm.freq = frequency/sum.sample)


pnas.correction.data.location  <- 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\clean_data\\pnas_correction_14-05-15-1551_2_mismatches.csv'
pnas.correction.data = read.csv(pnas.correction.data.location)
pnas.correction.data  <- pnas.correction.data %>%   gather(sample, frequency, A_24:F_4)
pnas.correction.data  <- pnas.correction.data %>% separate(sample, c("lineage", "day"), '_',convert = T)
pnas.correction.data <- pnas.correction.data %>%
  group_by(day,lineage) %>% 
  mutate(sum.sample= sum(frequency),norm.freq = frequency/sum.sample)


# joined_data <- left_join(post.umi.data,pnas.correction.data,by = c('X','lineage', 'day') )
# gini <- joined_data %>% 
#   group_by(day,lineage) %>%
#   summarise(post.umi.freq =  ineq(frequency.x,type="Gini"),pnas.correction.freq =  ineq(frequency.y,type="Gini"))
# gini <-  gini %>%
#   gather(data.set, gini,post.umi.freq:pnas.correction.freq)
# gini <-  gini %>%
#   filter(data.set == 'gini,post.umi.freq')
gini <- post.umi.data %>% 
    group_by(day,lineage) %>%
    summarise(gini =  ineq(frequency,type="Gini"))

ggplot(filter(gini,lineage =='B'), aes(y = gini,x = factor(day),fill = 'cyan')) + geom_bar(stat = 'identity',) +
   ylab('Gini index') + xlab('Day') + guides(fill = guide_legend(title = "Lineage")) + theme_aviv + 
  ggtitle('Gini index over time lineage B') + guides(fill=FALSE)
  

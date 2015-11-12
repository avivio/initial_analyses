
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

post.umi.data.location  <- 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\clean_data\\final_match_count_14-05-15-1551_2_mismatches.csv'
post.umi.data = read.csv(post.umi.data.location)
post.umi.data <- post.umi.data %>% select(-anc_2)
post.umi.data  <- post.umi.data %>%   gather(sample, frequency, A_24:F_4)
post.umi.data  <- post.umi.data %>% separate(sample, c("lineage", "generation"), '_',convert = T)
post.umi.data <- post.umi.data %>%
  group_by(generation,lineage) %>% 
  mutate(sum.sample= sum(frequency),norm.freq = frequency/sum.sample)


pre.umi.data.location  <- 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\match_count_2_mismatches_generations.csv'
pre.umi.data = read.csv(pre.umi.data.location)
pre.umi.data  <- pre.umi.data %>%   gather(sample, frequency, A_168:F_28)
pre.umi.data  <- pre.umi.data %>% separate(sample, c("lineage", "generation"), '_',convert = T)
pre.umi.data <- pre.umi.data %>%
  group_by(generation,lineage) %>% 
  mutate(sum.sample= sum(frequency),norm.freq = frequency/sum.sample)

post.umi.data$generation =pre.umi.data$generation


joined_data <- left_join(post.umi.data,pre.umi.data,by = c('X','lineage', 'generation') )
gini <- joined_data %>% 
  group_by(generation,lineage) %>%
  summarise(post.umi.freq =  ineq(frequency.x,type="Gini"),pre.umi.freq =  ineq(frequency.y,type="Gini")) 

# gini <- filter(gini,generation !=0)
# 
# gini[gini$generation == 1,'lineage'] <- 'A'
# gini[gini$generation == 1,'generation'] <- 0
# 
gini <- filter(gini,lineage !='anc')

gini <-  gini %>%
  gather(data.set, gini,post.umi.freq:pre.umi.freq)
# gini <-  gini %>%
#   filter(data.set == 'gini,post.umi.freq')
# gini <- post.umi.data %>% 
#     group_by(generation,lineage) %>%
#     summarise(gini =  ineq(frequency,type="Gini"))

png('C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\results\\for_thesis\\gini_over_time_lienages.png',
    type="cairo",    units="in", width=22, height=6, pointsize=12, res=500)    
ggplot(gini %>% filter(data.set == 'pre.umi.freq'), aes(y = gini,x = factor(generation))) + 
  geom_bar(stat = 'identity',position = 'dodge', fill = hue_pal()(2)[2]) +
   ylab('Gini index\n') + 
  xlab('\nGeneration') + 
#   scale_fill_discrete(name="Pre or Post\nUMI processing",
#                       breaks=c("pre.umi.freq", "post.umi.freq"),
#                       labels=c("Pre UMI", "Post UMI")) +
  # guides(fill = guide_legend(title = "Pre or Post\nUMI processing")) +
  theme_aviv +
  facet_wrap(~lineage,nrow = 1)
  
dev.off()
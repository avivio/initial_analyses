
rare  <- fitseq.data.tidy %>%
  select(Name,new.name, Promoter.Display,RBS.Display,Gene, CDS.type, day,lineage,frequency,Prot) %>%
  group_by(lineage,day) %>% 
  mutate(sum.freq = sum(frequency+1), normal.frequency=(1+frequency)/sum.freq)%>% 
  filter(CDS.type == 'Max Rare') 


common  <- fitseq.data.tidy %>%
  select(Name,new.name, Promoter.Display,RBS.Display,Gene, CDS.type, day,lineage,frequency,Prot) %>%
  group_by(lineage,day) %>% 
  mutate(sum.freq = sum(frequency+1), normal.frequency=(1+frequency)/sum.freq)%>%  
  filter(CDS.type == 'Min Rare')


joined  <- left_join(common,rare,by=c('Gene','RBS.Display','Promoter.Display','day','lineage'))

common.rare <- joined %>% select(Promoter = Promoter.Display, RBS = RBS.Display, Gene, day, lineage,
                                 rare.name = Name.y,rare.new.name = new.name.y,rare.frequency = frequency.y,
                                 rare.normal.frequency = normal.frequency.y,rare.prot =Prot.y,
                                 common.name = Name.x, common.new.name = new.name.x,common.frequency = frequency.x, 
                                 common.normal.frequency = normal.frequency.x,common.prot =Prot.x) %>% 
  mutate(rare_over_common = rare.normal.frequency/common.normal.frequency ) %>% 
  filter(common.prot != '#VALUE!',rare.prot != '#VALUE!')

common.rare  <- transform(common.rare,
                          common.prot = as.numeric(common.prot),
                          rare.prot = as.numeric(rare.prot))


common_rare_summary   <- common.rare %>% 
  group_by(Promoter,RBS,day,lineage) %>%
  summarise(median_rare_over_common = median(rare_over_common),mean_rare_over_common = mean(rare_over_common,
            sd_rare_over_common = sd(rare_over_common)))

ggplot(common_rare_summary,aes(x=day,y=median_rare_over_common,color = lineage)) + 
  geom_line() + 
  geom_point() + 
  facet_grid(Promoter~RBS) +
  coord_cartesian(ylim=c(0,3))


common.rare.C.0  <- common.rare %>% filter(day == 0, lineage == 'C') %>% arrange(Promoter,RBS,Gene)

common.rare.C.16  <- common.rare %>% filter(day == 16, lineage == 'C') %>% arrange(Promoter,RBS,Gene)


common.rare.C.0$RBS <- factor(common.rare.C.0$RBS,
                       levels = c("Strong", "Mid", "Weak", "WT"))

common.rare.C.16$RBS <- factor(common.rare.C.16$RBS,
                     levels = c("Strong", "Mid", "Weak", "WT"))

# 
# max_log = max(log2(test_4$rare.normal.frequency),y=log2(test_4$common.normal.frequency),
#               log2(test_28$rare.normal.frequency),y=log2(test_28$common.normal.frequency))
# 
# min_log = min(log2(test_4$rare.normal.frequency),y=log2(test_4$common.normal.frequency),
#                     log2(test_28$rare.normal.frequency),y=log2(test_28$common.normal.frequency))

#log

  ggplot(common.rare.C.0,aes(x=log2(rare.normal.frequency),y=log2(common.normal.frequency),
                    color = log2(rare.prot/common.prot))) + 
    geom_point(alpha=0.3) +
    facet_grid(RBS~Promoter) +
    coord_fixed(ratio=1
                ,xlim=c(-5,-20),ylim=c(-5,-20)
                ) +
  geom_abline(intercept = 0,slope = 1)  + 
  scale_color_gradient(low="red", high="blue")

  
  
#not log

ggplot(test,aes(x=rare.frequency,y=common.frequency)) + 
  geom_point() +
  facet_grid(RBS~Promoter) +
  coord_fixed(ratio = 1) + 
  expand_limits(y=17500,x=17500)


#compare rare to common over time


rare_more_4  <- test_4 %>% 
  mutate(ratio = rare.normal.frequency/common.normal.frequency) %>%
  group_by(Promoter,RBS) %>% 
  summarise(sum(ratio>1))

rare_more_28  <- test_28 %>%
  mutate(ratio = rare.normal.frequency/common.normal.frequency) %>% 
  group_by(Promoter,RBS) %>%
  summarise(day_28 = sum(ratio>1))


rare_more_c  <- left_join(rare_more_4,rare_more_28,by = c('Promoter','RBS')) %>%  
  gather(day, rare_above_common, day_4:day_28) %>% 
  separate(day, c('word', "day"), '_',convert = T) %>% 
  select(-word)

ggplot(rare_more_c, aes(x = interaction(Promoter,RBS), y = rare_above_common,fill= factor(day))) + geom_bar(stat = "identity",position = 'dodge')

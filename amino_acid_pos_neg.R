
require(seqinr)

names(aacost)[1]  <- 'amino.acid'

pos.neg.amino.acids <- fitseq.data.residuals.above.14 %>% 
  select(day,lineage,Gene,pos.neg,pep.Ala:pep.Tyr) %>%
  mutate_each(funs(.*11),-pos.neg)%>%
  gather(amino.acid,frequency,pep.Ala:pep.Tyr) %>%
  separate(amino.acid,into = c('pep','amino.acid'),'[.]',extra = 'drop') %>%
  select(-pep) %>%
  group_by(day,lineage,pos.neg,amino.acid) %>%
  summarize(gene.frequency = n(),amino.acid.sum = sum(frequency),
            amino.acid.frequency = amino.acid.sum/gene.frequency) %>% 
  select(-gene.frequency,-amino.acid.sum) %>% 
  spread(pos.neg,amino.acid.frequency) %>%
  left_join(aacost %>% select(amino.acid,tot),by = 'amino.acid')
names(pos.neg.amino.acids)[6]  <- 'cost'
  

pos.neg.amino.acids.rbs <- fitseq.data.residuals.above.14 %>%  
  select(day,lineage,Gene,RBS.Display,pos.neg,pep.Ala:pep.Tyr) %>%
  group_by(day,lineage,Gene,RBS.Display,pos.neg) %>%
  mutate_each(funs(.*11),-pos.neg)%>%
  gather(amino.acid,frequency,pep.Ala:pep.Tyr) %>%
  separate(amino.acid,into = c('pep','amino.acid'),'[.]',extra = 'drop') %>%
  select(-pep) %>%
  group_by(day,lineage,RBS.Display,pos.neg,amino.acid) %>%
  summarize(gene.frequency = n(),amino.acid.sum = sum(frequency),
            amino.acid.frequency = amino.acid.sum/gene.frequency) %>% 
  select(-gene.frequency,-amino.acid.sum) %>% 
  spread(pos.neg,amino.acid.frequency) %>%
  filter( RBS.Display != 'WT') %>%
  left_join(aacost %>% select(amino.acid,tot),by = 'amino.acid')
names(pos.neg.amino.acids.rbs)[7]  <- 'cost'

ggplot(pos.neg.amino.acids %>% filter(day == 12),
       aes(y=Positive,x=Negative,color = cost)) +
  geom_point(size = 5) + 
  scale_color_gradient2(low = hue_pal()(3)[3],
                        high = hue_pal()(3)[1],mid = hue_pal()(3)[2],midpoint =45) +
  theme_aviv + 
  geom_abline(slope = 1, intercept = 0) +
  geom_text(aes(label=amino.acid),hjust=0, vjust=0) + 
  coord_equal() + 
  facet_wrap(~lineage)

ggplot(pos.neg.amino.acids.rbs %>% filter(day == 12),
       aes(y=Positive,x=Negative,color = cost)) +
  geom_point() +
  scale_color_gradient2(low = hue_pal()(3)[3],
                        high = hue_pal()(3)[1],mid = hue_pal()(3)[2],midpoint =45)+
  theme_aviv + 
  geom_abline(slope = 1, intercept = 0) +
  geom_text(aes(label=amino.acid),hjust=0, vjust=0) + 
  coord_equal() + 
  facet_grid(RBS.Display~lineage)



#above n pos.neg
pos.neg.amino.acids <- fitseq.data.residuals.above.14 %>% 
  select(day,Gene,pos.neg,pep.Ala:pep.Tyr) %>%
  group_by(day,Gene,pos.neg) %>%
  mutate_each(funs(.*11))%>%
  gather(amino.acid,frequency,pep.Ala:pep.Tyr) %>%
  separate(amino.acid,into = c('pep','amino.acid'),'[.]',extra = 'drop') %>%
  select(-pep) %>%
  group_by(day,pos.neg,amino.acid) %>%
  summarize(gene.frequency = n(),amino.acid.sum = sum(frequency),
            amino.acid.frequency = amino.acid.sum/gene.frequency) %>% 
  select(-gene.frequency,-amino.acid.sum) %>% 
  spread(pos.neg,amino.acid.frequency) %>%
  left_join(aacost %>% select(amino.acid,tot),by = 'amino.acid')
names(pos.neg.amino.acids)[5]  <- 'cost'


pos.neg.amino.acids.rbs <- fitseq.data.residuals.above.14 %>%  
  select(day,Gene,RBS.Display,pos.neg,pep.Ala:pep.Tyr) %>%
  group_by(day,Gene,RBS.Display,pos.neg) %>%
  mutate_each(funs(.*11),-pos.neg)%>%
  gather(amino.acid,frequency,pep.Ala:pep.Tyr) %>%
  separate(amino.acid,into = c('pep','amino.acid'),'[.]',extra = 'drop') %>%
  select(-pep) %>%
  group_by(day,RBS.Display,pos.neg,amino.acid) %>%
  summarize(gene.frequency = n(),amino.acid.sum = sum(frequency),
            amino.acid.frequency = amino.acid.sum/gene.frequency) %>% 
  select(-gene.frequency,-amino.acid.sum) %>% 
  spread(pos.neg,amino.acid.frequency) %>%
  filter( RBS.Display != 'WT') %>%
  left_join(aacost %>% select(amino.acid,tot),by = 'amino.acid')
names(pos.neg.amino.acids.rbs)[6]  <- 'cost'

ggplot(pos.neg.amino.acids %>% filter(day == 12),
       aes(y=Positive,x=Negative,color = cost)) +
  geom_point(size = 5) + 
  scale_color_gradient2(low = hue_pal()(3)[3],
                        high = hue_pal()(3)[1],mid = hue_pal()(3)[2],midpoint =45) +
  theme_aviv + 
  geom_abline(slope = 1, intercept = 0) +
  geom_text(aes(label=amino.acid),hjust=0, vjust=0) + 
  coord_equal() 

ggplot(pos.neg.amino.acids.rbs %>% filter(day == 12),
       aes(y=Positive,x=Negative,color = cost)) +
  geom_point() +
  scale_color_gradient2(low = hue_pal()(3)[3],
                        high = hue_pal()(3)[1],mid = hue_pal()(3)[2],midpoint =45)+
  theme_aviv + 
  geom_abline(slope = 1, intercept = 0) +
  geom_text(aes(label=amino.acid),hjust=0, vjust=0) + 
  coord_equal() + 
  facet_grid(~RBS.Display)



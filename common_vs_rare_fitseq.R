
rare  <- fitseq.data.tidy %>%
  select(Name,new.name, Promoter.Display,RBS.Display,Gene, CDS.type, day,lineage,frequency) %>%
  filter(CDS.type == 'Max Rare')


common  <- fitseq.data.tidy %>%
  select(Name,new.name, Promoter.Display,RBS.Display,Gene, CDS.type, day,lineage,frequency) %>% 
  filter(CDS.type == 'Min Rare')


joined  <- left_join(common,rare,by=c('Gene','RBS.Display','Promoter.Display'))

common.rare <- joined %>% select(Promoter = Promoter.Display, RBS = RBS.Display, Gene, day, lineage, 
                                 rare.name = Name.y,rare.new.name = new.name.y,rare.frequency = frequency.y,
                                 common.name = Name.x,common.new.name = new.name.x,common.frequency = frequency.x)


test  <- common.rare %>% filter(day = 4, lineage = 'C')


ggplot(test,aes(x=log2(1+rare.frequency),y=log2(1+ common.frequency))) + 
  geom_point() +
  facet_grid(RBS~Promoter) +
  coord_fixed(ratio = 1) 

ggplot(test,aes(x=rare.frequency,y=common.frequency)) + 
  geom_point() +
  facet_grid(RBS~Promoter) +
  coord_fixed(ratio = 1) + 
  expand_limits(y=1000,x=1000)
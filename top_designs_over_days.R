#clean data

require(tidyr)
require(ggplot2)
require(dplyr)

result.dir = 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\results\\'

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
                                           -Fltr.SetGood,-full.seq)
fitseq.data.tidy  <- fitseq.data.tidy %>%   gather(sample, frequency, A_24:F_0)
fitseq.data.tidy  <- fitseq.data.tidy %>% separate(sample, c("lineage", "day"), '_',convert = T)
fitseq.data.tidy$RBS.Display <- factor(fitseq.data.tidy$RBS.Display,
                                       levels = c("Strong", "Mid", "Weak", "WT"))

fitseq.data.tidy  <- fitseq.data.tidy %>% 
  group_by(day,lineage) %>%
  mutate(sum.freq_1 = sum(frequency+1),norm.freq_1 = (frequency+1)/sum.freq_1,
         sum.freq = sum(frequency),norm.freq = (frequency)/sum.freq)


top_designs_of_lineage_over_days <- function(min_day,n,lineage_letter){
  
  top.n.fitseq  <-  fitseq.data.tidy %>%
    group_by(day,lineage) %>% 
    mutate(rank = rank(desc(norm.freq),t= "average")) %>%
    filter(rank <= n,lineage==lineage_letter,day >=min_day) %>% select(Name)
  top.n.fitseq.names <- unique(top.n.fitseq$Name)
  fitseq.top.n.all.days <- filter(fitseq.data.tidy,Name %in%top.n.fitseq.names,lineage ==lineage_letter )
  p <-  ggplot(fitseq.top.n.all.days, aes(y=norm.freq,x=day,color = new.name)) +
         geom_line() +
         geom_point()
  return(p)
}

top_designs_of_lineage_over_days(8,10,'D')


 

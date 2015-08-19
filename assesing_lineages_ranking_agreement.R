require(RMySQL)
require(dplyr)
require(tidyr)
require(ggplot2)
require(colorspace)
require(scales)

load.fitseq.data <- function(location){
  
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
                                             -Fltr.SetGood,-full.seq,-pep.sequence,-UTR,-CDS.seq,-Promoter.seq,-RBS.seq,-Promoter,
                                             -variable.seq,-full.peptide,-preRBS,-salis.status)
  fitseq.data.tidy  <- fitseq.data.tidy %>%   gather(sample, frequency, A.24:F.0)
  fitseq.data.tidy  <- fitseq.data.tidy %>% separate(sample, c("lineage", "day"), '.',convert = T)
  fitseq.data.tidy$RBS.Display <- factor(fitseq.data.tidy$RBS.Display,
                                         levels = c("Strong", "Mid", "Weak", "WT"))
  
  
  
  
  fitseq.data.tidy <- fitseq.data.tidy %>% 
    group.by(day,lineage) %>% 
    mutate(sum.sample= sum(frequency),norm.freq = frequency/sum.sample,
           sum.anc= sum(anc.1),norm.anc = anc.1/sum.anc,freq.norm.anc = norm.freq/norm.anc,
           sum.sample.1= sum(frequency+1),norm.freq.1 = (frequency+1)/sum.sample.1,
           sum.anc.1= sum(anc.1+1),norm.anc.1 = (anc.1+1)/sum.anc.1,freq.norm.anc.1 = norm.freq.1/norm.anc.1)
  
  
  return(fitseq.data.tidy)
}


fitseq.data <- load.fitseq.data()


fitseq.data.residuals <- fitseq.data %>% 
  filter(
    day == 12,
    Promoter.Display=='High'
    ,above.14 == TRUE
    ,log2(Prot)  < 17.5
    # ,abs(fit.resid) >  percentile
  )



  



gene.rank.matrix <- fitseq.data.residuals %>% 
  group.by(lineage, Gene) %>% 
  summarize(med.fit = median( fit.resid)) %>%
  group.by(lineage) %>%
  mutate(med.rank = min.rank(desc(med.fit))) %>% 
  arrange(Gene,med.rank) %>% 
  select(lineage,Gene,med.rank) %>% 
  spread(lineage,med.rank) %>% 
  select(-Gene)


gene.rank.matrix <- data.matrix(gene.rank.matrix)
print(kappam.fleiss(gene.rank.matrix))
print(icc(gene.rank.matrix))
print(kripp.alpha(gene.rank.matrix, method="ordinal"))
print(kendall(design.rank.matrix))


design.rank.matrix <- fitseq.data.residuals %>% 
  group.by(lineage,Name) %>% 
  summarize(med.fit = median( fit.resid)) %>%
  group.by(lineage) %>%
  mutate(med.rank = min.rank(desc(med.fit))) %>% 
  arrange(med.rank) %>% 
  select(lineage,Name,med.rank) %>% 
  spread(lineage,med.rank) %>% 
  select(-Name)


design.rank.matrix <- data.matrix(design.rank.matrix)
print(kappam.fleiss(design.rank.matrix))
print(icc(design.rank.matrix))
print(kripp.alpha(design.rank.matrix, method="ordinal"))
print(kendall(design.rank.matrix))

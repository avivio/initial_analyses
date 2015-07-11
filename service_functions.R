

load.packages <- function(){
  require(dplyr)
  require(tidyr)
  require(ggplot2)
  require(colorspace)
  require(scales)
}
get.wilcox.string <-  function(pos,neg){
  
  greater <- wilcox.test(pos,neg,alternative = 'greater')$p.value
  less <- wilcox.test(pos,neg,alternative = 'less')$p.value
  if(greater < less){
    return(paste('\n Positive > Negative\n p value =',formatC(greater,format='e',digits = 2)))
  } else{
    return(paste('\n Positive < Negative\n p value =',formatC(less,format='e',digits = 2)))
  }
  
  
}

load.packages()

theme_aviv <-     theme_minimal() + 
  theme( axis.line = element_line(colour = "black"),
         axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold"),
         strip.text.x = element_text(size=16,face="bold"), strip.text.y = element_text(size=16,face="bold"), 
         plot.title = element_text(size = 25,face = "bold"),
         legend.text = element_text(size=16),legend.title = element_text(size=18,face="bold")	)



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
  fitseq.data.tidy  <- fitseq.data.tidy %>%   gather(sample, frequency, A_24:F_0)
  fitseq.data.tidy  <- fitseq.data.tidy %>% separate(sample, c("lineage", "day"), '_',convert = T)
  fitseq.data.tidy$RBS.Display <- factor(fitseq.data.tidy$RBS.Display,
                                         levels = c("Strong", "Mid", "Weak", "WT"))
  
  
  
  
  fitseq.data.tidy <- fitseq.data.tidy %>% 
    group_by(day,lineage) %>% 
    mutate(sum.sample= sum(frequency),norm.freq = frequency/sum.sample,
           sum.anc= sum(anc_1),norm.anc = anc_1/sum.anc,freq.norm.anc = norm.freq/norm.anc,
           sum.sample.1= sum(frequency+1),norm.freq.1 = (frequency+1)/sum.sample.1,
           sum.anc.1= sum(anc_1+1),norm.anc.1 = (anc_1+1)/sum.anc.1,freq.norm.anc.1 = norm.freq.1/norm.anc.1)
  
  
  return(fitseq.data.tidy)
}

clean.filename <- function(filename){
  filename  <- gsub('\\(','_',filename)
  filename  <- gsub('\\)','',filename)
  filename  <- gsub('\\.','_',filename)
  filename  <- gsub('\\+ 1','_plus_1',filename)
  return(filename)
  
}





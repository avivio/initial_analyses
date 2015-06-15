require(dplyr)
require(tidyr)
require(ggplot2)


result.dir = 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\results\\fitness_vs_x_lineages\\'

fitseq.data.location  <- 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\clean_data\\goodman_salis_tuller_fitseq_2_mismatch_with_0.csv'
fitseq.data = read.csv(fitseq.data.location)
fitseq.data = transform(fitseq.data, Prot = as.numeric(Prot),Bin.Pct.1 = as.numeric(Bin.Pct.1),
                        Count.RNA = as.numeric(Count.RNA))
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

#zero in anc 1 come back later
View(fitseq.data.tidy %>% filter(anc_1 == 0 & day != 0 & frequency >0) %>% 
       select(new.name,anc_1,anc_2,day,lineage,frequency) %>%
       arrange(new.name,day,lineage))
#zero in anc 1 or anc 2
View(fitseq.data.tidy %>% filter(anc_1 == 0 | anc_2 ==0) %>% 
       select(new.name,anc_1,anc_2,day,lineage,frequency) %>%
       arrange(new.name,day,lineage))




fitseq.data.tidy <- fitseq.data.tidy %>% 
  group_by(day,lineage) %>% 
  mutate(sum.sample= sum(frequency),norm.freq = frequency/sum.sample,
         sum.anc= sum(anc_1),norm.anc = anc_1/sum.anc,freq.norm.anc = norm.freq/norm.anc,
         sum.sample.1= sum(frequency+1),norm.freq.1 = (frequency+1)/sum.sample.1,
         sum.anc.1= sum(anc_1+1),norm.anc.1 = (anc_1+1)/sum.anc.1,freq.norm.anc.1 = norm.freq.1/norm.anc.1)



x.strings <- 
  c('Rel.Codon.Freq', 'log2(RNA)','Prot','log2(Trans)','CAI','tAI','CDS.GC','GC','dG','dG.noutr','dG.unif',
    'RNA.FCC',	'Prot.FCC',	'Trans.FCC',	'tAI.DC',	'GC.DC',	'CDS.GC.DC'	,'CAI.DC',
    'RSCU.DC',	'dG.DC',	'dG.noutr.DC',	'dG.unif.DC',	'log2(salis.init)','salis.dG_total',
    'salis.dG_mRNA_rRNA',	'salis.dG_mRNA','TASEP.avgRate',	'TASEP.avgRiboNum',	'TASEP.density.avg',
    'TASEP.density.0',	'TASEP.bottle.neck.position',	'TASEP.bottle.neck.depth')
x.labels <- 
  c('Relative codon frequency', 'Log 2 RNA level','Protein level','Log 2 translation efficiency','CAI','tAI'
    ,'Coding sequence GC %','GC %','Delta G','Delta G noutr','Delta G unif',
    'RNA level log fold change from family mean ','Protein level log fold change from family mean',	
    'Translational effeciency log fold change from family mean',	'tAI DC',	'GC % DC',
    'Coding sequence GC % DC'	,'CAI DC', 'RSCU DC',	'Delta G DC',	'Delta G noutr DC',	'Delta G unif DC',
    'Log 2 RBS calculator initition rate','RBS calculator dG_total','RBS calculator dG_mRNA_rRNA',
    'RBS calculator dG_mRNA','Average translation rate (calculated using TASEP)',	
    'Average ribosome number (calculated using TASEP)',	
    'Average ribosome density per codon (calculated using TASEP)',
    'Density at start codon (calculated using TASEP)',	'Bottle neck position (calculated using TASEP)',
    'Bottle neck depth (calculated using TASEP)')

names(x.labels) <- x.strings

y.string <-  'log2(freq.norm.anc.1)'
y.label <- 'Log 2 of frequency in sample over frequency in ancestor'
# y.string <-  'log2(RNA)'
# y.label <- 'Log 2 of frequency '

x.string  <- 'TASEP.avgRiboNum'
x.label <- 'Ribosome number per mRNA (simulated by TASEP)'


for (i in 1:length(x.labels)){
  x.string <- names(x.labels)[i]
  x.label <- x.labels[i]
  print(x.string)
  print.all.plots.for.x.and.y(fitseq.data.tidy,x.string,y.string,x.label,y.label)
  
}



print.all.plots.for.x.and.y <- function(fitseq.data.tidy,x.string,y.string,x.label,y.label){
  new.dir <- clean.filename(x.string)
  result.dir <- paste0(result.dir,new.dir,'\\')
  dir.create(result.dir)
  
  fitseq.xy <-  fitseq.data.tidy %>%
    mutate_(y = y.string,x = x.string)
  
  fitseq.xy  <- fitseq.xy %>% 
    ungroup() %>%
    group_by(day,lineage) %>%
    mutate(all.r = cor.test(x,y)$estimate,all.p = cor.test(x,y)$p.value, 
           all.cor.text = paste0('r = ', round(all.r,2),'\np-value = ' ,formatC(all.p,format='e',digits = 2)))
  
  fitseq.xy <-  fitseq.xy %>% 
    group_by(day,lineage,Promoter.Display,RBS.Display ) %>%
    mutate(group.r = cor.test(x,y)$estimate,group.p = cor.test(x,y)$p.value,
           group.cor.text = paste0('r = ', round(group.r,2),'\np-value = ' ,formatC(group.p,format='e',digits = 2)))
  
  fitseq.xy <-  fitseq.xy %>% 
    select(x,y,all.cor.text,group.cor.text)
  
  # y.limits <-  fitseq.xy %>%
  #   ungroup() %>%
  #   summarise(min.y = min(y,na.rm = T),max.y = max(y,na.rm = T))
  # y.limits <- c(y.limits$min.y,y.limits$max.y)
  y.limits.all <- c(-5.5,4)
  
  y.limits.rbs.promoter <- c(-4,2)
  
  
  text.loc.y.all <- y.limits.all[2] - ((y.limits.all[2]-y.limits.all[1])/20)
  
  text.loc.y.rbs.promoter <- y.limits.rbs.promoter[2] - ((y.limits.rbs.promoter[2]-y.limits.rbs.promoter[1])/20)
  
  
  x.limits <-  fitseq.xy %>%
    ungroup() %>%
    summarise(min.x = min(x,na.rm = T),max.x = max(x,na.rm = T))
  x.limits <- c(x.limits$min.x,x.limits$max.x)
  
  text.loc.x <- x.limits[1] + ((x.limits[2]-x.limits[1])/100)
  
  fitseq.xy <-  fitseq.xy %>% 
    mutate(text.loc.x = text.loc.x,text.loc.y.all = text.loc.y.all,
           text.loc.y.rbs.promoter = text.loc.y.rbs.promoter )
  
  do.plots.for.all.samples(result.dir,fitseq.xy,x.limits,y.limits.all,y.limits.rbs.promoter,x.label,y.label,x.string,y.string,
                           text.loc.x,text.loc.y.all,text.loc.y.rbs.promoter)
  
}

do.plots.for.all.samples <- function(result.dir,fitseq.xy,x.limits,
                                     y.limits.all,y.limits.rbs.promoter,x.label,y.label,x.string,y.string,
                                     text.loc.x,text.loc.y.all,text.loc.y.rbs.promoter){
  for (lineage.letter in c('A','B','C','D','E','F')){
    for (day.number in seq(4,28,4)){
      get.plots.for.day.lineage(result.dir,fitseq.xy,day.number,lineage.letter,
                                x.limits,y.limits.all,y.limits.rbs.promoter,x.label,y.label,x.string,y.string,
                                text.loc.x,text.loc.y.all,text.loc.y.rbs.promoter)
    }
  }
}


get.plots.for.day.lineage <- function(result.dir,fitseq.xy,day.number,lineage.letter,
                                      x.limits,y.limits.all,y.limits.rbs.promoter,x.label,y.label,x.string,y.string,
                                      text.loc.x,text.loc.y.all,text.loc.y.rbs.promoter){
  
  print(day.number)
  print(lineage.letter)
  fitseq.current.sample <-  fitseq.xy %>% filter(day == day.number, lineage == lineage.letter, RBS.Display != 'WT')
  filename.all  <- filename.promoter.rbs  <- paste('lineage', lineage.letter,'all','day', day.number,'x',x.string ,'vs',y.string,sep = '_')
  filename.all  <-clean.filename(filename.all)
  
  png(paste0(result.dir,filename.all,'.png'),units="in",  width=15, height=12, res=70)
  
  p.all <- plot.scatter.all(fitseq.current.sample,day.number, lineage.letter,
                            x.limits,y.limits.all,x.label,y.label)
  print(p.all)
  dev.off()
  
  filename.promoter.rbs  <- filename.promoter.rbs  <- paste('lineage', lineage.letter,'promoter_rbs','day', day.number,'x',x.string ,'vs',y.string, sep = '_')

  filename.promoter.rbs  <-clean.filename(filename.promoter.rbs)
  
  png(paste0(result.dir,filename.promoter.rbs,'.png'),units="in",  width=15, height=12, res=70)
  
  p.promoter.rbs <- plot.scatter.rbs.promoter(fitseq.current.sample,day.number, lineage.letter,
                                              x.limits,y.limits.rbs.promoter,x.label,y.label,
                                              text.loc.x,text.loc.y.rbs.promoter)
  print(p.promoter.rbs)
  dev.off()
  
}







plot.scatter <-   function(fitseq.current.sample,day.number, lineage.letter,
                           x.limits,y.limits,x.label,y.label) {
  title <- paste(x.label ,'vs\n',y.label,'\nfor day', day.number, 'lineage', lineage.letter)
  
  p <- ggplot(fitseq.current.sample,aes(x= x, y = y)) +
    geom_point(alpha = 0.2)   +
    theme_minimal() + 
    theme( axis.line = element_line(colour = "black"),
           axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold"),
           strip.text.x = element_text(size=16,face="bold"), strip.text.y = element_text(size=16,face="bold"), 
           plot.title = element_text(size = 25,face = "bold"),legend.position="none") +
    # ylim(y.limits) + xlim(x.limits) +
    xlim(x.limits) +
    coord_cartesian(ylim = y.limits) +
    ylab(paste0(y.label,'\n')) + 
    xlab(paste0('\n',x.label)) +
    ggtitle(title) + 
    geom_smooth(color = 'red', size = 1.2, method = 'lm')
  
  
  return(p)
  
}

plot.scatter.rbs.promoter <-   function(fitseq.current.sample,day.number, lineage.letter,x.limits,y.limits,
                                        x.label,y.label,text.loc.x,text.loc.y) {
  title <- paste(x.label ,'vs\n',y.label,'\nfor day', day.number, 'lineage', lineage.letter)
  p <- plot.scatter(fitseq.current.sample,day.number, lineage.letter,x.limits,y.limits,x.label,y.label) +
    facet_grid(Promoter.Display ~ RBS.Display) +
    geom_text(aes(x=text.loc.x, y=text.loc.y.rbs.promoter, colour="red",face="bold", inherit.aes=FALSE, parse=FALSE,hjust = 0,label=group.cor.text), )
  
  return(p)
  
}
plot.scatter.all <-   function(fitseq.current.sample,day.number, lineage.letter,x.limits,y.limits,
                               x.label,y.label) {
  title <- paste(x.label ,'vs\n',y.label,'\nfor day', day.number, 'lineage', lineage.letter)
  p <- plot.scatter(fitseq.current.sample,day.number, lineage.letter,x.limits,y.limits,x.label,y.label) +
    geom_text(aes(x=text.loc.x, y=text.loc.y.all, colour="red",face="bold", 
                  inherit.aes=FALSE, parse=FALSE,hjust = 0,label=all.cor.text), )
  
  return(p)
  
}

clean.filename <- function(filename){
  filename  <- gsub('\\(','_',filename)
  filename  <- gsub('\\)','',filename)
  filename  <- gsub('\\.','_',filename)
  filename  <- gsub('\\+ 1','_plus_1',filename)
  return(filename)
  
}







require(dplyr)
require(tidyr)
require(ggplot2)

theme_aviv <- theme_minimal() + 
  theme( axis.line = element_line(colour = "black"),
         axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold"),
         strip.text.x = element_text(size=16,face="bold"), strip.text.y = element_text(size=16,face="bold"), 
         plot.title = element_text(size = 25,face = "bold"),legend.position="none")

result.dir = 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\results\\fitness_vs_x_lineages\\'

fitseq.data.location  <- 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\clean_data\\goodman_salis_tuller_fitseq_2_mismatch_with_0.csv'
fitseq.data = read.csv(fitseq.data.location)
fitseq.data = transform(fitseq.data, Prot = as.numeric(as.character(Prot))
                        ,Bin.Pct.1 = as.numeric(as.character((Bin.Pct.1))),
                        Count.RNA = as.numeric(as.character((Count.RNA))))
fitseq.data$RBS.Display <- factor(fitseq.data$RBS.Display,
                                       levels = c("Strong", "Mid", "Weak", "WT"))

fitseq.data.no.wt <- fitseq.data %>% filter(RBS.Display != 'WT')

#Count.DNA Count.RNA Count.A.DNA Count.A.RNA Count.B.DNA Count.B.RNA MismatchAvg.DNA MismatchAvg.RNA
# Count.Prot Offset.RBS RNA.A RNA.B RNA Prot Trans RNA.FCC DNA.FCC Prot.FCC Trans.FCC




result.dir = 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\results\\goodman_vs_goodman\\'
y.string <-  'log2(RNA)'
y.label <- 'Log 2 RNA level'


x.strings <- 
  c('Rel.Codon.Freq', 'CAI','tAI','CDS.GC','GC','dG','dG.noutr','dG.unif','tAI.DC',	'GC.DC',	'CDS.GC.DC'	
    ,'CAI.DC','RSCU.DC',	'dG.DC',	'dG.noutr.DC',	'dG.unif.DC',	'log2(salis.init)','salis.dG_total',
    'salis.dG_mRNA_rRNA',	'salis.dG_mRNA','TASEP.avgRate',	'TASEP.avgRiboNum',	'TASEP.density.avg',
    'TASEP.density.0',	'TASEP.bottle.neck.position',	'TASEP.bottle.neck.depth')
x.labels <- 
  c('Relative codon frequency', 'CAI','tAI','Coding sequence GC %','GC %','Delta G','Delta G noutr',
    'Delta G unif','tAI DC',	'GC % DC','Coding sequence GC % DC'	,'CAI DC', 'RSCU DC',	'Delta G DC',
    'Delta G noutr DC',	'Delta G unif DC', 'Log 2 RBS calculator initition rate','RBS calculator dG_total',
    'RBS calculator dG_mRNA_rRNA','RBS calculator dG_mRNA','Average translation rate (calculated using TASEP)',	
    'Average ribosome number (calculated using TASEP)',
    'Average ribosome density per codon (calculated using TASEP)','Density at start codon (calculated using TASEP)',	'Bottle neck position (calculated using TASEP)',
    'Bottle neck depth (calculated using TASEP)')

names(x.labels) <- x.strings

for (i in 1:length(x.labels)){
  x.string <- names(x.labels)[i]
  x.label <- x.labels[i]
  print(x.string)
  print.all.plots.for.x.and.y.goodman(fitseq.data,x.string,y.string,x.label,y.label)
  
}




print.all.plots.for.x.and.y.goodman <- function(fitseq.data.tidy,x.string,y.string,x.label,y.label){
  dir.create(result.dir)
  
  fitseq.xy <-  fitseq.data %>%
    mutate_(y = y.string,x = x.string)
  
  fitseq.xy  <- fitseq.xy %>% 
    ungroup() %>%
    mutate(all.r = cor.test(x,y)$estimate,all.p = cor.test(x,y)$p.value, 
           all.cor.text = paste0('r = ', round(all.r,2),'\np-value = ' ,formatC(all.p,format='e',digits = 2)))
  
  fitseq.xy <-  fitseq.xy %>% 
    group_by(Promoter.Display,RBS.Display ) %>%
    mutate(group.r = cor.test(x,y)$estimate,group.p = cor.test(x,y)$p.value,
           group.cor.text = paste0('r = ', round(group.r,2),'\np-value = ' ,formatC(group.p,format='e',digits = 2)))
  
  fitseq.xy <-  fitseq.xy %>% 
    select(x,y,all.cor.text,group.cor.text)
  
  y.limits <-  fitseq.xy %>%
    ungroup() %>%
    summarise(min.y = min(y,na.rm = T),max.y = max(y,na.rm = T))
  y.limits <- c(y.limits$min.y,y.limits$max.y)
#   y.limits.all <- c(-5.5,4)
#   
#   y.limits.rbs.promoter <- c(-4,2)
  y.limits.all <- y.limits
  y.limits.rbs.promoter <- y.limits
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
  
  get.plots.goodman.data(result.dir,fitseq.xy,
                            x.limits,y.limits.all,y.limits.rbs.promoter,x.label,y.label,x.string,y.string,
                            text.loc.x,text.loc.y.all,text.loc.y.rbs.promoter)
  
}



get.plots.goodman.data <- function(result.dir,fitseq.xy,
                                      x.limits,y.limits.all,y.limits.rbs.promoter,x.label,y.label,x.string,y.string,
                                      text.loc.x,text.loc.y.all,text.loc.y.rbs.promoter){
  
  fitseq.current.sample <-  fitseq.xy %>% filter(RBS.Display != 'WT')
  filename.all  <- filename.promoter.rbs  <- paste('x',x.string ,'vs',y.string,sep = '_')
  filename.all  <-clean.filename(filename.all)
  
  png(paste0(result.dir,filename.all,'.png'),units="in",  width=15, height=12, res=70)
  
  p.all <- plot.scatter.all.goodman(fitseq.current.sample,x.limits,y.limits.all,x.label,y.label)
  print(p.all)
  dev.off()
  
  filename.promoter.rbs  <- filename.promoter.rbs  <- paste('promoter_rbs','x',x.string ,'vs',y.string, sep = '_')
  
  filename.promoter.rbs  <-clean.filename(filename.promoter.rbs)
  
  png(paste0(result.dir,filename.promoter.rbs,'.png'),units="in",  width=15, height=12, res=70)
  
  p.promoter.rbs <- plot.scatter.rbs.promoter.goodman(fitseq.current.sample,x.limits,y.limits.rbs.promoter,x.label,y.label,
                                              text.loc.x,text.loc.y.rbs.promoter)
  print(p.promoter.rbs)
  dev.off()
  
}

plot.scatter.rbs.promoter.goodman <-   function(fitseq.current.sample,x.limits,y.limits,
                                        x.label,y.label,text.loc.x,text.loc.y) {
  title <- paste(x.label ,'vs\n',y.label)
  p <- plot.scatter(fitseq.current.sample,x.limits,y.limits,x.label,y.label,title) +
    facet_grid(Promoter.Display ~ RBS.Display) +
    geom_text(aes(x=text.loc.x, y=text.loc.y.rbs.promoter, colour="red",face="bold", inherit.aes=FALSE, parse=FALSE,hjust = 0,label=group.cor.text), )
  
  return(p)
  
}
plot.scatter.all.goodman <-   function(fitseq.current.sample,x.limits,y.limits,
                               x.label,y.label) {
  title <- paste(x.label ,'vs\n',y.label)
  p <- plot.scatter(fitseq.current.sample,x.limits,y.limits,x.label,y.label,title) +
    geom_text(aes(x=text.loc.x, y=text.loc.y.all, colour="red",face="bold", 
                  inherit.aes=FALSE, parse=FALSE,hjust = 0,label=all.cor.text), )
  
  return(p)
  
}

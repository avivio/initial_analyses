require(dplyr)
require(tidyr)
require(ggplot2)

theme_aviv <-     theme_minimal() + 
  theme( axis.line = element_line(colour = "black"),
         axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold"),
         strip.text.x = element_text(size=16,face="bold"), strip.text.y = element_text(size=16,face="bold"), 
         plot.title = element_text(size = 25,face = "bold"),
         legend.text = element_text(size=16),legend.title = element_text(size=18,face="bold")	)


get.wilcox.string <-  function(pos,neg){
  
  greater <- wilcox.test(pos,neg,alternative = 'greater')$p.value
  less <- wilcox.test(pos,neg,alternative = 'less')$p.value
  if(greater < less){
    return(paste('Positive > Negative\np value =',formatC(greater,format='e',digits = 2)))
  } else{
    return(paste('Positive < Negative\np value =',formatC(less,format='e',digits = 2)))
  }
  
  
}


fitseq.data.location  <- 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\clean_data\\goodman_salis_tuller_fitseq_2_mismatch_with_0.csv'
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
                                           -Fltr.SetGood,-full.seq)
fitseq.data.tidy  <- fitseq.data.tidy %>%   gather(sample, frequency, A_24:F_0)
fitseq.data.tidy  <- fitseq.data.tidy %>% separate(sample, c("lineage", "day"), '_',convert = T)
fitseq.data.tidy$RBS.Display <- factor(fitseq.data.tidy$RBS.Display,
                                       levels = c("Strong", "Mid", "Weak", "WT"))



day.number = 12
lineage.letter = 'B'
for (lineage.letter in c('A','B','C','D','E','F')){
  print(lineage.letter)

fitseq.data.tidy.letter <-  fitseq.xy %>% 
  filter(above.14 == T,Promoter.Display == 'High',log2(Prot) < 17.5)  %>% 
  filter(day == day.number, lineage == lineage.letter)

fitseq.data.tidy.letter <- fitseq.data.tidy.letter %>% 
  group_by(day,lineage) %>% 
  mutate(sum.sample= sum(frequency),norm.freq = frequency/sum.sample,
         sum.anc= sum(anc_1),norm.anc = anc_1/sum.anc,freq.norm.anc = norm.freq/norm.anc,
         sum.sample.1= sum(frequency+1),norm.freq.1 = (frequency+1)/sum.sample.1,
         sum.anc.1= sum(anc_1+1),norm.anc.1 = (anc_1+1)/sum.anc.1,freq.norm.anc.1 = norm.freq.1/norm.anc.1)


data.lm.h = lm(log2(freq.norm.anc.1)~log2(Prot), data =filter(fitseq.current.sample, above.14 == TRUE) )
  
slope.h <- signif(data.lm.h$coef[[2]], 5)
p.value.h <- summary(data.lm.h)$coef[2,4]
rmsd.h = signif(sqrt(mean(data.lm.h$residuals^2)))

lm.text.h <- paste('Slope:',slope.h,'\n','p value:',p.value.h,'\n','RMSD:',rmsd.h )

data.lm.l = lm(log2(freq.norm.anc.1)~log2(Prot), data =filter(fitseq.current.sample, above.14 == FALSE) )

slope.l <- signif(data.lm.l$coef[[2]], 5)
p.value.l <- summary(data.lm.l)$coef[2,4]
rmsd.l = signif(sqrt(mean(data.lm.l$residuals^2)))

lm.text.l <- paste('Slope:',slope.l,'\n','p value:',p.value.l,'\n','RMSD:',rmsd.l )

fitseq.data.tidy.letter <- fitseq.data.tidy.letter %>% mutate(low.text = lm.text.l ,high.text = lm.text.h,
                                                high.x = 14,high.y =2,low.x = 14,low.y =2 )


result.dir = 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\results\\minerva\\test\\'
dir.create(result.dir)
title <- 'Log 2 expression level vs log 2 frequency compared to ancesotr\n for day 12 lineage C'
png(paste0(result.dir,'prot_vs_fit_col_prom.png'),units="in",  width=15, height=12, res=70)

# filter(fitseq.data.tidy.letter, Promoter.Display == 'High')
p <- ggplot(fitseq.data.tidy.letter,aes(x= log2(Prot), y = log2(freq.norm.anc.1)
                                                                     ,color  = Promoter.Display)) +
  geom_point(alpha = 0.3)   +
  coord_cartesian(ylim = c(-5.5,4)) +
  xlab(paste0('Log 2 Protein level','\n')) + 
  ylab(paste0('\n','Log 2 of frequency in sample over frequency in ancestor')) +
  ggtitle(title) + 
  guides(color = guide_legend(title = "Promoter"),fill =guide_legend(title = "Promoter")) +
  geom_smooth( size = 1.2, method = 'lm',fill = NA) +
#   geom_text(aes(x=high.x, y=high.y, face="bold", inherit.aes=FALSE,
#                 parse=FALSE,hjust = 0,label=high.text))  + 
#   geom_text(aes(x=low.x, y=low.y, face="bold", inherit.aes=FALSE,
#                   parse=FALSE,hjust = 0,label=low.text)) +
  theme_aviv

print(p)

dev.off()


fitseq.data.high <- fitseq.data.tidy.letter %>% filter(Promoter.Display == 'High')


fitseq.data.models <- fitseq.data.high %>% 
  group_by(lineage) %>% 
  do(lin =  lm(freq.norm.anc.1 ~ Prot, data=.)) %>%
  mutate(intercept = coefficients(lin)[1],slope = coefficients(lin)[2] ) %>%
  select(-lin)


fitseq.data.prediction <- left_join(fitseq.data.high,fitseq.data.models,by= 'lineage') %>% 
  mutate(prediction =intercept + slope*Prot,
         fit.resid = freq.norm.anc.1 -  prediction, 
         pos.neg = ifelse(fit.resid > 0, 'Positive', 'Negative') ) %>%
  filter(!is.na(pos.neg))


# x.strings <- 
#   c('Rel.Codon.Freq','CAI','tAI','CDS.GC','dG','log2(salis.init)',
#     'TASEP.avgRiboNum',	'TASEP.density.0',	'TASEP.bottle.neck.position',	'TASEP.bottle.neck.depth')
# x.labels <- 
#   c('Relative codon frequency','CAI','tAI','Coding sequence GC %','Delta G','Log 2 RBS calculator initition rate',	
#     'Average ribosome number (calculated using TASEP)',	'Density at start codon (calculated using TASEP)',
#     'Bottle neck position (calculated using TASEP)','Bottle neck depth (calculated using TASEP)')

# x.strings <- 
#   c('CDS.GC','GC','dG','dG.noutr','dG.unif')
# x.labels <- 
#   c('Coding sequence\nGC content','GC content','Delta G','Delta G\nno UTR','Delta G\nsequence cut at -5')

x.strings <- 
  c('tAI')
x.labels <- 
  c('tAI')



names(x.labels) <- x.strings
for (i in 1:length(x.labels)){
  x.string  <- names(x.labels)[i]
  x.label <- x.labels[i]
  print(x.label)
title <- paste('Comparing',x.label, 'between positive and negative fitness residuals','\nLineage', lineage.letter)
png(paste0(result.dir,'all_pos_neg_hist_' ,x.string,'_lineage_',lineage.letter, '.png'),units="in",  width=15, height=12, res=70)

neg  <- fitseq.data.prediction %>%
  ungroup %>% 
  filter(pos.neg== 'Negative') %>%
  select_(x.string) 
neg <- as.numeric( unlist(neg))  
pos  <- fitseq.data.prediction %>%
  ungroup %>% 
  filter(pos.neg== 'Positive') %>%
  select_(x.string) 
pos <- as.numeric( unlist(pos)) 

wilcox.string <- get.wilcox.string(pos,neg)

x.limits <-  fitseq.xy %>%
  ungroup() %>%
  mutate_(x = x.string) %>%
  summarise(min.x = min(x,na.rm = T),max.x = max(x,na.rm = T))
x.limits <- c(x.limits$min.x,x.limits$max.x)
text.loc.x <- x.limits[1] + ((x.limits[2]-x.limits[1])/100)


  p <- ggplot(fitseq.data.prediction, aes_string(x = x.string,fill = 'pos.neg')) +
    geom_density(alpha = 0.4) +
    xlab(NULL)+
    ylab(NULL) +  
    #     ggtitle(title) +
    annotate("text", y = Inf, x = -Inf, label = wilcox.string, color="black",hjust = 0,vjust = 1, size = 6,
             parse=FALSE,fontface = 'bold',colour="red") +
    guides(fill=guide_legend(title='Fitness residual')) +
    theme_aviv
  print(p)
  
dev.off()  
  
png(paste0(result.dir,'rbs_pos_neg_hist_' ,x.string,'_lineage_',lineage.letter, '.png'),units="in",  width=15, height=12, res=70)


p <- ggplot(filter(fitseq.data.prediction,RBS.Display != 'WT'), aes_string(x = x.string,fill = 'pos.neg')) +
geom_density(alpha = 0.4) +
facet_grid(RBS.Display ~ .) +
xlab(NULL)+
ylab(NULL) +  
# ggtitle(title)+
guides(fill=guide_legend(title='Fitness residual')) +
theme_aviv
print(p)

dev.off()
}
}

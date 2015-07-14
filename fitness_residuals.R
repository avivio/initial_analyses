require(dplyr)
require(tidyr)
require(ggplot2)
theme_aviv <-     theme_minimal() + 
  theme( axis.line = element_line(colour = "black"),
         axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold"),
         strip.text.x = element_text(size=16,face="bold"), strip.text.y = element_text(size=16,face="bold"), 
         plot.title = element_text(size = 25,face = "bold"))

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




fitseq.data.tidy <- fitseq.data.tidy %>% 
  group_by(day,lineage) %>% 
  mutate(sum.sample= sum(frequency),norm.freq = frequency/sum.sample,
         sum.anc= sum(anc_1),norm.anc = anc_1/sum.anc,freq.norm.anc = norm.freq/norm.anc,
         sum.sample.1= sum(frequency+1),norm.freq.1 = (frequency+1)/sum.sample.1,
         sum.anc.1= sum(anc_1+1),norm.anc.1 = (anc_1+1)/sum.anc.1,freq.norm.anc.1 = norm.freq.1/norm.anc.1)


x.strings <- 
  c('Rel.Codon.Freq','CAI','tAI','CDS.GC','dG','log2(salis.init)',
    'TASEP.avgRiboNum',	'TASEP.density.0',	'TASEP.bottle.neck.position',	'TASEP.bottle.neck.depth')
x.labels <- 
  c('Relative codon frequency','CAI','tAI','Coding sequence GC %','Delta G','Log 2 RBS calculator initition rate',	
    'Average ribosome number (calculated using TASEP)',	'Density at start codon (calculated using TASEP)',
    'Bottle neck position (calculated using TASEP)','Bottle neck depth (calculated using TASEP)')

names(x.labels) <- x.strings


# protein.fitness.lm <-  lm(freq.norm.anc.1 ~ Prot, data=fitseq.data.tidy)
# protein.fitness.coeffs <- coefficients(protein.fitness.lm)
# 
# fitseq.data.prediction <- fitseq.data.tidy  %>% 
#   mutate(prediction =protein.fitness.coeffs[1] + protein.fitness.coeffs[2]*Prot,
#          fit.resid = freq.norm.anc.1 -  prediction)

day.number = 12
fitseq.data.tidy.day <- fitseq.data.tidy %>% filter(day == day.number)

fitseq.data.models <- fitseq.data.tidy.day %>% 
  group_by(lineage) %>% 
  do(lin =  lm(freq.norm.anc.1 ~ Prot, data=.)) %>%
  mutate(intercept = coefficients(lin)[1],slope = coefficients(lin)[2] ) %>%
  select(-lin)


fitseq.data.prediction <- left_join(fitseq.data.tidy,fitseq.data.models,by= 'lineage') %>% 
  mutate(prediction =intercept + slope*Prot,
         fit.resid = freq.norm.anc.1 -  prediction)



fitseq.data.prediction <- fitseq.data.prediction %>% filter(fit.resid < 2.6) 
y.limits <-  fitseq.data.prediction %>%
  ungroup() %>%
  summarise(min.y = min(fit.resid,na.rm = T),max.y = max(fit.resid,na.rm = T))
y.limits <- c(y.limits$min.y,y.limits$max.y)

result.dir = 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\results\\fitness_residuals\\'
dir.create(result.dir)


for (lineage.letter in c('A','B','C','D','E','F')){
  print(lineage.letter)


fitseq.current.sample <-  fitseq.data.prediction %>% filter(day == day.number, lineage == lineage.letter)

fitseq.current.sample.no.wt <- fitseq.current.sample %>% 
  filter(RBS.Display!='WT')



filename.all <- paste('all_fitness_residual_histogram_lineage', lineage.letter,'day', day.number,sep = '_' )

title <- paste('Fitness residual histogram \nfor lineage', lineage.letter,'day', day.number)

png(paste0(result.dir,filename.all,'.png'),units="in",  width=15, height=12, res=70)

p <- ggplot(fitseq.current.sample, aes(fit.resid)) + 
  geom_histogram() + 
  theme_aviv +
  xlab('Fitness residuals') +
  ggtitle(title) +
  expand_limits(y = c(0,3000))
  

print(p)
dev.off()

filename.rbs.promoter <- paste('rbs_promoter_fitness_residual_histogram_lineage', lineage.letter,'day', day.number,sep = '_' )

png(paste0(result.dir,filename.rbs.promoter,'.png'),units="in",  width=15, height=12, res=70)


p <- ggplot(fitseq.current.sample.no.wt, aes(fit.resid)) + 
  geom_histogram() + 
  theme_aviv +
  xlab('Fitness residuals') +
  ggtitle(paste('Fitness residual histogram \nfor lineage', lineage.letter,'day', day.number)) +
  facet_grid(Promoter.Display~RBS.Display)+
  expand_limits(y = c(0,650))


print(p)
dev.off()

y.string <-  'fit.resid'
y.label <- 'Fitness Residual'




for (i in 1:length(x.labels)){
  x.string <- names(x.labels)[i]
  x.label <- x.labels[i]
  print(x.label)
  filename <- paste(x.string ,'vs',y.string,'lineage',lineage.letter,sep = '_')
  title <- paste(x.label ,'vs',y.label,'\n','lineage',lineage.letter,'day', day.number)
  print.all.plots.for.x.and.y.goodman(result.dir,fitseq.current.sample,x.string,y.string,x.label,y.label,
                                      filename,title,y.limits)
  
}

}
#load packages and data
load.packages()

get.data.for.correlation.category <- function(source.data,y.string,y.label,x.string,x.label,col.string,col.label){
  
  
  
  
  fitseq.xy <-  source.data %>%
    mutate_(y = y.string,x = x.string,col = col.string)
  
  fitseq.xy <- fitseq.xy %>%
    filter(!is.na(y) & !is.na(x) & !is.na(col))
  
  fitseq.xy.lm.all <- fitseq.xy %>% 
    ungroup() %>%
    group_by(day,lineage,col) %>% 
    do(mod = lm(y ~ x, data = .)) %>%
    mutate(slope = summary(mod)$coeff[2],
           intercept = summary(mod)$coeff[1],
           rmsd = sqrt(mean(mod$residuals^2)),
           p.value = summary(mod)$coeff[2,4],
           lm.string =paste('\n',  'slope = ', signif(slope,2),'\n',
                            'p-value = ' ,formatC(p.value,format='e',digits = 2),'\n',
                            'RMSD = ', signif(rmsd,2),'\n')) %>% 
    select(-mod)
  
  fitseq.xy <- left_join(fitseq.xy,fitseq.xy.lm.all,by = c('day','lineage','col'))
  
  fitseq.xy.lm.group <- fitseq.xy %>% 
    ungroup() %>%
    group_by(day,lineage,col,Promoter.Display,RBS.Display ) %>%
    do(mod = lm(y ~ x, data = .)) %>%
    mutate(slope = summary(mod)$coeff[2],
           intercept = summary(mod)$coeff[1],
           rmsd = sqrt(mean(mod$residuals^2)),
           p.value = summary(mod)$coeff[2,4],
           lm.string.group = paste('\n',  'slope = ', signif(slope,2),'\n',
                                   'p-value = ' ,formatC(p.value,format='e',digits = 2),'\n',
                                   'RMSD = ', signif(rmsd,2),'\n')) %>% 
    select(day,lineage,col,Promoter.Display,RBS.Display ,lm.string.group)
  
  
  fitseq.xy <- left_join(fitseq.xy,fitseq.xy.lm.group,by = c('day','lineage','col','Promoter.Display','RBS.Display'))
  
  
  fitseq.xy <-  fitseq.xy %>% 
    mutate(text.loc.y = ifelse(above.14,-Inf,Inf),
           text.loc.x= ifelse(above.14,Inf,-Inf),
           hor = ifelse(above.14,1,0),
           ver = ifelse(above.14,0,1))
  
  fitseq.xy <-  fitseq.xy %>% 
    mutate(prediction = intercept + slope*log2(Prot),
           fit.resid = log2(freq.norm.anc.1) -  prediction, 
           pos.neg = ifelse(fit.resid > 0, 'Positive', 'Negative') )
  return(fitseq.xy)
  
}
fitseq.data.location  <- 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\clean_data\\goodman_salis_tuller_fitseq_2_mismatch_with_0.csv'
fitseq.data.tidy  <- load.fitseq.data(fitseq.data.location)
fitseq.data.tidy <- fitseq.data.tidy %>% select(-CDS.seq,-Promoter.seq,-RBS.seq,-Promoter,-RBS,
                                                -variable.seq,-full.peptide,salis.status)

fitseq.data.tidy <- fitseq.data.tidy %>% 
  mutate(above.14 = log2(Prot)> 14)


#define parameters for correlation 

y.string <-  'log2(freq.norm.anc.1)'
y.label <- 'Log 2 of frequency in sample over frequency in ancestor'
x.string  <- 'log2(Prot)'
x.label <- 'Log 2 Protein level'
col.string  <- 'above.14'
# col.label <- 'Number of ribosomes\nper mRNA \n(calulated by TASEP)'
col.label <- 'Protein level\nabove 14'

# get correlation data into data set
fitseq.data.residuals <- get.data.for.correlation.category(fitseq.data.tidy,y.string,y.label,x.string,x.label,col.string,col.label)




base.result.dir = 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\results\\fitness_residuals\\'
for (lineage.letter in c('A','B','C','D','E','F')){
  
png(paste0(base.result.dir,'correlation_over_time_lineage_',lineage.letter,'.png'),units="in",  width=15, height=12, res=70)
title <- paste(x.label ,'vs',y.label,'\nfor day', day.number,'lineage', lineage.letter)
p <- ggplot(fitseq.data.residuals %>% filter(lineage == lineage.letter,day > 0),aes(x= x, y = y,color  = col)) +
  geom_point(alpha = 0.3)   + facet_wrap(~day) +
  ylab(paste0(y.label,'\n')) + 
  xlab(paste0('\n',x.label)) +
  ggtitle(title) + 
  geom_smooth(aes_string(linetype = col.string),show_guide  = F,
              color = 'red', size = 1.2, method = 'lm') +
  geom_text(aes(x=text.loc.x, y=text.loc.y, face="bold", inherit.aes=FALSE,
                parse=FALSE,label=lm.string,colour=col,hjust=hor, vjust=ver),show_guide  = F) +
  guides(colour = guide_legend(override.aes = list(alpha = 1),title = col.label)) +
  theme_aviv
print(p)
dev.off()
}
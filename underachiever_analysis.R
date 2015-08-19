  get.wilcox.string.underachievers  <-  function(label.1,label.2,legend.string.1,legend.string.2){
    
    greater <- wilcox.test(legend.string.1,legend.string.2,alternative = 'greater')$p.value
    less <- wilcox.test(legend.string.1,legend.string.2,alternative = 'less')$p.value
    if(greater < less){
      return(paste(label.1, ' > ', label.2,'\n p value =',formatC(greater,format='e',digits = 2)))
    } else{
      return(paste(label.1, ' < ', label.2,'\n p value =',formatC(less,format='e',digits = 2)))
    }
    
    
  }
  
  compare.pos.neg.all.lineages.underachievers  <- function(fitseq.data.prediction,variable,label,day.number,data.set.name
                                            ,result.dir,legend.string,legend.label){
    
    fitseq.current.sample <-  fitseq.data.prediction %>% 
      mutate_(legend.string = legend.string)%>%
      filter(day == day.number)
    fitseq.current.sample.no.wt <-  fitseq.current.sample %>% 
      filter(RBS.Display != 'WT')
    
    
    legend.string.negative  <- fitseq.current.sample %>%
      ungroup %>% 
      filter(legend.string== 'Negative') %>%
      mutate_(new_variable = variable)  %>% select(new_variable) 
    legend.string.negative <- as.numeric( unlist(legend.string.negative))  
    legend.string.positive  <- fitseq.current.sample %>%
      ungroup %>% 
      filter(legend.string== 'Positive') %>%
      mutate_(new_variable = variable)  %>% select(new_variable) 
    legend.string.positive <- as.numeric( unlist(legend.string.positive)) 
    
    legend.string.under  <- fitseq.current.sample %>%
      ungroup %>% 
      filter(legend.string== 'Underachiever') %>%
      mutate_(new_variable = variable)  %>% select(new_variable) 
    legend.string.under <- as.numeric( unlist(legend.string.under)) 
    
    
    wilcox.string.all.pos.neg <- get.wilcox.string.underachievers('Positive','Negative',legend.string.positive,legend.string.negative)
    wilcox.string.all.neg.under <- get.wilcox.string.underachievers('Negative','Underachiever',legend.string.negative,legend.string.under)
    
    
    fitseq.current.sample <- fitseq.current.sample %>%
      mutate(wilcox.all = paste('',wilcox.string.all.pos.neg,wilcox.string.all.neg.under,sep = '\n '))
    
    rbs.list = c('Strong','Mid','Weak')
    wilcox.string.rbs.pos.neg = vector(mode = 'character',length = 3)
    wilcox.string.rbs.neg.under = vector(mode = 'character',length = 3)
    names(wilcox.string.rbs.pos.neg) <- rbs.list
    names(wilcox.string.rbs.neg.under) <- rbs.list
    
    
    for (rbs in rbs.list){
      legend.string.negative.rbs  <- fitseq.current.sample.no.wt %>%
        ungroup %>% 
        filter(legend.string== 'Negative',RBS.Display == rbs) %>%
        mutate_(new_variable = variable)  %>% select(new_variable)  
      legend.string.negative.rbs <- as.numeric( unlist(legend.string.negative.rbs))  
      legend.string.positive.rbs  <- fitseq.current.sample.no.wt %>%
        ungroup %>% 
        filter(legend.string== 'Positive',RBS.Display == rbs) %>%
        mutate_(new_variable = variable)  %>% select(new_variable) 
      legend.string.positive.rbs <- as.numeric( unlist(legend.string.positive.rbs)) 
      legend.string.under.rbs  <- fitseq.current.sample.no.wt %>%
        ungroup %>% 
        filter(legend.string== 'Underachiever',RBS.Display == rbs) %>%
        mutate_(new_variable = variable)  %>% select(new_variable) 
      legend.string.under.rbs <- as.numeric( unlist(legend.string.under.rbs)) 
      
      wilcox.string.rbs.pos.neg[rbs] <- get.wilcox.string.underachievers('Positive','Negative',legend.string.positive.rbs,legend.string.negative.rbs)
      wilcox.string.rbs.neg.under[rbs] <- get.wilcox.string.underachievers('Negative','Underachiever',legend.string.negative.rbs,legend.string.under.rbs)
      
      
    }
  
    
      fitseq.current.sample.no.wt <- fitseq.current.sample.no.wt %>%
      mutate(wilcox.rbs = paste('',wilcox.string.rbs.pos.neg[RBS.Display],
                                wilcox.string.rbs.neg.under[RBS.Display],sep= '\n '))
    
    title <- paste('Comparing',label, 'between positive, negative and underachiever',legend.label,
                   '\nDay', day.number, gsub('_',' ',data.set.name))
    
    png(paste0(result.dir,'all_pos_neg_hist_' ,variable,'_',data.set.name, '.png'),
        units="in",  width=15, height=12, res=70)
    p <- ggplot(fitseq.current.sample, aes_string(x = variable,fill =legend.string )) +
      geom_density(alpha = 0.4) +
      xlab(paste0(label,'\n'))+
      # expand_limits(y = y.lim.all  ) +
      ylab('\nDensity') +  
      ggtitle(title) +
      geom_text(aes( y = Inf, x = -Inf, label = wilcox.all), color="black",hjust = 0,vjust = 1, size = 6,
                parse=FALSE,fontface = 'bold',colour="red") +
      guides(fill=guide_legend(title=legend.label)) +
      theme_aviv
    print(p)
    dev.off()
    
    png(paste0(result.dir,'promoter_rbs_pos_neg_hist_' ,variable,'_',data.set.name, '.png'),
        units="in",  width=15, height=12, res=70)
    p <- ggplot(fitseq.current.sample.no.wt, aes_string(x = variable,fill =legend.string)) +
      geom_density(alpha = 0.4) +
      facet_grid(RBS.Display ~ .) +
      # expand_limits(y = y.lim.rbs  ) +
      xlab(paste0(label,'\n')) +
      ylab('\nDensity') +  
      ggtitle(title) +
      geom_text(aes( y = Inf, x = -Inf, label = wilcox.rbs), color="black",hjust = 0,vjust = 1, size = 6,
                parse=FALSE,fontface = 'bold',colour="red") +
      guides(fill=guide_legend(title=legend.label)) +
      theme_aviv
    print(p)
    dev.off()
    
    
    
  }
  
  

filter.n.fitness.residuals.underachiever <- function(n,fitseq.data.residuals){
  
  fitseq.data.pos.neg.lineages <- 
    fitseq.data.residuals %>%
    select(Name, day, lineage,pos.neg) %>% 
    spread(lineage,pos.neg)
  
  n.pos <- fitseq.data.pos.neg.lineages %>%
    mutate(A= as.numeric(A=='Positive'),B= as.numeric(B=='Positive'),C= as.numeric(C=='Positive'),
           D= as.numeric(D=='Positive'),E=as.numeric( E=='Positive'),F= as.numeric(F=='Positive')) %>%
    mutate(sum.lineages = A+B+C+D+E+F,n.pos.neg = ifelse(sum.lineages >= n,'Positive',NA) ) %>%
    select(Name,day,n.pos.neg) %>%
    filter(!is.na(n.pos.neg))
  
  n.neg <- fitseq.data.pos.neg.lineages %>%
    mutate(A= as.numeric(A=='Negative'),B= as.numeric(B=='Negative'),C= as.numeric(C=='Negative'),
           D= as.numeric(D=='Negative'),E=as.numeric( E=='Negative'),F= as.numeric(F=='Negative')) %>%
    mutate(sum.lineages = A+B+C+D+E+F,n.pos.neg = ifelse(sum.lineages >= n,'Negative',NA) ) %>%
    select(Name,day,n.pos.neg) %>%
    filter(!is.na(n.pos.neg))
  
  n.pos.neg <-  bind_rows(n.pos,n.neg)
  
  
  fitseq.data.underacheiver.lineages <- 
    fitseq.data.residuals %>%
    select(Name, lineage,underachiever)  %>%  
    spread(lineage,underachiever)  
  
  n.under <- fitseq.data.underacheiver.lineages %>%
    mutate(A= as.numeric(A==TRUE),B= as.numeric(B==TRUE),C= as.numeric(C==TRUE),
           D= as.numeric(D==TRUE),E=as.numeric( E==TRUE),F= as.numeric(F==TRUE)) %>%
    mutate(sum.lineages = A+B+C+D+E+F,n.underachiever = ifelse(sum.lineages >= n,TRUE,FALSE)) %>% 
    select(Name,day,n.underachiever)
    
  
    n.pos.neg <- inner_join(n.pos.neg,n.under, by = c('day','Name'))

    fitseq.data.residuals <-  inner_join(fitseq.data.residuals,n.pos.neg, by = c('day','Name'))
    
    fitseq.data.residuals <- fitseq.data.residuals %>% 
      filter(lineage =='A') %>% 
      ungroup() %>% 
      select(-lineage)
    
  
  return(fitseq.data.residuals)
  
  
}

load.packages()

fitseq.data.location  <- 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\fitseq_data_with_meta_day_12_no_umi_05_08.csv'
fitseq.data.residuals  <- load.fitseq.data(fitseq.data.location)

fitseq.data.residuals <- fitseq.data.residuals %>% 
  mutate(above.14 = log2(Prot)> 14)

#adding something tiny for commit
#define parameters for correlation 

y.string <-  'log2(freq.norm.anc.1)'
y.label <- 'Log 2 of frequency in sample over frequency in ancestor'
x.string  <- 'log2(Prot)'
x.label <- 'Log 2 Protein level'
col.string  <- 'above.14'
# col.label <- 'Number of ribosomes\nper mRNA \n(calulated by TASEP)'
col.label <- 'Protein level\nabove 14'




fitseq.data.residuals <- fitseq.data.residuals %>% 
  filter(     log2(Prot) < 17.5,
              
              Promoter.Display=='High')

# get correlation data into data set
fitseq.data.residuals <- get.data.for.correlation.category(fitseq.data.residuals,y.string,y.label,x.string,x.label,col.string,col.label)




fitseq.data.residuals.below.14 <- fitseq.data.residuals %>% filter(above.14 == F)


fit.summary  <- fitseq.data.residuals.below.14 %>% 
  group_by(lineage) %>% 
  summarise(percentile= quantile(fit.resid,probs = c(0.025))) 
fitseq.data.residuals.below.14 <-   left_join(fitseq.data.residuals.below.14,fit.summary,by = c('lineage'))


fitseq.data.residuals.below.14 <- fitseq.data.residuals.below.14 %>% 
  mutate(underachiever = fit.resid< percentile)
    # mutate(underachiever = log2(freq.norm.anc.1)< -1.5)

above.n <- filter.n.fitness.residuals.underachiever(5,fitseq.data.residuals.below.14)
above.n <- above.n %>% mutate(n.pos.neg.under = ifelse(n.underachiever==T,'Underachiever',n.pos.neg))

variables <- 
  c('CAI','tAI','CDS.GC','dG','dG.noutr',	'TASEP.bottle.neck.position',	'TASEP.bottle.neck.depth','pep.cost','pep.boman','pep.hydro','pep.polar',
    'pep.charge','pep.basic','pep.acidic', 'sd.max.score','sd.mean', 'sd.sdev', 'sd.median', 'sd.count','log2(Trans)','log2(RNA)')
labels <- 
  c('CAI','tAI','Coding sequence GC %','Delta G','Delta G no UTR','Bottle neck position (calculated using TASEP)',
    'Bottle neck depth (calculated using TASEP)','Peptide cost','Peptide Boman index','Peptide hydrophobicity','Peptide polar amino acid content',
    'Peptide charge','Peptide basic amino acid content','Peptide acidic amino acid content',
    'Variable region max SD affinity score ','Variable region mean SD affinity score', 'Variable region mean SD affinity score stdev', 
    'Variable region median SD affinity score','Variable region SD affinity score count','Log 2 Translation efficeincy (protein/RNA)',
    'Log 2 RNA level ' )
names(labels) <- variables
result.dir <- 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\results\\underachievers\\'
day.number <-  12
data.set <- 'under.14.high.promoter'
legend.string <- 'n.pos.neg.under'
legend.label <-'Fitness residual'
for (i in 1:length(labels)){
  variable  <- names(labels)[i]
  label <- labels[i]
  print(label)
  compare.pos.neg.all.lineages.underachievers(above.n,variable,label,day.number,data.set,
                                              result.dir,legend.string,legend.label)
}
  

  
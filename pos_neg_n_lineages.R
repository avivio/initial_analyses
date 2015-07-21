load.packages()

get.wilcox.string.all.lineages  <-  function(legend.label,legend.string.true,legend.string.false){
  
  greater <- wilcox.test(legend.string.true,legend.string.false,alternative = 'greater')$p.value
  less <- wilcox.test(legend.string.true,legend.string.false,alternative = 'less')$p.value
  if(greater < less){
    return(paste('\n',legend.label ,'\n True > False\n p value =',formatC(greater,format='e',digits = 2)))
  } else{
    return(paste('\n',legend.label ,'\n True < False\n p value =',formatC(less,format='e',digits = 2)))
  }
  
  
}
#something tiny for commit

compare.pos.neg.all.lineages  <- function(fitseq.data.prediction,variable,label,day.number,data.set.name
                            ,result.dir,y.lim.all,y.lim.rbsm,legend.string,legend.label){
  
  fitseq.current.sample <-  fitseq.data.prediction %>% 
    mutate_(legend.string = legend.string)%>%
    filter(day == day.number)
  fitseq.current.sample.no.wt <-  fitseq.current.sample %>% 
    filter(RBS.Display != 'WT')
  
  
  legend.string.negative  <- fitseq.current.sample %>%
    ungroup %>% 
    filter(legend.string== 'Negative') %>%select_(variable) 
  legend.string.negative <- as.numeric( unlist(legend.string.negative))  
  legend.string.positive  <- fitseq.current.sample %>%
    ungroup %>% 
    filter(legend.string== 'Positive') %>%
    select_(variable) 
  legend.string.positive <- as.numeric( unlist(legend.string.positive)) 
  
  wilcox.string.all <- get.wilcox.string(legend.string.positive,legend.string.negative)
  
  
  fitseq.current.sample <- fitseq.current.sample %>%
    mutate(wilcox.all = wilcox.string.all)
  
  rbs.list = c('Strong','Mid','Weak')
  wilcox.string.rbs = vector(mode = 'character',length = 3)
  names(wilcox.string.rbs) <- rbs.list
  
  for (rbs in rbs.list){
    legend.string.negative.rbs  <- fitseq.current.sample.no.wt %>%
      ungroup %>% 
      filter(legend.string== 'Negative',RBS.Display == rbs) %>%
      select_(variable) 
    legend.string.negative.rbs <- as.numeric( unlist(legend.string.negative.rbs))  
    legend.string.positive.rbs  <- fitseq.current.sample.no.wt %>%
      ungroup %>% 
      filter(legend.string== 'Positive',RBS.Display == rbs) %>%
      select_(variable) 
    legend.string.positive.rbs <- as.numeric( unlist(legend.string.positive.rbs)) 
    
    wilcox.string.rbs[rbs] <- get.wilcox.string(legend.string.positive.rbs,legend.string.negative.rbs)
    
    
    
  }
  fitseq.current.sample.no.wt <- fitseq.current.sample.no.wt %>%
    mutate(wilcox.rbs = wilcox.string.rbs[RBS.Display])
  
  title <- paste('Comparing',label, 'between designs with and without',legend.label,
                 '\nDay', day.number, gsub('_',' ',data.set.name))
  
  png(paste0(result.dir,'all_pos_neg_hist_' ,variable,'_',data.set.name, '.png'),
      units="in",  width=15, height=12, res=70)
  p <- ggplot(fitseq.current.sample, aes_string(x = variable,fill =legend.string )) +
    geom_density(alpha = 0.4) +
    xlab(paste0(label,'\n'))+
    expand_limits(y = y.lim.all  ) +
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
    expand_limits(y = y.lim.rbs  ) +
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

col.plots.factor.all.lineages <- function(fitseq.xy,day.number,data.set.name,result.dir,
                             y.string,y.label,x.string,x.label,col.string,col.label){
  
  y.limits.all <- c(-5.5,4)
  y.limits.rbs.promoter <- c(-4.7,2)
  
  fitseq.current.sample <-  fitseq.xy %>% filter(day == day.number)
  fitseq.current.sample.no.wt <-  fitseq.xy %>% filter(day == day.number,RBS.Display != 'WT')
  
  
  filename.all  <- paste('all','day', day.number,'x',x.string ,'vs',y.string, 'col', col.string,data.set.name,sep = '_')
  filename.all  <-clean.filename(filename.all)
  title <- paste(x.label ,'vs',y.label,'\nfor day', day.number,gsub('_',' ',data.set.name))
  
  png(paste0(result.dir,filename.all,'.png'),units="in",  width=15, height=12, res=70)
  
  
  p <- ggplot(fitseq.current.sample,aes(x= x, y = y,color  = col)) +
    geom_point(alpha = 0.3)   +
    # xlim(x.limits) +
    coord_cartesian(ylim = y.limits.all) +
    ylab(paste0(y.label,'\n')) + 
    xlab(paste0('\n',x.label)) +
    ggtitle(title) + 
    geom_line(aes(x = x,y= prediction, color = col,inherit.aes=FALSE),
              size = 1.2,show_guide  = F)+
    # geom_abline(aes(slope = slope, intercept = intercept), color = 'dodgerblue', size = 1.2) +
    #     geom_smooth(aes_string(linetype = col.string),show_guide  = F,
    #                 color = 'red', size = 1.2, method = 'lm') +
    geom_text(aes(x=text.loc.x, y=text.loc.y, face="bold", inherit.aes=FALSE,
                  parse=FALSE,label=lm.string,colour=col,hjust=hor, vjust=ver),show_guide  = F) +
    guides(colour = guide_legend(override.aes = list(alpha = 1),title = col.label)) +
    theme_aviv
  print(p)
  
  dev.off()
  
  filename.promoter.rbs  <- paste('promoter_rbs','day', day.number,'x',x.string ,'vs',y.string, 'col',data.set.name, col.string,sep = '_')
  
  filename.promoter.rbs  <-clean.filename(filename.promoter.rbs)
  
  png(paste0(result.dir,filename.promoter.rbs,'.png'),units="in",  width=15, height=12, res=70)
  p <- ggplot(fitseq.current.sample.no.wt,aes(x= x, y = y,color  = col)) +
    geom_point(alpha = 0.3)   +
    # xlim(x.limits) +
    coord_cartesian(ylim = y.limits.rbs.promoter) +
    ylab(paste0(y.label,'\n')) + 
    xlab(paste0('\n',x.label)) +
    ggtitle(title) + 
    geom_line(aes(x = x,y= prediction, color = col, inherit.aes=FALSE),
              size = 1.2,show_guide  = F)+
    #     geom_smooth(
    #       aes_string(linetype = col.string),
    #       color = 'red', size = 1.2, method = 'lm',show_guide  = F) +
    facet_grid(Promoter.Display ~ RBS.Display)  +
    geom_text(aes(x=text.loc.x, y=text.loc.y, face="bold", inherit.aes=FALSE,
                  parse=FALSE,label=lm.string.group,colour=col,hjust=hor, vjust=ver),show_guide  = F) +
    guides(colour = guide_legend(override.aes = list(alpha = 1),title = col.label)) +
    theme_aviv
  
  print(p)
  dev.off()
}


filter.n.fitness.residuals <- function(n,fitseq.data.residuals){
  
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
  
  return(inner_join(fitseq.data.residuals,n.pos.neg))
  
  
}


fitseq.data.location  <- 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\clean_data\\goodman_salis_tuller_fitseq_2_mismatch_with_0.csv'
fitseq.data.tidy  <- load.fitseq.data(fitseq.data.location)
fitseq.data.tidy <- fitseq.data.tidy %>% select(-CDS.seq,-Promoter.seq,-RBS.seq,-Promoter,
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



fitseq.data.residuals <- fitseq.data.residuals %>% 
  filter( 
    Promoter.Display=='High'
    ,log2(Prot) < 17.5
    # ,abs(fit.resid) >  percentile
  )

fitseq.data.residuals <- filter.n.fitness.residuals(5,fitseq.data.residuals)

#filter above 14 and choose random lineage and drop lineage column because they don't matter anymore

fitseq.data.residuals.above.14 <- fitseq.data.residuals %>% filter(above.14 == TRUE) %>% 
  filter(lineage =='A') %>% ungroup() %>% select(-lineage)




data.set.name <- 'high_promoter_with_wall_5_positive'
legend.string <- 'n.pos.neg'
legend.label <- '5 lineages positive'

day.number <- 12


print(paste('fitseq data above 14',gsub('_',' ',data.set.name),'day',day.number))

fitseq.data.residuals.above.14 %>% 
  filter(day==day.number)%>%
  ungroup() %>%
  group_by(n.pos.neg) %>% 
  summarise(n())

print(paste('fitseq data above 14',gsub('_',' ',data.set.name),'day',day.number,'by rbs'))

fitseq.data.residuals.above.14  %>% 
  filter(day==day.number,RBS.Display != 'WT') %>%
  ungroup() %>%
  group_by(n.pos.neg,RBS.Display) %>% 
  summarise(n())



base.result.dir = 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\results\\fitness_residuals\\'
base.result.dir= paste0(base.result.dir,data.set.name,'\\')
dir.create(base.result.dir)
for (lineage.letter in c('A','B','C','D','E','F')){
  print(lineage.letter)
  
  #   for (day.number in seq(4,28,4)){
  #     print(day.number)
new.dir <- clean.filename(paste(x.string,col.string,sep = '_'))
  result.dir <- paste0(base.result.dir,new.dir,'\\')
  dir.create(result.dir)
  col.plots.factor(fitseq.data.residuals ,day.number,lineage.letter,data.set.name,
                   result.dir,y.string,y.label,x.string,x.label,col.string,col.label)
}
  
  result.dir <- paste0(base.result.dir,'pos_neg_histograms\\')
  dir.create(result.dir)
  
  
  variables <- 
    c('Rel.Codon.Freq','CAI','tAI','CDS.GC','GC','dG','dG.noutr','dG.unif','log2(salis.init)',
      'TASEP.avgRiboNum',	'TASEP.density.0',	'TASEP.bottle.neck.position',	'TASEP.bottle.neck.depth')
  labels <- 
    c('Relative codon frequency','CAI','tAI','Coding sequence GC %',' GC %','Delta G','Delta G no UTR',
      'Delta G cut at -5','Log 2 RBS calculator initition rate',
      'Average ribosome number (calculated using TASEP)','Density at start codon (calculated using TASEP)',
      'Bottle neck position (calculated using TASEP)','Bottle neck depth (calculated using TASEP)')
  names(labels) <- variables
  for (i in 1:length(labels)){
    variable  <- names(labels)[i]
    label <- labels[i]
    print(label)
    fitseq.data.residuals.above.14 <- fitseq.data.residuals.above.14 %>%
      mutate_(var = variable)
    dens.all <- fitseq.data.residuals.above.14 %>%
      filter(day == day.number) %>%
      group_by(n.pos.neg) %>%
      summarise(group.dens.all = max(unlist(density(var)[2]))) %>%
      ungroup() %>%
      summarise(max.dens.all = max(group.dens.all))
    y.lim.all <- as.numeric( unlist(dens.all))
    dens.rbs = vector(mode = 'numeric',length = 3)
    
    dens.rbs <- fitseq.data.residuals.above.14 %>%
      filter(day == day.number,RBS.Display != 'WT') %>%
      group_by(n.pos.neg,RBS.Display) %>%
      summarise(group.dens.all = max(unlist(density(var)[2]))) %>%
      ungroup() %>%
      summarise(max.dens.all = max(group.dens.all))
    y.lim.rbs <- max(dens.rbs)
    
    
    compare.pos.neg.all.lineages(fitseq.data.residuals.above.14,variable,label,day.number,data.set.name,
                    result.dir,y.lim.all,y.lim.rbs,legend.string,legend.label)
  }

# }




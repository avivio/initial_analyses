


#load packages and data
load.packages()




compare.pos.neg <- function(fitseq.data.prediction,variable,label,day.number,lineage.letter,data.set.name
                            ,result.dir,y.lim.all,y.lim.rbs){
  
  fitseq.current.sample <-  fitseq.data.prediction %>% 
    filter(day == day.number, lineage == lineage.letter)
  fitseq.current.sample.no.wt <-  fitseq.current.sample %>% 
    filter(RBS.Display != 'WT')
  
  
  neg  <- fitseq.current.sample %>%
    ungroup %>% 
    filter(pos.neg== 'Negative') %>%
    mutate_(new_variable = variable)  %>% select(new_variable)  
  neg <- as.numeric( unlist(neg))  
  pos  <- fitseq.current.sample %>%
    ungroup %>% 
    filter(pos.neg== 'Positive') %>%
    mutate_(new_variable = variable)  %>% select(new_variable)  
  pos <- as.numeric( unlist(pos)) 
  
  wilcox.string.all <- get.wilcox.string(pos,neg)
  
  
  fitseq.current.sample <- fitseq.current.sample %>%
    mutate(wilcox.all = wilcox.string.all)
  
  rbs.list = c('Strong','Mid','Weak')
  wilcox.string.rbs = vector(mode = 'character',length = 3)
  names(wilcox.string.rbs) <- rbs.list
  
  for (rbs in rbs.list){
    neg.rbs  <- fitseq.current.sample.no.wt %>%
      ungroup %>% 
      filter(pos.neg== 'Negative',RBS.Display == rbs) %>%
      mutate_(new_variable = variable)  %>% select(new_variable) 
    neg.rbs <- as.numeric( unlist(neg.rbs))  
    pos.rbs  <- fitseq.current.sample.no.wt %>%
      ungroup %>% 
      filter(pos.neg== 'Positive',RBS.Display == rbs) %>%
      mutate_(new_variable = variable)  %>% select(new_variable) 
    pos.rbs <- as.numeric( unlist(pos.rbs)) 
    
    wilcox.string.rbs[rbs] <- get.wilcox.string(pos.rbs,neg.rbs)
    
    
    
  }
  fitseq.current.sample.no.wt <- fitseq.current.sample.no.wt %>%
    mutate(wilcox.rbs = wilcox.string.rbs[RBS.Display])
  
  title <- paste('Comparing',label, 'between positive and negative fitness residuals',
                 '\nLineage', lineage.letter,'day', day.number, gsub('_',' ',data.set.name))
  
  png(paste0(result.dir,'all_pos_neg_hist_' ,variable,'_lineage_',lineage.letter,'_',data.set.name, '.png'),
      units="in",  width=15, height=12, res=70)
  p <- ggplot(fitseq.current.sample, aes_string(x = variable,fill = 'pos.neg')) +
    geom_density(alpha = 0.4) +
    xlab(paste0(label,'\n'))+
    expand_limits(y = y.lim.all  ) +
    ylab('\nDensity') +  
    ggtitle(title) +
    geom_text(aes( y = Inf, x = -Inf, label = wilcox.all), color="black",hjust = 0,vjust = 1, size = 6,
              parse=FALSE,fontface = 'bold',colour="red") +
    guides(fill=guide_legend(title='Fitness\nresidual')) +
    theme_aviv
  print(p)
  dev.off()
  
  png(paste0(result.dir,'promoter_rbs_pos_neg_hist_' ,variable,'_lineage_',lineage.letter,'_',data.set.name, '.png'),
      units="in",  width=15, height=12, res=70)
  p <- ggplot(fitseq.current.sample.no.wt, aes_string(x = variable,fill = 'pos.neg')) +
    geom_density(alpha = 0.4) +
    facet_grid(RBS.Display ~ .) +
    expand_limits(y = y.lim.rbs  ) +
    xlab(paste0(label,'\n')) +
    ylab('\nDensity') +  
    ggtitle(title) +
    geom_text(aes( y = Inf, x = -Inf, label = wilcox.rbs), color="black",hjust = 0,vjust = 1, size = 6,
              parse=FALSE,fontface = 'bold',colour="red") +
    guides(fill=guide_legend(title='Fitness\nresidual')) +
    theme_aviv
  print(p)
  dev.off()
  
  
  
}




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

col.plots.factor <- function(fitseq.xy,day.number,lineage.letter,data.set.name,result.dir,
                             y.string,y.label,x.string,x.label,col.string,col.label){
  
  y.limits.all <- c(-5.5,4)
  y.limits.rbs.promoter <- c(-4.7,2)
  
  fitseq.current.sample <-  fitseq.xy %>% filter(day == day.number, lineage == lineage.letter)
  fitseq.current.sample.no.wt <-  fitseq.xy %>% filter(day == day.number, lineage == lineage.letter, RBS.Display != 'WT')
  
  
  filename.all  <- paste('lineage', lineage.letter,'all','day', day.number,'x',x.string ,'vs',y.string, 'col', col.string,data.set.name,sep = '_')
  filename.all  <-clean.filename(filename.all)
  title <- paste(x.label ,'vs',y.label,'\nfor day', day.number, 'lineage', lineage.letter
                 ,gsub('_',' ',data.set.name))
  
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
  
  filename.promoter.rbs  <- paste('lineage', lineage.letter,'promoter_rbs','day', day.number,'x',x.string ,'vs',y.string, 'col',data.set.name, col.string,sep = '_')
  
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
  
  
  
  #decide what data set you're working on
  
  fit.summary  <- fitseq.data.residuals %>% 
    group_by(day,lineage,pos.neg) %>% 
    summarise(percentile= quantile(abs(fit.resid),probs = c(0.2))) 
  
  fitseq.data.residuals <-   left_join(fitseq.data.residuals,fit.summary,by = c('day','lineage','pos.neg'))
  # 
  # fitseq.data.residuals.wall.sample <- fitseq.data.residuals %>% 
  #   ungroup %>%
  #   filter( 
  #     Promoter.Display=='High',
  #     log2(Prot) > 17.5
  #     ,abs(fit.resid) >  percentile,
  #     lineage=='B',
  #     day == 12
  #   ) %>% sample_n(120) %>% select(Name)
  # 
  # fitseq.data.residuals.wall.sample <- inner_join(fitseq.data.residuals,fitseq.data.residuals.wall.sample,by= 'Name')
  
  
  
  # fitseq.data.residuals <-bind_rows(fitseq.data.residuals,fitseq.data.residuals.wall.sample)
  
  fitseq.data.residuals <- fitseq.data.residuals %>% 
    filter( 
      abs(fit.resid) >  percentile
    )
  
  # fitseq.data.residuals <- filter.n.fitness.residuals(6,'Positive',fitseq.data.residuals)
  
  
  fitseq.data.residuals.above.14 <- fitseq.data.residuals %>% filter(above.14 == TRUE)
  
  
  
  
  data.set.name <- 'high_promoter_no_wall_80_percentile_no_umi'
  
  
  
  lineage.letter = 'C'
  day.number = 12
  base.result.dir = 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\results\\no_umi_rerun\\fitness_residuals\\'
  result.dir= paste0(base.result.dir,data.set.name,'\\')
  dir.create(result.dir)
  
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
  
  
  
  # variables <- 
  #   c('pep.pi','pep.cost','pep.mw','pep.aindex','pep.boman','pep.charge','pep.hmoment','pep.hydro','pep.instability','pep.Ala','pep.Cys','pep.Asp','pep.Glu',
  #     'pep.Phe','pep.Gly','pep.His','pep.Ile','pep.Lys','pep.Leu','pep.Met','pep.Asn','pep.Pro','pep.Gln','pep.Arg','pep.Ser','pep.Thr','pep.Val','pep.Trp','pep.Tyr',
  #     'pep.tiny','pep.small','pep.aliphatic','pep.aromatic','pep.non.polar','pep.polar','pep.charged','pep.basic','pep.acidic')
  # labels <- 
  #   c('Peptide pI','Peptide cost','Peptide molecular weight','Peptide aliphatic index','Peptide Boman index','Peptide charge','Peptide hydrophyicl moment',
  #     'Peptide hydrophobicity','Peptide instability index','Peptide Alanine content','Peptide Cystine content','Peptide Aspartate content',
  #     'Peptide glutamine content','Peptide Phenylalanine content','Peptide Glycine content','Peptide histidine content','Peptide Isoleucine content',
  #     'Peptide Lysine content','Peptide Leucine content','Peptide Methionine content','Peptide Asparginine content','Peptide Proline content',
  #     'Peptide Glutamine content','Peptide Arginine content','Peptide Serine content','Peptide Threonine content','Peptide Valine content',
  #     'Peptide Tryptophan content','Peptide Tyrosine content','Peptide tiny amino acid content','Peptide small amino acid content',
  #     'Peptide aliphatic amino acid content','Peptide aromitc amino acid content','Peptide non polar amino acid content','Peptide polar amino acid content',
  #     'Peptide charged amino acid content','Peptide basic amino acid content','Peptide acidic amino acid content')
  # 
  #   variables <-     c('sd_rbs', 'sd_max_score', 'sd_max_position', 'sd_mean', 'sd_sdev', 'sd_median', 'sd_count', 'sd_gfp_max_score',
  #                      'sd_gfp_max_position', 'sd_gfp_mean', 'sd_gfp_sdev', 'sd_gfp_median', 'sd_gfp_count')
  #   
  #   labels <-     c('RBS SD affinity score', 'Variable region max SD affinity score ', 'Variable region max SD affinity score position',
  #                   'Variable region mean SD affinity score', 'Variable region mean SD affinity score stdev', 'Variable region median SD affinity score',
  #                   'Variable region SD affinity score count','GFP region max SD affinity score ', 'GFP region max SD affinity score position',
  #                   'GFP region mean SD affinity score', 'GFP region mean SD affinity score stdev', 'GFP region median SD affinity score',
  #                   'GFP region SD affinity score count')
  #   
    
  
  #   variables <- 
  #     c('Rel.Codon.Freq','CAI','tAI','CDS.GC','GC','dG','dG.noutr','dG.unif','log2(salis.init)',
  #       'TASEP.avgRiboNum',	'TASEP.density.0',	'TASEP.bottle.neck.position',	'TASEP.bottle.neck.depth')
  #   labels <- 
  #     c('Relative codon frequency','CAI','tAI','Coding sequence GC %',' GC %','Delta G','Delta G no UTR',
  #       'Delta G cut at -5','Log 2 RBS calculator initition rate',
  #       'Average ribosome number (calculated using TASEP)','Density at start codon (calculated using TASEP)',
  #       'Bottle neck position (calculated using TASEP)','Bottle neck depth (calculated using TASEP)')
    names(labels) <- variables
  
  
  for (lineage.letter in c('A','B','C','D','E','F')){
    print(lineage.letter)
    #   for (day.number in seq(4,28,4)){
    #     print(day.number)
    new.dir <- clean.filename(paste(x.string,col.string,sep = '_'))
    result.dir <- paste0(base.result.dir,new.dir,'\\')
    dir.create(result.dir)
    col.plots.factor(fitseq.data.residuals,day.number,lineage.letter,data.set.name,result.dir,
                     y.string,y.label,x.string,x.label,col.string,col.label)
    result.dir <- paste0(base.result.dir,'pos_neg_histograms\\')
    dir.create(result.dir)
    
    for (i in 1:length(labels)){
      variable  <- names(labels)[i]
      label <- labels[i]
      print(label)
      fitseq.data.residuals.above.14 <- fitseq.data.residuals.above.14 %>%
        mutate_(var = variable)
      dens.all <- fitseq.data.residuals.above.14 %>%
        filter(day == day.number) %>%
        group_by(pos.neg,lineage) %>%
        summarise(group.dens.all = max(unlist(density(var)[2]))) %>%
        ungroup() %>%
        summarise(max.dens.all = max(group.dens.all))
      y.lim.all <- as.numeric( unlist(dens.all))
      dens.rbs = vector(mode = 'numeric',length = 3)
      
      dens.rbs <- fitseq.data.residuals.above.14 %>%
        filter(day == day.number,RBS.Display != 'WT') %>%
        group_by(pos.neg,lineage,RBS.Display) %>%
        summarise(group.dens.all = max(unlist(density(var)[2]))) %>%
        ungroup() %>%
        summarise(max.dens.all = max(group.dens.all))
      y.lim.rbs <- max(dens.rbs)
      
      
      compare.pos.neg(fitseq.data.residuals.above.14,variable,label,day.number,lineage.letter,data.set.name,
                      result.dir,y.lim.all,y.lim.rbs)
    }
  }
  # }
  

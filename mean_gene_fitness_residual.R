
fitseq.data.location  <- 'C:\\Users\\shlomo\\Documents\\aviv\\workspace\\data\\pep_data_goodman_salis_tuller_fitseq_2_mismatch_with_0.csv'
fitseq.data.tidy  <- load.fitseq.data(fitseq.data.location)
# fitseq.data.tidy <- fitseq.data.tidy %>% select(-CDS.seq,-Promoter.seq,-RBS.seq,-Promoter,
#                                                 -variable.seq,-full.peptide,salis.status)

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


load.packages()

# get correlation data into data set
fitseq.data.residuals <- get.data.for.correlation.category(fitseq.data.tidy,y.string,y.label,x.string,x.label,col.string,col.label)
fit.summary  <- fitseq.data.residuals %>% 
  group_by(day,lineage,pos.neg) %>% 
  summarise(percentile= quantile(abs(fit.resid),probs = c(0.2))) 

fitseq.data.residuals <-   left_join(fitseq.data.residuals,fit.summary,by = c('day','lineage','pos.neg'))


fitseq.data.residuals <- fitseq.data.residuals %>% 
  filter( 
    Promoter.Display=='High')


fitseq.data.residuals <- fitseq.data.residuals %>% 
  filter( 
    log2(Prot) < 17.5,
    abs(fit.resid) >  percentile
  )



fitseq.data.residuals.above.14 <- fitseq.data.residuals %>% filter(above.14 == TRUE)








fitseq.data.residuals.gene  <- fitseq.data.residuals.above.14 %>% 
  group_by(day, lineage, Gene, Promoter.Display, RBS.Display) %>% 
  summarize(mean.resid = mean(fit.resid),
            median.resid = median(fit.resid),
            sd.resid = sd(fit.resid),
            sum.resid = sum(fit.resid)) 


fitseq.data.residual.rand  <- fitseq.data.residuals.above.14 %>% 
  group_by(day, lineage, Gene, Promoter.Display, RBS.Display) %>% 
  sample_n(1) %>% 
  select(-Name,-new.name)
fitseq.data.residuals.gene   <- left_join(fitseq.data.residual.rand ,
                                          fitseq.data.residuals.gene, 
                                          by = c('day','lineage','Gene','Promoter.Display','RBS.Display'))
fitseq.data.residuals.gene  <- mutate(fitseq.data.residuals.gene,
                                      pos.neg = ifelse(mean.resid > 0, 'Positive', 'Negative'))
data.set.name <- 'high_promoter_no_wall_80_percentile_gene_mean'

lineage.letter = 'C'
day.number = 12
base.result.dir = 'C:\\Users\\shlomo\\Documents\\aviv\\workspace\\results\\fitness_residuals_gene\\'
base.result.dir= paste0(base.result.dir,data.set.name,'\\')
dir.create(base.result.dir)

for (lineage.letter in c('A','B','C','D','E','F')){
  print(lineage.letter)
  result.dir <- paste0(base.result.dir,'pos_neg_histograms\\')
  dir.create(result.dir)
  
  
      variables <- 
        c('pep.pi','pep.cost','pep.mw','pep.aindex','pep.boman','pep.charge','pep.hmoment','pep.hydro','pep.instability','pep.Ala','pep.Cys','pep.Asp','pep.Glu',
          'pep.Phe','pep.Gly','pep.His','pep.Ile','pep.Lys','pep.Leu','pep.Met','pep.Asn','pep.Pro','pep.Gln','pep.Arg','pep.Ser','pep.Thr','pep.Val','pep.Trp','pep.Tyr',
          'pep.tiny','pep.small','pep.aliphatic','pep.aromatic','pep.non.polar','pep.polar','pep.charged','pep.basic','pep.acidic')
      labels <- 
        c('Peptide pI','Peptide cost','Peptide molecular weight','Peptide aliphatic index','Peptide Boman index','Peptide charge','Peptide hydrophobic moment',
          'Peptide hydrophobicity','Peptide instability index','Peptide Alanine content','Peptide Cystine content','Peptide Aspartate content',
          'Peptide glutamate content','Peptide Phenylalanine content','Peptide Glycine content','Peptide histidine content','Peptide Isoleucine content',
          'Peptide Lysine content','Peptide Leucine content','Peptide Methionine content','Peptide Asparginine content','Peptide Proline content',
          'Peptide Glutamine content','Peptide Arginine content','Peptide Serine content','Peptide Threonine content','Peptide Valine content',
          'Peptide Tryptophan content','Peptide Tyrosine content','Peptide tiny amino acid content','Peptide small amino acid content',
          'Peptide aliphatic amino acid content','Peptide aromitc amino acid content','Peptide non polar amino acid content','Peptide polar amino acid content',
          'Peptide charged amino acid content','Peptide basic amino acid content','Peptide acidic amino acid content')
      
      names(labels) <- variables
    
  
  for (i in 1:length(labels)){
    variable  <- names(labels)[i]
    label <- labels[i]
    print(label)
    fitseq.data.residuals.gene <- fitseq.data.residuals.gene %>%
      mutate_(var = variable)
    dens.all <- fitseq.data.residuals.gene %>%
      filter(day == day.number) %>%
      group_by(pos.neg,lineage) %>%
      summarise(group.dens.all = max(unlist(density(var)[2]))) %>%
      ungroup() %>%
      summarise(max.dens.all = max(group.dens.all))
    y.lim.all <- as.numeric( unlist(dens.all))
    dens.rbs = vector(mode = 'numeric',length = 3)
    
    dens.rbs <- fitseq.data.residuals.gene %>%
      filter(day == day.number,RBS.Display != 'WT') %>%
      group_by(pos.neg,lineage,RBS.Display) %>%
      summarise(group.dens.all = max(unlist(density(var)[2]))) %>%
      ungroup() %>%
      summarise(max.dens.all = max(group.dens.all))
    y.lim.rbs <- max(dens.rbs)
    
    
    compare.pos.neg(fitseq.data.residuals.gene,variable,label,day.number,lineage.letter,data.set.name,
                    result.dir,y.lim.all,y.lim.rbs)
  }
}
# }


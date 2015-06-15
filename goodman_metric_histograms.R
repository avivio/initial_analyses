result.dir = 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\results\\metric_histograms\\'

metrics <- c('Rel.Codon.Freq', 'log2(RNA)','Prot','log2(Trans)','CAI','tAI','CDS.GC','GC','dG','dG.noutr','dG.unif',
             'RNA.FCC',	'Prot.FCC',	'Trans.FCC',	'tAI.DC',	'GC.DC',	'CDS.GC.DC'	,'CAI.DC',
             'RSCU.DC',	'dG.DC',	'dG.noutr.DC',	'dG.unif.DC',	'log2(salis.init)','salis.dG_total',
             'salis.dG_mRNA_rRNA',	'salis.dG_mRNA','TASEP.avgRate',	'TASEP.avgRiboNum',	'TASEP.density.avg',
             'TASEP.density.0',	'TASEP.bottle.neck.position',	'TASEP.bottle.neck.depth')

x.labels <- 
  c('Relative codon frequency', 'Log 2 RNA level','Protein level','Log 2 translation efficiency','CAI','tAI'
    ,'Coding sequence GC %','GC %','Delta G','Delta G no UTR','Delta G from -5',
    'RNA level log fold change from family mean ','Protein level log fold change from family mean',	
    'Translational effeciency log fold change from family mean',	'tAI difference from family mean',	'GC % difference from family mean',
    'Coding sequence GC % difference from family mean'	,'CAI difference from family mean', 'RSCU difference from family mean',	'Delta G difference from family mean',	'Delta G no UTR difference from family mean',	'Delta G from -5 difference from family mean',
    'Log 2 RBS calculator initition rate','RBS calculator total delta G','RBS calculator delta G of mRNA rRNA interaction',
    'RBS calculator delta G mRNA','Average translation rate (calculated using TASEP)',	
    'Average ribosome number (calculated using TASEP)',	
    'Average ribosome density per codon (calculated using TASEP)',
    'Density at start codon (calculated using TASEP)',	'Bottle neck position (calculated using TASEP)',
    'Bottle neck depth (calculated using TASEP)')

names(x.labels) <- x.strings


for (i in 1:length(x.labels)){
  x.string <- names(x.labels)[i]
  x.label <- x.labels[i]
  
  
  print(x.string)
  filename  <- paste(x.string, 'historgram',sep = '_')
  filename  <-clean.filename(filename)
  
  png(paste0(result.dir,filename,'.png'),units="in",  width=15, height=12, res=70)
  
  # title <- paste(x.label, 'histogram')
  
  p <- ggplot(fitseq.data,aes_string(x.string)) + 
    geom_histogram() + 
    theme_aviv + 
    xlab(x.label)
  print(p)
  dev.off()
  
}

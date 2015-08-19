


library(scales)


for (lineage_letter in c('A','B','C','D','E','F')){
  

png(paste0('C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\results\\no_umi_rerun\\prot_vs_fitness_rbs\\fitness_prot_rbs',lineage_letter,'.png'),
    units="in",  width=15, height=12, res=70)

p <- ggplot(fitseq.data.residuals %>% filter(day == 12,lineage == lineage_letter ),
            aes(x=log2(Prot),y=log2(freq.norm.anc.1),color = as.logical(RBS.Display)) + geom_point(alpha = 0.3)   +    
  ylab(paste0('Frequency compared to ancestor frequency log 2','\n')) + 
  xlab(paste0('\n','Log 2 translation efficiency (protein/RNA)')) + 
  ggtitle(paste0('Fitness vs translation efficiency level day 12 lineage ',lineage_letter,'\nHigh promoter no wall')) +   
  coord_cartesian(ylim = c(-6,8)) +
  guides(colour = guide_legend(override.aes = list(alpha = 1),title = 'RBS')) +
  theme_aviv
  

print(p)

dev.off()
}
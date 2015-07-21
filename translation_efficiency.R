library(scales)


for (lineage.letter in c('A','B','C','D','E','F')){
  

png(paste0('C:\\Users\\shlomo\\Documents\\aviv\\workspace\\results\\fit_resid_fitness_trans_',lineage.letter,'.png'),
    units="in",  width=15, height=12, res=70)

p <- ggplot(fitseq.data.residuals.above.14 %>% filter(day == 12,lineage == lineage.letter ),
            aes(x=log2(Trans),y=log2(freq.norm.anc.1),color = pos.neg)) + geom_point(alpha = 0.3)   +    ylab(paste0('Frequency compared to ancestor frequency log 2','\n')) + 
  xlab(paste0('\n','Log 2 translation effeciency (protein/RNA)')) + 
  ggtitle(paste0('Fitness vs Translation effeciency day 12 lineage ',lineage.letter)) +   
            # geom_smooth(show_guide  = F,method = 'lm') +   
  coord_cartesian(ylim = c(-6,8)) +
  
            guides(colour = guide_legend(override.aes = list(alpha = 1),
                                         title = 'Fitness residual\npositive or negative')) +
  theme_aviv
print(p)
dev.off()
}
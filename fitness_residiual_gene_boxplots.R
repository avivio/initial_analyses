


load_packages()

fitseq_data <- load_db()


fitseq_data_residuals <- fitseq_data %>% 
  filter(
    day == 12,
    Promoter_Display=='High'
    ,above_14 == TRUE
    ,below_wall == TRUE
    # ,abs(fit_resid) >  percentile
  )



  


fitseq_data_residuals <- collect(fitseq_data_residuals)

fit_resid_limits <-c(min(fitseq.data.residuals.above.14$fit.resid), max(fitseq.data.residuals.above.14$fit.resid))

result_dir = 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\results\\no_umi_rerun\\gene_boxplots\\'
dir.create(result_dir)

gene_rank_matrix <- fitseq_data_residuals %>% 
  group_by(lineage, Gene) %>% 
  summarize(med_fit = median( fit_resid)) %>%
  group_by(lineage) %>%
  mutate(med_rank = min_rank(desc(med_fit))) %>% 
  arrange(Gene,med_rank) %>% 
  select(lineage,Gene,med_rank) %>% 
  spread(lineage,med_rank) %>% 
  select(-Gene)


gene_rank_matrix <- data.matrix(gene_rank_matrix)
print(kappam.fleiss(gene_rank_matrix))
print(icc(gene_rank_matrix))
print(kripp.alpha(gene_rank_matrix, method="ordinal"))
print(kendall(design_rank_matrix))


design_rank_matrix <- fitseq_data_residuals %>% 
  group_by(lineage,Name) %>% 
  summarize(med_fit = median( fit_resid)) %>%
  group_by(lineage) %>%
  mutate(med_rank = min_rank(desc(med_fit))) %>% 
  arrange(med_rank) %>% 
  select(lineage,Name,med_rank) %>% 
  spread(lineage,med_rank) %>% 
  select(-Name)


design_rank_matrix <- data.matrix(design_rank_matrix)
print(kappam.fleiss(design_rank_matrix))
print(icc(design_rank_matrix))
print(kripp.alpha(design_rank_matrix, method="ordinal"))
print(kendall(design_rank_matrix))

for (lineage_letter in c('A','B','C','D','E','F')){
  
  fitseq_current_data <- fitseq.data.residuals.above.14 %>% 
    filter(lineage == lineage_letter)

  sorted_genes <- fitseq_current_data %>% 
    group_by(lineage, Gene) %>% 
    summarize(med_fit = median( fit.resid)) %>%
    arrange(lineage,med_fit)
  
  
  fitseq_current_data$Gene <- factor(fitseq_current_data$Gene,
                                       levels = sorted_genes$Gene)
  

  filename  <- paste('fitness_gene_boxplot','lineage', lineage_letter,sep = '_')
  
png(paste0(result_dir,filename,'.png'),units="in",  width=25, height=12, res=70)
  
p <- ggplot(fitseq_current_data,
       aes(x = Gene, y = fit.resid)) +
  geom_boxplot() +
  expand_limits(y = fit_resid_limits) +
   # coord_cartesian(ylim = c(-1.5,6)) +
  theme_aviv +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab('Frequency normalized to ancestor frequency \n') +
  ggtitle(paste('Fitness by gene lineage',lineage_letter))
print(p)
dev.off()

}
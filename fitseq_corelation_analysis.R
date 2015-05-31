require(dplyr)
require(ineq)

results.dir = 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\results\\'

raw_data_location = 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\final_result_24-04_changed_label.csv'
fitseq_raw_data = read.csv(raw_data_location,check.names=FALSE)


by_evo_time_fitseq_data  <- fitseq_raw_data %>% 
  select(-design ,anc_1, anc_2, A_4, B_4,C_4,D_4,E_4,F_4, A_8, B_8,C_8,D_8,E_8,F_8
         , A_12, B_12,C_12,D_12,E_12,F_12, A_16, B_16,C_16,D_16,E_16,F_16, A_20, B_20,C_20,D_20,E_20,F_20
         , A_24, B_24,C_24,D_24,E_24,F_24, A_28, B_28,C_28,D_28,E_28,F_28)

by_line_fitseq_data  <- 
  fitseq_raw_data %>% 
  select(-design ,anc_1, anc_2, A_4,A_8,A_12,A_16,A_20,A_24,A_28,B_4,B_8,B_12,B_16,B_20,B_24,B_28,
         C_4,C_8,C_12,C_16,C_20,C_24,C_28,D_4,D_8,D_12,D_16,D_20,D_24,D_28,
         E_4,E_8,E_12,E_16,E_20,E_24,E_28,F_4,F_8,F_12,F_16,F_20,F_24,F_28)

fitseq_raw_data_no_design  <- fitseq_raw_data  %>% select(-design)
design_list = vector("list")
freq_list = vector("list")
fitseq_a_series_no_design <- fitseq_a_series  %>% select(-design)
for (i in 1:length(fitseq_a_series_no_design)){
  name  <- colnames(fitseq_a_series_no_design)[i]
  current_column  <- fitseq_a_series[,c('design',name)]
  ordered_current_column = current_column[order(current_column[,2],decreasing = T),]
  design_list = c(design_list,as.list(as.character(head(ordered_current_column,3)[,1]))))
freq_list = c(design_list,as.list(as.numeric(head(ordered_current_column,3)[,2])))
}
print 



fitseq_a_series  <- fitseq_raw_data %>% select(design,A_4,A_8,A_12,A_16,A_20,A_24,A_28)
fitseq_b_series  <- fitseq_raw_data %>% select(design,B_4,B_8,B_12,B_16,B_20,B_24,B_28)
fitseq_c_series  <- fitseq_raw_data %>% select(design,C_4,C_8,C_12,C_16,C_20,C_24,C_28)
fitseq_d_series  <- fitseq_raw_data %>% select(design,D_4,D_8,D_12,D_16,D_20,D_24,D_28)
fitseq_e_series  <- fitseq_raw_data %>% select(design,E_4,E_8,E_12,E_16,E_20,E_24,E_28)
fitseq_f_series  <- fitseq_raw_data %>% select(design,F_4,F_8,F_12,F_16,F_20,F_24,F_28)
fitseq_ancestors  <- fitseq_raw_data %>% select(design,anc_1, anc_2)

a_variants_over_time   <- 
  fitseq_a_series %>% select(-design)  %>%  summarise_each(funs(sum(.>0)/14234))

a_variants_over_time   <- 
  fitseq_a_series %>% select(-design)  %>%  summarise_each(funs(sum(.>0)/14234))

a_variants_over_time   <- 
  fitseq_a_series %>% select(-design)  %>%  summarise_each(funs(sum(.>0)/14234))

a_variants_over_time   <- 
  fitseq_a_series %>% select(-design)  %>%  summarise_each(funs(sum(.>0)/14234))

a_variants_over_time   <- 
  fitseq_a_series %>% select(-design)  %>%  summarise_each(funs(sum(.>0)/14234))

a_variants_over_time   <- 
  fitseq_a_series %>% select(-design)  %>%  summarise_each(funs(sum(.>0)/14234))



png(paste0(results.dir,'a_series_scatter_matrix.png'),units="in", width=20, height=20, res=100)

pairs(merge(fitseq_ancestors,fitseq_a_series,by='design')[,-1])
dev.off()


fitseq_tp_4  <- fitseq_raw_data %>% select(design,A_4,B_4,C_4,D_4,E_4,F_4)
fitseq_tp_8  <- fitseq_raw_data %>% select(design,A_8,B_8,C_8,D_8,E_8,F_8)
fitseq_tp_12  <- fitseq_raw_data %>% select(design,A_12,B_12,C_12,D_12,E_12,F_12)
fitseq_tp_16  <- fitseq_raw_data %>% select(design,A_16,B_16,C_16,D_16,E_16,F_16)
fitseq_tp_20  <- fitseq_raw_data %>% select(design,A_20,B_20,C_20,D_20,E_20,F_20)
fitseq_tp_24  <- fitseq_raw_data %>% select(design,A_24,B_24,C_24,D_24,E_24,F_24)
fitseq_tp_28  <- fitseq_raw_data %>% select(design,A_28,B_28,C_28,D_28,E_28,F_28)

#clustering
#spearman

fitseq_spearman_distance_matrix  <-  dist(1-cor(fitseq_raw_data %>% select(-design), method = "spearman"))

fitseq_spearman_clustering  <- hclust(fitseq_spearman_distance_matrix,method = 'average')

png(paste0(results.dir,'fitseq_spearman_upgma_clustering.png'),units="in", width=20, height=10, res=100)

plot(fitseq_spearman_clustering)

dev.off()


write.csv(cor(fitseq_raw_data %>% select(-design), method = "spearman"),
          paste0(results.dir,'spearman_correlation_matrix.csv'))

#pearson

fitseq_pearson_distance_matrix  <-  dist(1-cor(fitseq_raw_data %>% select(-design), method = "pearson"))

fitseq_pearson_clustering  <- hclust(fitseq_pearson_distance_matrix,method = 'average')

png(paste0(results.dir,'fitseq_pearson_upgma_clustering.png'),units="in", width=20, height=10, res=100)

plot(fitseq_pearson_clustering)

dev.off()
write.csv(cor(fitseq_raw_data %>% select(-design), method = "pearson"),
          paste0(results.dir,'pearson_correlation_matrix.csv'))



#scatter plots and matrices

library(lattice )
set = '28'

log10(fitseq_raw_data %>% select(contains(set)))


# scatter_plot_matrix_for  <-  function(set){

max =  max(fitseq_raw_data %>% select(contains(set)) )


png(paste0(results.dir,'scatter_plot_matrix_',set,'.png'),units="in", width=20, height=10, res=100)


splom(log10(fitseq_raw_data %>% select(contains(set)))
      #   fitseq_raw_data %>% select(contains(set))
      #       ,xlim =c(0,max),ylim =c(0,max)
      #      ,scales= list(
      #              relation = 'same',
      #                    log = TRUE,equispaced.log = TRUE
      
      # #                   ,at = list(max/1000,max/100,max/10,max/2,max)
      # #                     , labels = list(as.character(max/1000),as.character(max/100),as.character(max/10),
      # #                                  as.character(max/2),as.character(max)
      #                                   )
)
dev.off()
# }

max(fitseq_raw_data %>% select(anc_2,anc_1) )

anc_max = max(fitseq_raw_data %>% select(anc_2,anc_1) )
plot(unlist(fitseq_raw_data %>% select(anc_1)),
     unlist(fitseq_raw_data %>% select(anc_2)),
     log ='xy', xlim = c(1,anc_max),ylim = c(1,anc_max))
abline(0,1)
cor(fitseq_raw_data %>% select(anc_1),
    fitseq_raw_data %>% select(B_8))

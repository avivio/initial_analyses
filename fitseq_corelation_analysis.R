require(dplyr)
require(ineq)
require(dendextend)
results.dir = 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\results\\no_umi_rerun\\'

raw_data_location_2_mismatches = 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\match_count_2_mismatches_generations.csv'
fitseq_raw_data_2_mismatches = read.csv(raw_data_location_2_mismatches,check.names=FALSE)
colnames(fitseq_raw_data_2_mismatches)[1] = 'design'
fitseq_raw_data_mat_2_mismatches  <- data.matrix(fitseq_raw_data_2_mismatches %>% select(-design))
fitseq_raw_data_norm_2_mismatches  <- data.frame(sweep(fitseq_raw_data_mat_2_mismatches,2,colSums(fitseq_raw_data_mat_2_mismatches),`/`))
fitseq_raw_data_norm_log1p_2_mismatches  <- data.frame(log10(sweep(fitseq_raw_data_mat_2_mismatches+1,2,colSums(fitseq_raw_data_mat_2_mismatches+1),`/`)))



raw_data_location_0_mismatches = 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\clean_data\\final_result_14-05-15-1551_0_mismatches.csv'
fitseq_raw_data_0_mismatches = read.csv(raw_data_location_0_mismatches,check.names=FALSE)
colnames(fitseq_raw_data_0_mismatches)[1] = 'design'
fitseq_raw_data_0_mismatches  <- fitseq_raw_data_0_mismatches  %>% select(-anc_2)
fitseq_raw_data_mat_0_mismatches  <- data.matrix(fitseq_raw_data_0_mismatches %>% select(-design))
fitseq_raw_data_norm_0_mismatches  <- data.frame(sweep(fitseq_raw_data_mat_0_mismatches,2,colSums(fitseq_raw_data_mat_0_mismatches),`/`))
fitseq_raw_data_norm_log1p_0_mismatches  <- data.frame(log10(sweep(fitseq_raw_data_mat_0_mismatches+1,2,colSums(fitseq_raw_data_mat_0_mismatches+1),`/`)))


lineage_colors <- hue_pal()(7)
names(lineage_colors) <- c('Ancestor', 'A', 'B', 'C', 'D', 'E', 'F')

generation_colors <- hue_pal()(8)
names(generation_colors) <- seq(0,196,28)
generation_colors["0"] <- lineage_colors['Ancestor']

#clustering 2 mismatches
#spearman

fitseq_spearman_distance_matrix_2_mismatches  <-  dist(1-cor(fitseq_raw_data_2_mismatches %>% select(-design) , method = "spearman"))

fitseq_spearman_clustering_2_mismatches  <- hclust(fitseq_spearman_distance_matrix_2_mismatches,method = 'average')
spearman.dend <- as.dendrogram(fitseq_spearman_clustering_2_mismatches)

branch.col <- lineage_colors[as.character(t(as.data.frame(strsplit(names(fitseq_raw_data_2_mismatches)[-1][order.dendrogram(spearman.dend)],'_'))[1,]))]
spearman.dend <- color_branches(spearman.dend,col =branch.col )

labels_colors(spearman.dend) <- generation_colors[as.character(t(as.data.frame(strsplit(names(fitseq_raw_data_2_mismatches)[-1][order.dendrogram(spearman.dend)],'_'))[2,]))]

spearman.dend <- set(spearman.dend, 'branches_lwd',3)

png('C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\results\\for_thesis\\fitseq_spearman_upgma_clustering_2_mismatches.png',
        type="cairo",    units="in", width=14, height=6, pointsize=12, res=500)

plot(spearman.dend)

dev.off()


write.csv(cor(fitseq_raw_data %>% select(-design), method = "spearman"),
          paste0(results.dir,'spearman_correlation_matrix_2_mismatches.csv'))

#pearson

fitseq_pearson_distance_matrix_2_mismatches  <-  dist(1-cor(fitseq_raw_data_2_mismatches %>% select(-design), method = "pearson"))

fitseq_pearson_clustering_2_mismatches  <- hclust(fitseq_pearson_distance_matrix_2_mismatches,method = 'average')

pearson.dend <- as.dendrogram(fitseq_pearson_clustering_2_mismatches)

branch.col <- lineage_colors[as.character(t(as.data.frame(strsplit(names(fitseq_raw_data_2_mismatches)[-1][order.dendrogram(pearson.dend)],'_'))[1,]))]
pearson.dend <- color_branches(pearson.dend,col =branch.col )

labels_colors(pearson.dend) <- generation_colors[as.character(t(as.data.frame(strsplit(names(fitseq_raw_data_2_mismatches)[-1][order.dendrogram(pearson.dend)],'_'))[2,]))]
pearson.dend <- set(pearson.dend, 'branches_lwd',3)

png('C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\results\\for_thesis\\fitseq_pearson_upgma_clustering_2_mismatches.png',
    type="cairo",    units="in", width=14, height=6, pointsize=12, res=500)

plot(pearson.dend)

dev.off()

write.csv(cor(fitseq_raw_data %>% select(-design), method = "pearson"),
          paste0(results.dir,'pearson_correlation_matrix_2_mismatches.csv'))




#clustering 0 mismatches
#spearman

fitseq_spearman_distance_matrix_0_mismatches  <-  dist(1-cor(fitseq_raw_data_0_mismatches %>% select(-design), method = "spearman"))

fitseq_spearman_clustering_0_mismatches  <- hclust(fitseq_spearman_distance_matrix_0_mismatches,method = 'average')

png(paste0(results.dir,'fitseq_spearman_upgma_clustering_0_mismatches.png'),units="in", width=20, height=10, res=100)

plot(fitseq_spearman_clustering,main = 'Spearman correlation UPGMA dendrogram based on 0 mismatch data')

dev.off()


write.csv(cor(fitseq_raw_data %>% select(-design), method = "spearman"),
          paste0(results.dir,'spearman_correlation_matrix_0_mismatches.csv'))

#pearson

fitseq_pearson_distance_matrix_0_mismatches  <-  dist(1-cor(fitseq_raw_data_0_mismatches %>% select(-design), method = "pearson"))

fitseq_pearson_clustering_0_mismatches  <- hclust(fitseq_pearson_distance_matrix_0_mismatches,method = 'average')

png(paste0(results.dir,'fitseq_pearson_upgma_clustering_0_mismatches.png'),units="in", width=20, height=10, res=100)

plot(fitseq_pearson_clustering,main = 'Pearson correlation UPGMA dendrogram based on 0 mismatch data')

dev.off()

write.csv(cor(fitseq_raw_data %>% select(-design), method = "pearson"),
          paste0(results.dir,'pearson_correlation_matrix_0_mismatches.csv'))




png(paste0(results.dir,'scatter_plot_matrix_lineage_C.png'),units="in", width=20, height=10, res=100)


  ggpairs_identical_scales(select(fitseq_raw_data_norm_log1p_2_mismatches,C_4,C_8,C_12,C_16,C_20,C_24,C_28))


dev.off()


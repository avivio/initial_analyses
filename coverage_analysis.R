

require(dplyr)
require(ggplot2)
require(scales)
require(tidyr)

  summary_data_location = 'D:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\final_summary_08-05-15-2225.csv'
summary_data = t(read.csv(summary_data_location,check.names=DALSE,row.names = 1,as.is = T))
summary_data <- transform(summary_data,
                          all_reads = as.numeric(as.character(all_reads)),
                          merged_reads= as.numeric(as.character(merged_reads)),
                          trimmed_reads= as.numeric(as.character(trimmed_reads)),
                          found_designs_0_mismatches = as.numeric(as.character(found_designs_0_mismatches)),
                          design_count_pre_umi_0_mismatches = as.numeric(as.character(design_count_pre_umi_0_mismatches)),
                          design_count_post_umi_0_mismatches = as.numeric(as.character(design_count_post_umi_0_mismatches)),
                          found_designs_1_mismatches = as.numeric(as.character(found_designs_1_mismatches)),
                          design_count_pre_umi_1_mismatches = as.numeric(as.character(design_count_pre_umi_1_mismatches)), 
                          design_count_post_umi_1_mismatches = as.numeric(as.character(design_count_post_umi_1_mismatches)),
                          found_designs_2_mismatches = as.numeric(as.character(found_designs_2_mismatches)),   
                          design_count_pre_umi_2_mismatches = as.numeric(as.character(design_count_pre_umi_2_mismatches)),
                          design_count_post_umi_2_mismatches = as.numeric(as.character(design_count_post_umi_2_mismatches)),
                          found_designs_3_mismatches = as.numeric(as.character(found_designs_3_mismatches)),
                          design_count_pre_umi_3_mismatches = as.numeric(as.character(design_count_pre_umi_3_mismatches)), 
                          design_count_post_umi_3_mismatches = as.numeric(as.character(design_count_post_umi_3_mismatches)),
                          found_designs_all = as.numeric(as.character(found_designs_all)),   
                          design_count_pre_umi_all = as.numeric(as.character(design_count_pre_umi_all)),
                          design_count_post_umi_all = as.numeric(as.character(design_count_post_umi_all)),
                          line = as.character(line),
                          day = as.numeric(as.character(day))
)




results_dir = 'D:\\Users\\dell7\\Documents\\Tzachi\\workspace\\results\\'


fitseq_data_location = 'D:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\final_result_08-05-15-2225_all_mismatches.csv'
fitseq_data = read.csv(fitseq_data_location,row.names = 1)

n_terminal_data_location = 'D:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\goodman_salis_tuller_new_name_data.csv'
n_terminal_data = read.csv(n_terminal_data_location,row.names = 1)

fitseq_n_term_data = merge(fitseq_data,n_terminal_data,by = 0)


breaks = c(0,100,250,500,1000,2000,5000,10000,50000,100000,200000)
max_breaks = max(breaks)

p_anc = ggplot(fitseq_n_term_data, aes(x=anc_1/sum(anc_1))) + geom_histogram() +  facet_wrap(~Promoter.Display)
p_anc = p_anc  +  scale_y_continuous(trans=log1p_trans(),breaks = lseq(1,4e-04,10)) +  scale_x_continuous(trans=log1p_trans(),breaks = breaks) + expand_limits(y=2000,x=max_breaks)
p_anc = p_anc + xlab('log(frequency+1)') + ylab('log(count + 1)') + ggtitle('Ancestor 1 log frequency histogram')

p_D_16 = ggplot(fitseq_n_term_data, aes(x=D_16)) + geom_histogram() +  facet_wrap(~Promoter.Display)  
p_D_16 = p_D_16 +  scale_y_continuous(trans=log1p_trans(),breaks = c(0,500,1000,2000)) +  scale_x_continuous(trans=log1p_trans(),breaks = breaks) + expand_limits(y=2000,x=max_breaks)
p_D_16 = p_D_16 + xlab('log(frequency+1)') + ylab('log(count + 1)') + ggtitle('Line D day 16  log frequency histogram')


p_D_28 = ggplot(fitseq_n_term_data, aes(x=D_28)) + geom_histogram() +  facet_wrap(~Promoter.Display)
p_D_28 = p_D_28 + scale_y_continuous(trans=log1p_trans(),breaks = c(0,500,1000,2000)) +  scale_x_continuous(trans=log1p_trans(),breaks = breaks) + expand_limits(y=2000,x=max_breaks)
p_D_28 = p_D_28 + xlab('log(frequency+1)') + ylab('log(count + 1)') + ggtitle('Line D day 28  log frequency histogram')


multiplot(p_anc,p_D_16, p_D_28)


p_anc_2 = ggplot(fitseq_n_term_data, aes(x=anc_2)) + geom_histogram() +  facet_wrap(~Promoter.Display)
p_anc_2 = p_anc_2  +  scale_y_continuous(trans=log1p_trans(),breaks = c(0,500,1000,2000)) +  scale_x_continuous(trans=log1p_trans(),breaks = breaks) + expand_limits(y=2000,x=max_breaks)
p_anc_2 = p_anc_2 + xlab('log(frequency+1)') + ylab('log(count + 1)') + ggtitle('Ancestor 2 log frequency histogram')

p_anc_1 = ggplot(fitseq_n_term_data, aes(x=anc_1)) + geom_histogram() +  facet_wrap(~Promoter.Display)
p_anc_1 = p_anc_1  +  scale_y_continuous(trans=log1p_trans(),breaks = c(0,500,1000,2000)) +  scale_x_continuous(trans=log1p_trans(),breaks = breaks) + expand_limits(y=2000,x=max_breaks)
p_anc_1 = p_anc_1 + xlab('log(frequency+1)') + ylab('log(count + 1)') + ggtitle('Ancestor 1 log frequency histogram')

multiplot(p_anc_1,p_anc_2)

#sqrt
ggplot(fitseq_n_term_data, aes(x=log1p(anc_1))) + geom_histogram() +  facet_wrap(~Promoter.Display) + scale_y_sqrt()




ggplot(fitseq_n_term_data,aes(x = c('anc_1','anc_2','A_4','A_8','A_12','A_16','A_20',,'A_24','A_28'))) + facet_wrap(~Promoter.Display) + geom_boxplot()
ggplot(fitseq_n_term_data, aes(x=anc_2)) + geom_boxplot() +  facet_grid(RBS.Display~Promoter.Display)
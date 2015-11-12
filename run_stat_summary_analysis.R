
require(dplyr)
require(ggplot2)
require(ggvis)
require(scales)
require(tidyr)

results.dir = 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\results\\no_umi_rerun'

summary_data_location = 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\final_summary_generations.csv'
summary_data = t(read.csv(summary_data_location,check.names=FALSE,row.names = 1,as.is = T))
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
          found_designs_4_mismatches = as.numeric(as.character(found_designs_4_mismatches)),
          design_count_pre_umi_4_mismatches = as.numeric(as.character(design_count_pre_umi_4_mismatches)), 
          design_count_post_umi_4_mismatches = as.numeric(as.character(design_count_post_umi_4_mismatches)),          
          found_designs_all = as.numeric(as.character(found_designs_all)),   
          design_count_pre_umi_all = as.numeric(as.character(design_count_pre_umi_all)),
          design_count_post_umi_all = as.numeric(as.character(design_count_post_umi_all)),
          lineage = as.character(lineage),
          generation = as.numeric(as.character(generation))
)
          

# total reads over days 
summary_data %>% ggvis(~day,~all_reads, fill = ~line) %>% layer_points() %>% add_axis( "x", 
title = "Days") %>% add_axis( "y", title = "reads", title_offset = 75)


#  merge percent over days 
summary_data %>%mutate(merged_perc  <- merged_reads/all_reads )  %>% ggvis(
  ~day,~merged_perc, fill = ~lineage) %>% layer_points() %>% add_axis(
    "x", title = "Days") %>% add_axis( "y", title = "merged reads over all reads", title_offset = 75)

#  trim percent from all over days 
summary_data %>%mutate(trim_all_perc =trimmed_reads/all_reads )  %>% ggvis(
  ~day,~trim_all_perc, fill = ~line) %>% layer_points() %>% add_axis( 
    "x",title = "Days") %>% add_axis( "y", title = "trimmed reads over all reads", title_offset = 75)

#  trim percent from merged over days 
summary_data %>%mutate(trim_merge_perc = trimmed_reads/merged_reads )  %>% ggvis(
  ~day,~trim_merge_perc, fill = ~line) %>% layer_points() %>% add_axis( 
    "x",title = "Days") %>% add_axis( "y", title = "trimmed reads over merged reads", title_offset = 75)

#  found designs percent over days 0 mismatches
summary_data %>% mutate(found_design_perc = found_designs_0_mismatches/14234 )  %>%  ggvis(
  ~day,~found_design_perc, fill = ~line) %>% layer_points() %>% add_axis( 
    "x",title = "Days") %>% add_axis(
      "y", title = "found designs percent from library 0 mismatches", title_offset = 75)

#  found designs percent over days 1 mismatches
summary_data %>% mutate(found_design_perc = found_designs_1_mismatches/14234 )  %>%  ggvis(
  ~day,~found_design_perc, fill = ~line) %>% layer_points() %>% add_axis( 
    "x",title = "Days") %>% add_axis(
      "y", title = "found designs percent from library 1 mismatches", title_offset = 75)

#  found designs percent over days 2 mismatches

 summary_data %>% mutate(found_design_perc = found_designs_2_mismatches/14234 )  %>%  
  ggvis(~day,~found_design_perc, fill = ~lineage) %>% 
  layer_points() %>% 
  add_axis( "x",title = "Days") %>%
  add_axis("y", title = "found designs percent from library 2 mismatches)", title_offset = 75)

 summary_data_design_coverage_2_mismatches <- summary_data %>% mutate(found_design_perc = found_designs_2_mismatches/14234 )  
 summary_data_design_coverage_2_mismatches$lineage <- factor(summary_data_design_coverage_2_mismatches$lineage,
                                    levels = c("Ancestor", "A", "B", "C", "D", "E", "F"))
 
png('C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\results\\for_thesis\\percent_of_library_covered_over_time.png',
    type="cairo",    units="in", width=10, height=6, pointsize=12, res=500)    
ggplot(summary_data_design_coverage_2_mismatches, aes(x = generation,y = found_design_perc, color = lineage)) + 
  geom_point(size = 3) +
  xlab("Generations") +
  ylab("Percent of library covered") + 
  scale_x_continuous(breaks=c(0, 28, 56, 84, 112,140,168,196)) +
  theme_aviv

dev.off()

#  found designs percent over days 3 mismatches
summary_data %>% mutate(found_design_perc = found_designs_3_mismatches/14234 )  %>%  ggvis(
  ~day,~found_design_perc, fill = ~line) %>% layer_points() %>% add_axis( 
    "x",title = "Days") %>% add_axis(
      "y", title = "found designs percent from library 3 mismatches)", title_offset = 75)

#  found designs percent over days 4 mismatches
summary_data %>% mutate(found_design_perc = found_designs_4_mismatches/14234 )  %>%  ggvis(
  ~day,~found_design_perc, fill = ~line) %>% layer_points() %>% add_axis( 
    "x",title = "Days") %>% add_axis(
      "y", title = "found designs percent from library 4 mismatches)", title_offset = 75)


#  match percent from all over days 0 mismatch
summary_data %>%mutate(match_percent =design_count_pre_umi_0_mismatches/all_reads )  %>% ggvis(
  ~day,~match_percent, fill = ~line) %>% layer_points() %>% add_axis( 
    "x",title = "Days") %>% add_axis( "y", title = "matched reads over all reads 0 mismatches", title_offset = 75)

#  match percent from trimmed reads over days 0 mismatch
summary_data %>%mutate(match_percent =design_count_pre_umi_0_mismatches/trimmed_reads )  %>% ggvis(
  ~day,~match_percent, fill = ~line) %>% layer_points() %>% add_axis( 
    "x",title = "Days") %>% add_axis( "y", title = "matched reads over trimmed reads 0 mismatches", title_offset = 75)



#  match percent from all over days 1 mismatch
summary_data %>%mutate(match_percent =design_count_pre_umi_1_mismatches/all_reads )  %>% ggvis(
  ~day,~match_percent, fill = ~line) %>% layer_points() %>% add_axis( 
    "x",title = "Days") %>% add_axis( "y", title = "matched reads over all reads 1 mismatches", title_offset = 75)

#  match percent from trimmed reads over days 1 mismatch
summary_data %>%mutate(match_percent =design_count_pre_umi_1_mismatches/trimmed_reads )  %>% ggvis(
  ~day,~match_percent, fill = ~line) %>% layer_points() %>% add_axis( 
    "x",title = "Days") %>% add_axis( "y", title = "matched reads over trimmed reads 1 mismatches", title_offset = 75)


#  match percent from all over days 2 mismatch
summary_data %>%mutate(match_percent =design_count_pre_umi_2_mismatches/all_reads )  %>% ggvis(
  ~day,~match_percent, fill = ~line) %>% layer_points() %>% add_axis( 
    "x",title = "Days") %>% add_axis( "y", title = "matched reads over all reads 2 mismatches", title_offset = 75)

#  match percent from trimmed reads over days 2 mismatch
summary_data %>%mutate(match_percent =design_count_pre_umi_2_mismatches/trimmed_reads )  %>% ggvis(
  ~day,~match_percent, fill = ~line) %>% layer_points() %>% add_axis( 
    "x",title = "Days") %>% add_axis( "y", title = "matched reads over trimmed reads 2 mismatches", title_offset = 75)


#  match percent from all over days 3 mismatch
summary_data %>%mutate(match_percent =design_count_pre_umi_3_mismatches/all_reads )  %>% ggvis(
  ~day,~match_percent, fill = ~line) %>% layer_points() %>% add_axis( 
    "x",title = "Days") %>% add_axis( "y", title = "matched reads over all reads 3 mismatches", title_offset = 75)

#  match percent from trimmed reads over days 3 mismatch
summary_data %>%mutate(match_percent =design_count_pre_umi_3_mismatches/trimmed_reads )  %>% ggvis(
  ~day,~match_percent, fill = ~line) %>% layer_points() %>% add_axis( 
    "x",title = "Days") %>% add_axis( "y", title = "matched reads over trimmed reads 3 mismatches", title_offset = 75)


#  match percent from all over days 4 mismatch
summary_data %>%mutate(match_percent =design_count_pre_umi_4_mismatches/all_reads )  %>% ggvis(
  ~day,~match_percent, fill = ~line) %>% layer_points() %>% add_axis( 
    "x",title = "Days") %>% add_axis( "y", title = "matched reads over all reads 4 mismatches", title_offset = 75)

#  match percent from trimmed reads over days 4 mismatch
summary_data %>%mutate(match_percent =design_count_pre_umi_4_mismatches/trimmed_reads )  %>% ggvis(
  ~day,~match_percent, fill = ~line) %>% layer_points() %>% add_axis( 
    "x",title = "Days") %>% add_axis( "y", title = "matched reads over trimmed reads 4 mismatches", title_offset = 75)

#  umi frequency from all over days 0 mismatch
summary_data %>%mutate(umi_percent =design_count_post_umi_0_mismatches/all_reads )  %>% ggvis(
  ~day,~umi_percent, fill = ~line) %>% layer_points() %>% add_axis( 
    "x",title = "Days") %>% add_axis( "y", title = "umi frequency over all reads 0 mismatches", title_offset = 75)

#  umi frequency from trimmed reads over days 0 mismatch
summary_data %>%mutate(umi_percent =design_count_post_umi_0_mismatches/trimmed_reads )  %>% ggvis(
  ~day,~umi_percent, fill = ~line) %>% layer_points() %>% add_axis( 
    "x",title = "Days") %>% add_axis( "y", title = "umi frequency over trimmed reads 0 mismatches", title_offset = 75)

#  umi frequency from all over days 1 mismatch
summary_data %>%mutate(umi_percent =design_count_post_umi_all/all_reads )  %>% ggvis(
  ~day,~umi_percent, fill = ~line) %>% layer_points() %>% add_axis( 
    "x",title = "Days") %>% add_axis( "y", title = "umi frequency over all reads 3 mismatches", title_offset = 75)

#  umi frequency from trimmed reads over days 1 mismatch
summary_data %>%mutate(umi_percent =design_count_post_umi_1_mismatches/trimmed_reads )  %>% ggvis(
  ~day,~umi_percent, fill = ~line) %>% layer_points() %>% add_axis( 
    "x",title = "Days") %>% add_axis( "y", title = "umi frequency over trimmed reads 1 mismatches", title_offset = 75)

#  umi frequency from all over days 2 mismatch
summary_data %>%mutate(umi_percent =design_count_post_umi_2_mismatches/all_reads )  %>% ggvis(
  ~day,~umi_percent, fill = ~line) %>% layer_points() %>% add_axis( 
    "x",title = "Days") %>% add_axis( "y", title = "umi frequency over all reads 2 mismatches", title_offset = 75)

#  umi frequency from trimmed reads over days 2 mismatch
summary_data %>%mutate(umi_percent =design_count_post_umi_2_mismatches/trimmed_reads )  %>% ggvis(
  ~day,~umi_percent, fill = ~line) %>% layer_points() %>% add_axis( 
    "x",title = "Days") %>% add_axis( "y", title = "umi frequency over trimmed reads 2 mismatches", title_offset = 75)

#  umi frequency from all over days 3 mismatch
summary_data %>%mutate(umi_percent =design_count_post_umi_3_mismatches/all_reads )  %>% ggvis(
  ~day,~umi_percent, fill = ~line) %>% layer_points() %>% add_axis( 
    "x",title = "Days") %>% add_axis( "y", title = "umi frequency over all reads 3 mismatches", title_offset = 75)

#  umi frequency from trimmed reads over days 3 mismatch
summary_data %>%mutate(umi_percent =design_count_post_umi_3_mismatches/trimmed_reads )  %>% ggvis(
  ~day,~umi_percent, fill = ~line) %>% layer_points() %>% add_axis( 
    "x",title = "Days") %>% add_axis( "y", title = "umi frequency over trimmed reads 3 mismatches", title_offset = 75)

#  umi frequency from all over days 4 mismatch
summary_data %>%mutate(umi_percent =design_count_post_umi_4_mismatches/all_reads )  %>% ggvis(
  ~day,~umi_percent, fill = ~line) %>% layer_points() %>% add_axis( 
    "x",title = "Days") %>% add_axis( "y", title = "umi frequency over all reads 4 mismatches", title_offset = 75)

#  umi frequency from trimmed reads over days 4 mismatch
summary_data %>%mutate(umi_percent =design_count_post_umi_4_mismatches/trimmed_reads )  %>% ggvis(
  ~day,~umi_percent, fill = ~line) %>% layer_points() %>% add_axis( 
    "x",title = "Days") %>% add_axis( "y", title = "umi frequency over trimmed reads 4 mismatches", title_offset = 75)

lineage_letter = 'A'
for (lineage_letter in c('A','B','C','D','E','F')){ 
print(lineage_letter)

#match count 0,1,2,3,4 mismatches over days for a lineage
match  <- summary_data %>% 
  filter(grepl(lineage_letter,lineage)) %>% 
  mutate(match_perc_0 =design_count_pre_umi_0_mismatches/trimmed_reads ,
         match_perc_1 =design_count_pre_umi_1_mismatches/trimmed_reads,
         match_perc_2 =design_count_pre_umi_2_mismatches/trimmed_reads,
         match_perc_3 =design_count_pre_umi_3_mismatches/trimmed_reads,
         match_perc_4 =design_count_pre_umi_4_mismatches/trimmed_reads)%>%
  select(day,match_perc_0,match_perc_1,match_perc_2,match_perc_3,match_perc_4)

# png(paste0(results.dir,lineage_letter,'_mismatch_comparison_pre_umi.png'),units="in", width=30, height=15, res=50)


 p <- ggplot(match, aes(day)) +
  geom_line(aes(y = match_perc_0, colour = "0"),size =1.5) +
  geom_line(aes(y = match_perc_1, colour = "1"),size =1.5) +
  geom_line(aes(y = match_perc_2, colour = "2"),size =1.5) +
  geom_line(aes(y = match_perc_3, colour = "3"),size =1.5) +
  geom_line(aes(y = match_perc_4, colour = "4"),size =1.5) +
  geom_point(aes(y = match_perc_0, colour = "0"),size =5) +
  geom_point(aes(y = match_perc_1, colour = "1"),size =5) +
  geom_point(aes(y = match_perc_2, colour = "2"),size =5) +
  geom_point(aes(y = match_perc_3, colour = "3"),size =5) +
  geom_point(aes(y = match_perc_4, colour = "4"),size =5) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold"),
                     legend.text=element_text(size=14), legend.title=element_text(size=16,face="bold"),
                     plot.title = element_text(size = 25,face = "bold")) +
  ylab('Percent of matched reads out of trimmed reads\n') +
  guides(colour = guide_legend("Allowed mismatches")) + expand_limits(y=c(0,1))  + 
  ggtitle(paste0('lineage ',lineage_letter,' coverage of trimmed reads\nfor different mismatch values over time')) +
   theme_minimal() + 
   theme( axis.line = element_line(colour = "black"),
          axis.text=element_text(size=18), axis.title=element_text(size=20,face="bold"),
          strip.text.x = element_text(size=20,face="bold"), strip.text.y = element_text(size=20,face="bold"), 
          plot.title = element_text(size = 30,face = "bold"),
          legend.text = element_text(size=20),legend.title = element_text(size=22,face="bold")	)
print(p) 
# dev.off()
}

#umi count 0,1,2,3,4 mismatches over days for a lineage
umi  <- summary_data %>% 
  filter(grepl(lineage_letter,lineage)) %>% 
  mutate(umi_perc_0 =design_count_post_umi_0_mismatches/trimmed_reads ,
         umi_perc_1 =design_count_post_umi_1_mismatches/trimmed_reads,
         umi_perc_2 =design_count_post_umi_2_mismatches/trimmed_reads,
         umi_perc_3 =design_count_post_umi_3_mismatches/trimmed_reads,
         umi_perc_4 =design_count_post_umi_4_mismatches/trimmed_reads)%>%
  select(day,umi_perc_0,umi_perc_1,umi_perc_2,umi_perc_3,umi_perc_4)

png(paste0(results.dir,lineage_letter,'_mismatch_comparison_post_umi.png'),units="in", width=30, height=15, res=50)

ggplot(umi, aes(day)) +
  geom_line(aes(y = umi_perc_0, colour = "0"),size =1.5) +
  geom_line(aes(y = umi_perc_1, colour = "1"),size =1.5) +
  geom_line(aes(y = umi_perc_2, colour = "2"),size =1.5) +
  geom_line(aes(y = umi_perc_3, colour = "3"),size =1.5) +
  geom_line(aes(y = umi_perc_4, colour = "4"),size =1.5) +
  geom_point(aes(y = umi_perc_0, colour = "0"),size =5) +
  geom_point(aes(y = umi_perc_1, colour = "1"),size =5) +
  geom_point(aes(y = umi_perc_2, colour = "2"),size =5) +
  geom_point(aes(y = umi_perc_3, colour = "3"),size =5) +
  geom_point(aes(y = umi_perc_4, colour = "4"),size =5) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold"),
                     legend.text=element_text(size=14), legend.title=element_text(size=16,face="bold"),
                     plot.title = element_text(size = 25,face = "bold")) +
  ylab('Percent of umi filtered matched reads out of trimmed reads\n') +
  guides(colour = guide_legend("Allowed mismatches")) +
  expand_limits(y=c(0,1)) + 
  ggtitle(paste0('lineage ',lineage_letter,' post UMI coverage of trimmed reads for different mismatch values over time')) +
  theme_minimal() + 
  theme( axis.line = element_line(colour = "black"),
         axis.text=element_text(size=18), axis.title=element_text(size=20,face="bold"),
         strip.text.x = element_text(size=20,face="bold"), strip.text.y = element_text(size=20,face="bold"), 
         plot.title = element_text(size = 30,face = "bold"),
         legend.text = element_text(size=20),legend.title = element_text(size=22,face="bold")	)

dev.off()


ggplot(sum, aes(x=day,y=merged_perc,colour=lineage)) +geom_point(size = 5)+
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                     axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold"),
                     legend.text=element_text(size=14), legend.title=element_text(size=16,face="bold"),
                     plot.title = element_text(size = 25,face = "bold")) +
  ylab('Percent merged reads out of all reads') +
  guides(colour = guide_legend("Lineages")) + 
  ggtitle(paste0('Percent merged reads out of all reads over time for all lineages')) +
  theme_minimal() + 
  theme( axis.line = element_line(colour = "black"),
         axis.text=element_text(size=18), axis.title=element_text(size=20,face="bold"),
         strip.text.x = element_text(size=20,face="bold"), strip.text.y = element_text(size=20,face="bold"), 
         plot.title = element_text(size = 30,face = "bold"),
         legend.text = element_text(size=20),legend.title = element_text(size=22,face="bold")	)




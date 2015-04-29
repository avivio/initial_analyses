

require(dplyr)
#summary data

results.dir = 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\results\\'

summary_data_location = 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\final_summary_24-04.csv'
summary_data = read.csv(summary_data_location,check.names=FALSE,row.names = 1)

png(paste0(results.dir,'percent_of_library_covered.png'),units="in", width=10, height=8, res=50)
plot(as.numeric(as.character(unlist(summary_data['day',])))
     ,as.numeric(as.character(unlist(summary_data['found_designs',])))/14234*100,
     ylab = 'percent of designs in library found',xlab = 'days')
dev.off()

png(paste0(results.dir,'percent_mapped_reads_after_umi.png'),units="in", width=10, height=8, res=50)
plot(as.numeric(as.character(unlist(summary_data['day',])))
     ,as.numeric(as.character(unlist(summary_data['design_count_post_umi',])))/
       as.numeric(as.character(unlist(summary_data['design_count_pre_umi',])))*100,
     ylab = 'percent of mapped reads that pass umi filter',xlab = 'days')
dev.off()

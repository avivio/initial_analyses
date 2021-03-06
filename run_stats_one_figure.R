theme_aviv <-       theme_minimal() + 
  theme( axis.line = element_line(colour = "black"),
         axis.text=element_text(size=14), ,axis.title.x = element_text(vjust=-0.5),
         axis.title.y = element_text(vjust=1.5),strip.text.x = element_text(size=16,face="bold"), strip.text.y = element_text(size=16,face="bold"), 
         plot.title = element_text(size = 25,face = "bold"),
         legend.text = element_text(size=16),legend.title = element_text(size=18,face="bold")	)


require(dplyr)
require(tidyr)
require(ggplot2)


summary_data_location = 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\final_summary_generations.csv'
summary_data = data.frame(t(read.csv(summary_data_location,check.names=TRUE,row.names = 1)))
summary_data <- cbind(row.names(summary_data),summary_data)
names(summary_data)[1] = 'name'
names(summary_data)[12]= 'matched_reads'
summary_data <- summary_data %>% 
  select(name,all_reads,merged_reads,trimmed_reads,matched_reads) %>% 
  gather(stage,reads,-name) 
summary_data <- transform(summary_data, reads = as.numeric(reads))
  
summary_data$name = factor(summary_data$name, levels = c('Ancestor','A_28','A_56','A_84','A_112','A_140','A_168','A_196'
                                                         ,'B_28','B_56','B_84','B_112','B_140','B_168','B_196'
                                                         ,'C_28','C_56','C_84','C_112','C_140','C_168','C_196'
                                                         ,'D_28','D_56','D_84','D_112','D_140','D_168','D_196'
                                                         ,'E_28','E_56','E_84','E_112','E_140','E_168','E_196'
                                                         ,'F_28','F_56','F_84','F_112','F_140','F_168','F_196'))
png('C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\results\\for_thesis\\pipeline_stats_cool_hue_3.png',
    type="cairo",    units="in", width=39, height=15, pointsize=8, res=500) 

ggplot(summary_data,aes(x = name, y = reads, fill = stage)) + 
  geom_bar(stat = "identity",position="identity") + 
  scale_fill_hue(h=c(300,180),name="Pipeline\nStage",labels=c("All", "Merged", "Trimmed",'Matched'))  +
  ylab('Reads') + 
  xlab('Sample') + 
  theme_aviv +
  theme( axis.title=element_text(size=20 ,face="bold" ),         axis.text=element_text(size=18),
           legend.text = element_text(size=20),legend.title = element_text(size=22,face="bold"))
dev.off()


anc.1.freq = read.csv("C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\prelim_fitseq_results_130515\\fitseq_sample_data_result_1_anc.csv")


anc.1.freq[which(as.numeric(anc.1.freq[,2]) == 0),]
hist(log10(as.numeric(anc.1.freq[,2])))

fitseq.all.raw = read.csv("C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\prelim_fitseq_results_130515\\all_fitseq.csv",row.names=1)
fitseq.all.raw[fitseq.all.raw$X1_anc>50000,] <-  0


plot(fitseq.all.raw[,'X1_anc'],fitseq.all.raw[,'X2_anc']*0.22, ylim = c(0,300), xlim = c(0,300))

plot(fitseq.all[which(as.numeric(fitseq.all[,3]) <10  ) ,2],fitseq.all[which(as.numeric(fitseq.all[,3]) <10  ) ,3]* 0.22 )


as.numeric(fitseq.all[13365,2])

fitseq.all[13365,2] = 0

fitseq.all[which(as.numeric(fitseq.all[,2]) >300) ,2] = 0

results.dir = 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\results\\initial_data_analyses\\'


png(paste0(results.dir,'initial_timepoint_vs_ancestor.png'),units="in", width=15, height=15, res=100)

par(mfrow = c(7,2))

for (i in 1:7){
  
  plot(fitseq.all.raw[,'X1_anc']/sum(fitseq.all.raw[,'X1_anc']),fitseq.all.raw[,3+i]/sum(fitseq.all.raw[,3+i]), 
       xlab = 'Ancestor 1', ylab = paste('Timepoint',colnames(fitseq.all.raw)[3+i] ))
  abline(0,1)
  
  plot(fitseq.all.raw[,'X2_anc']/sum(fitseq.all.raw[,'X2_anc']),fitseq.all.raw[,3+i]/sum(fitseq.all.raw[,3+i]), 
       xlab = 'Ancestor 2', ylab = paste('Timepoint',colnames(fitseq.all.raw)[3+i] ))
  
  abline(0,1)
  
  
}
dev.off()



png(paste0(results.dir,'initial_timepoint_vs_timepoint.png'),units="in", width=15, height=15, res=100)

par(mfrow = c(7,7))

for (i in 1:7){
  for (j in 1:7){
    
    plot(fitseq.all.raw[,1+i]/sum(fitseq.all.raw[,1+i]),fitseq.all.raw[,1+j]/sum(fitseq.all.raw[,1+j]), 
         xlab = paste('Timepoint',colnames(fitseq.all.raw)[1+i] ), ylab = paste('Timepoint',colnames(fitseq.all.raw)[1+j] ))
    abline(0,1)
    
    
  }
}
dev.off()

days = c(0,4,8,12,16,20,24,28)
day_values = vector(length = 10)

for (i in c(2,4:10)){
  print(i)
  print(fitseq.all.raw[1,i])
  print(sum(fitseq.all.raw[,i]))
  day_values[i] = as.numeric(fitseq.all.raw[1,i])/as.numeric(sum(fitseq.all.raw[,i]))
  print(day_values[i])
  print(day_values)
  
}
plot(days,day_values[c(2,4:10)])





all_sums = apply(fitseq.all.raw[,3:9],2,sum)

days = c(0,4,8,12,16,20,24,28)
days_mat =matrix(days,nrow=100,ncol=length(days),byrow=TRUE)


random.normal  <-  fitseq.all.raw %>% 
  select(-X2_anc) %>% 
sample_n(100)/all_sums

normal.data  <- fitseq.all.raw %>% 
   select(-X2_anc)/all_sums

plot(days,normal.data[1,],type = "l", lty = 1)

matplot(t(days_mat),t(normal.data[1:100,]), 
         ylim = c(0,max(normal.data[1:100,]) ))

matlines(t(days_mat),t(normal.data[1:100,]), 
        ylim = c(0,max(normal.data[1:100,]) ))


random.normal  <-  fitseq.all.raw %>% 
  select(-X2_anc) %>% 
  sample_n(100)/all_sums

normal.data  <- fitseq.all.raw %>% 
  select(-X2_anc)/all_sums

fitseq_lm  <-
  normal.data %>%
  mutate(R2 = get.lm.to.days.r2(1:8),slope = get.lm.to.days.slope(1:8)) %>%
  arrange(slope,R2) %>%
  select(-slope,-R2)

glimpse(fitseq_lm)

plotnRows(paste0(results.dir,'bottom_100_slopes.png'),100,tail(fitseq_lm,100),
          "Top 100 R^2 varaints with ascending linear models")



plotnRows  <-  function(figure_name,n,rows,title){
  days_mat =matrix(days,nrow=n,ncol=length(days),byrow=TRUE)
  
  png(figure_name,units="in", width=10, height=7, res=100)
  
  print(nrow(days_mat))
  print(nrow(t(days_mat)))
  print( nrow(rows))
  print( nrow(t(rows)))
  
  matplot(t(days_mat),t(rows), 
          ylim = c(0,max(rows) ), ylab = 'DNA coverage normalized by all mapped reads in bin',
          xlab = 'days',main = title)
  
  matlines(t(days_mat),t(rows), 
           ylim = c(0,rows))
  
  dev.off()
  
}


get.lm.to.days.slope  <- function(data){
  linear = lm(unlist(data)~days)
  return(coef(linear)[2] )
}


get.lm.to.days.r2  <- function(data){
  linear = lm(unlist(data)~days)
  return(summary(linear)$adj.r.squared )
}

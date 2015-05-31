require(ggplot2)
require(dplyr)
require(stringdist)
results_dir = 'D:\\Users\\dell7\\Documents\\Tzachi\\workspace\\results\\'

fitseq.data.location  <- 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\clean_data\\goodman_salis_tuller_fitseq_2_mismatch_with_0.csv'
seq.data = read.csv(fitseq.data.location) %>% select(Name, new.name,Promoter = Promoter.Display,
                                                     RBS = RBS.Display,Gene,CDS.type,sequence = variable.seq)

seq.data  <- transform(seq.data,sequence = as.character(sequence))


distance.91  <- stringdist(seq.data.91$sequence[1],seq.data.91$sequence[-1],method = 'h')
min(distance.91)

get.min.distance  <- function(index,strings){
  
  distances  <- stringdist(seq.data.91$sequence[index],seq.data.91$sequence[-index],method = 'h')
  return(min(distances))
}

min.dists = vector(length = nrow(seq.data))
offset = 0
for (len in 91:94){
  seq.data.len <-  seq.data %>% filter(nchar(sequence)==len) %>% select(sequence)
  seq.data.len  <- unlist(seq.data.len)
  for (i in 1:length(seq.data.len)){
    min.dists[offset + i] = get.min.distance(i,seq.data.len)  
  }
  offset = (len-90) * length(seq.data.len)
  
}

png(paste0(results_dir,'min_distance_hist.png'),units= 'mm', height = 550, widt = 400)

ggplot(min.dists,aes(min.dists)) + 
  geom_histogram(breaks=c(0:10)) + 
  scale_x_continuous(breaks=c(0:10)) +
  theme_minimal() +
  ggtitle('Distribution of minimum distances between designs in libary') + 
  xlab('Minimum distance')


dev.off()

bottlenecks = read.csv('C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\bottleneck_data.csv')

positions = vector()
depths = vector()
for (i in 1:nrow(bottlenecks)){
  positions  <- append(positions, as.numeric(unlist(strsplit(names(which.max(bottlenecks[i,c(2:ncol(bottlenecks))])),split = '\\.'))[3]))
  depths  <- append(depths,max(bottlenecks[i,c(2:ncol(bottlenecks))]))
  
}
n.term.no.bn  <- read.csv('C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\n_terminal_data_no_bottleneck.csv')
TASEP.bottle.neck.position  <- positions
write.csv(cbind(n.term.no.bn,TASEP.bottle.neck.position,TASEP.bottle.neck.depth),file = 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\goodman_salis_tuller_new_name_data.csv')

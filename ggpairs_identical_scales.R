ggpairs_identical_scales  <- function(data){
  plot_matrix  <- ggpairs(data)
  
  data_min  <- floor(min(data))
  data_max  <- ceiling(max(data))
  data_limits  <- c(data_min,data_max)
  data_breaks = c(data_min:data_max)
  for (i in 1:ncol(data)){
    for (j in 1:i){
      current_plot = getPlot(plot_matrix,i,j)
      new_plot  <- current_plot +   
#         coord_fixed(ylim = data_limits,xlim = data_limits) +
        scale_y_continuous(breaks=data_breaks,limits = data_limits) +
        scale_x_continuous(breaks=data_breaks,limits = data_limits)
      if (i!= j){
        new_plot <- new_plot + geom_abline(intercept = 0,slope = 1,color='red')
      }
      plot_matrix <- putPlot(plot_matrix,new_plot,i,j)
    }
  }
  return(plot_matrix)
}





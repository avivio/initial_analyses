get.growth.rate.difference  <- function(read.ratio,generations){
  difference = log2(read.ratio)/generations
  return(difference)
}

print(get.growth.rate.difference(5*10^-5,99))
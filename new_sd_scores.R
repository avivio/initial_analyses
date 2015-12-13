t = 303.15
R = 8.3144

require(Biostrings)
require(psych)
require(dplyr)
require(ggplot)
require(tidyr)


get.sd.score <- function(seq,position){
  sd.score <- Inf
  for (i in 8:11){
    base <- position-i
    site <- substring(seq,base-7,base+2)
    sd <- tenmers[site,1]
    if (sd < sd.score){
      sd.score <- sd
      
    }
  }
  k.b = exp((sd.score/(R*t)))
  velocity = k.b/(1+k.b)
  return(velocity)
}



factor.2.rna.char <- function(fact){
  return( toupper( as.character(RNAString(DNAString(as.character(fact)))))) 
}
factor.2.rna.char <- Vectorize(factor.2.rna.char)


get.sd.propterties <- function(utr,cds){
  start.location <- nchar(utr)+1
  seq <- paste0(utr,cds)
  sd.vec <- vector(length = nchar(cds), mode = 'double')
  
  
  sd.vec <- vector(length = nchar(cds), mode = 'double')
  index <-  0
  for (i in start.location:(nchar(seq))){
    index <-  index +1
    sd <- get.sd.score(seq,i)
    sd.vec[index] <- sd
    
  }
  return(sd.vec)


}
get.sd.propterties <- Vectorize(get.sd.propterties)



print('load tenmers')

tenmers  <-   read.csv('C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\sd_tenmer_data.csv',row.names = 1)




print('get sd data for fitseq sequences')
fitseq.data <- read.csv('C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\fitseq_cds_utr.csv')

fitseq.data.test <- fitseq.data[1:5,]
sd.properties <- get.sd.propterties(factor.2.rna.char(fitseq.data$UTR),factor.2.rna.char(fitseq.data$CDS.seq))

print('transpose')
sd.properties <- data.frame(t(sd.properties))
print('numeric')

names(sd.properties) <- c('sd_pos_1','sd_pos_2','sd_pos_3','sd_pos_4','sd_pos_5','sd_pos_6','sd_pos_7','sd_pos_8','sd_pos_9','sd_pos_10',
                          'sd_pos_11','sd_pos_12','sd_pos_13','sd_pos_14','sd_pos_15','sd_pos_16','sd_pos_17','sd_pos_18','sd_pos_19','sd_pos_20',
                          'sd_pos_21','sd_pos_22','sd_pos_23','sd_pos_24','sd_pos_25','sd_pos_26','sd_pos_27','sd_pos_28','sd_pos_29','sd_pos_30',
                          'sd_pos_31','sd_pos_32','sd_pos_33')


sd.properties.summarized   <- sd.properties %>% 
  rowwise() %>% 
  do(sd_harm_mean = harmonic.mean(as.numeric(.)), sd_art_mean = mean(as.numeric(.)),
     sd_max = max(as.numeric(.)),sd_min = min(as.numeric(.)),
     sd_sum = sum(as.numeric(.)), sd_median = median(as.numeric(.)),sd_count = sum(. >  0.5 ))
     

sd.properties.summarized <- transform(sd.properties.summarized,sd_harm_mean = as.numeric(sd_harm_mean), sd_max = as.numeric(sd_max),
                                      sd_min = as.numeric(sd_min), sd_sum = as.numeric(sd_sum),sd_median = as.numeric(sd_median),
                                      sd_art_mean = as.numeric(sd_art_mean), sd_count = as.numeric(sd_count))

sd.properties.name <- bind_cols(fitseq.data %>% select(Name, new.name),sd.properties.summarized, sd.properties)

print('output data to file')
write.csv(sd.properties.name, row.names=FALSE,
          file = 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\sd_affinity_per_position_new_calc.csv')
print('Done')






# cds = 'ATGTCTCTCAACTTTCTGGATTTCGAACAACCG'
# cds <-toupper( as.character(RNAString(DNAString(cds))))
# 
#          
# utr = 'TAATAAAGAGGAGAAAtactag'
# utr <- toupper( as.character(RNAString(DNAString(utr))))
# 
# get.sd.propterties(utr,cds)

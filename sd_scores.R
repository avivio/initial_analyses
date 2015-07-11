


require(Biostrings)

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
  return(sd.score)
}



factor.2.rna.char <- function(fact){
  return( toupper( as.character(RNAString(DNAString(as.character(fact)))))) 
}
factor.2.rna.char <- Vectorize(factor.2.rna.char)


get.sd.propterties <- function(utr,cds){
  sd.properties <- list()  
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
  
  sd.properties$sd.rbs <- sd.vec[1]
  
  sd.vec.no.rbs <- sd.vec[-1]
  
  sd.properties$sd.max.score <- min(sd.vec.no.rbs,na.rm = T)
  sd.properties$sd.max.position <- which.min(sd.vec.no.rbs)
  sd.vec.no.rbs.all.negative <- sd.vec.no.rbs[sd.vec.no.rbs < 0]
  
  sd.properties$sd.mean <-  mean(sd.vec.no.rbs.all.negative,na.rm = T)
  sd.properties$sd.sdev <-  sd(sd.vec.no.rbs.all.negative,na.rm = T)
  sd.properties$sd.median <-  median(sd.vec.no.rbs.all.negative,na.rm = T)
  sd.properties$sd.count <- length(sd.vec.no.rbs.all.negative)
  
  
  full.seq <- paste0(seq,gfp)
  gfp.location <- 33
  
  sd.gfp.vec <- vector(length = 17, mode = 'double')
  index <-  0
  for (i in (gfp.location+start.location):(gfp.location+start.location+17)){
    index <-  index +1
    sd <- get.sd.score(full.seq,i)
    sd.gfp.vec[index] <- sd
    
  }
  sd.gfp.vec <- c(sd.gfp.vec,sd.late.gfp.vec)
  
  sd.properties$sd.gfp.max.score <- min(sd.gfp.vec,na.rm = T)
  sd.properties$sd.gfp.max.position <- which.min(sd.gfp.vec)
  sd.gfp.vec.no.rbs.all.negative <- sd.gfp.vec[sd.gfp.vec < 0]
  sd.properties$sd.gfp.mean <-  mean(sd.gfp.vec.no.rbs.all.negative,na.rm = T)
  sd.properties$sd.gfp.sdev <-  sd(sd.gfp.vec.no.rbs.all.negative,na.rm = T)
  sd.properties$sd.gfp.median <-  median(sd.gfp.vec.no.rbs.all.negative,na.rm = T)
  sd.properties$sd.gfp.count <- length(sd.gfp.vec.no.rbs.all.negative)
  
  return(as.data.frame(sd.properties))
  
}
get.sd.propterties <- Vectorize(get.sd.propterties)



print('load tenmers')

tenmers  <-   read.csv('C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\sd_tenmer_data.csv',row.names = 1)
print('calculate gfp sd')
gfp = 'CATATGCGTAAAGGCGAAGAGCTGTTCACTGGTTTCGTCACTATTCTGGTGGAACTGGATGGTGATGTCAACGGTCATAAGTTTTCCGTGCGTGGCGAGGGTGAAGGTGACGCAACTAATGGTAAACTGACGCTGAAGTTCATCTGTACTACTGGTAAACTGCCGGTACCTTGGCCGACTCTGGTAACGACGCTGACTTATGGTGTTCAGTGCTTTGCTCGTTATCCGGACCACATGAAGCAGCATGACTTCTTCAAGTCCGCCATGCCGGAAGGCTATGTGCAGGAACGCACGATTTCCTTTAAGGATGACGGCACGTACAAAACGCGTGCGGAAGTGAAATTTGAAGGCGATACCCTGGTAAACCGCATTGAGCTGAAAGGCATTGACTTTAAAGAAGACGGCAATATCCTGGGCCATAAGCTGGAATACAATTTTAACAGCCACAATGTTTACATCACCGCCGATAAACAAAAAAATGGCATTAAAGCGAATTTTAAAATTCGCCACAACGTGGAGGATGGCAGCGTGCAGCTGGCTGATCACTACCAGCAAAACACTCCAATCGGTGATGGTCCTGTTCTGCTGCCAGACAATCACTATCTGAGCACGCAAAGCGTTCTGTCTAAAGATCCGAACGAGAAACGCGATCACATGGTTCTGCTGGAGTTCGTAACCGCAGCGGGCATCACGCATGGTATGGATGAACTGTACAAATAA'
gfp <-toupper( as.character(RNAString(DNAString(gfp))))
sd.late.gfp.vec <- vector(length = nchar(gfp)-19, mode = 'double')
index <-  0
for (i in 19:(nchar(gfp))){
  index <-  index +1
  sd <- get.sd.score(gfp,i)
  sd.late.gfp.vec[index] <- sd
  
}


print('load fitseq data')
fitseq.data.location  <- 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\fitseq_subset.csv'
fitseq.data = read.csv(fitseq.data.location)
print('get sd data for fitseq sequences')
sd.properties <- t(get.sd.propterties(factor.2.rna.char(fitseq.data$UTR),factor.2.rna.char(fitseq.data$CDS.seq)))
sd.properties <- apply(sd.properties, 2, as.numeric)
print('binding fitseq and sd data')
fitseq.data <- cbind(fitseq.data,sd.properties)

print('output data to file')
write.csv(fitseq.data,, row.names=FALSE,
          file = 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\sd_pep_goodman_salis_tuller_fitseq_2_mismatch_with_0.csv')
print('Done')



       
# cds = 'ATGTCTCTCAACTTTCTGGATTTCGAACAACCG'
# cds <-toupper( as.character(RNAString(DNAString(cds))))
# 
#          
# utr = 'TAATAAAGAGGAGAAAtactag'
# utr <- toupper( as.character(RNAString(DNAString(utr))))
# 

# get.sd.propterties(utr,cds)

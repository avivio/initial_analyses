#add 8-11 logic to all positions
#calculate rbs score
#calculate constant region max score, max position, median score, median number of positions
#calculate variable region max score, max position, median score, median number of positions



require(Biostrings)
tenmers  <-   read.csv('C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\sd_tenmer_data.csv',row.names = 1)

gfp = 'CATATGCGTAAAGGCGAAGAGCTGTTCACTGGTTTCGTCACTATTCTGGTGGAACTGGATGGTGATGTCAACGGTCATAAGTTTTCCGTGCGTGGCGAGGGTGAAGGTGACGCAACTAATGGTAAACTGACGCTGAAGTTCATCTGTACTACTGGTAAACTGCCGGTACCTTGGCCGACTCTGGTAACGACGCTGACTTATGGTGTTCAGTGCTTTGCTCGTTATCCGGACCACATGAAGCAGCATGACTTCTTCAAGTCCGCCATGCCGGAAGGCTATGTGCAGGAACGCACGATTTCCTTTAAGGATGACGGCACGTACAAAACGCGTGCGGAAGTGAAATTTGAAGGCGATACCCTGGTAAACCGCATTGAGCTGAAAGGCATTGACTTTAAAGAAGACGGCAATATCCTGGGCCATAAGCTGGAATACAATTTTAACAGCCACAATGTTTACATCACCGCCGATAAACAAAAAAATGGCATTAAAGCGAATTTTAAAATTCGCCACAACGTGGAGGATGGCAGCGTGCAGCTGGCTGATCACTACCAGCAAAACACTCCAATCGGTGATGGTCCTGTTCTGCTGCCAGACAATCACTATCTGAGCACGCAAAGCGTTCTGTCTAAAGATCCGAACGAGAAACGCGATCACATGGTTCTGCTGGAGTTCGTAACCGCAGCGGGCATCACGCATGGTATGGATGAACTGTACAAATAA'

gfp <-toupper( as.character(RNAString(DNAString(gfp))))

cds = 'ATGAGTCTGAATTTCCTTGATTTTGAACAGCCG'
               
cds <-toupper( as.character(RNAString(DNAString(cds))))


utr = 'TTAATTCACACAGGAAAGtactag'

utr <- toupper( as.character(RNAString(DNAString(utr))))



start.location <- nchar(utr)+1

seq <- paste0(utr,cds)

substring(seq,start.location+1,nchar(seq)+1)

sd.vec <- vector(length = nchar(cds), mode = 'double')

sd.rbs <- Inf
offset  <-  -1
for (i in 8:11){
  print(i)
  base <- start.location-i
  site <- substring(seq,base-7,base+2)
  sd <- tenmers[site,1]
  print(sd)
  if (sd < sd.rbs){
    sd.rbs <- sd
    offset <-  i
  }
}

sd.vec <- vector(length = nchar(cds), mode = 'double')
sd.vec[1] = NA

for (i in (start.location + 1):(nchar(seq)+1)){
    base <- i-offset
    print(base)
    site <- substring(seq,base-7,base+2)
    sd <- tenmers[site,1]
    sd.vec[i-(start.location + 1)] <- sd
    print(sd)

}
print(sd.vec)

sd.mean <-  mean(sd.vec,na.rm = T)
sd.sdev <-  sd(sd.vec,na.rm = T)
sd.median. <-  median(sd.vec,na.rm = T)
sd.max.score <- min(sd.vec,na.rm = T)
sd.max.position <- which.min(sd.vec)

load_packages()
require(seqinr)
require(Biostrings)

goodman_data_raw <- read.csv('C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\n_terminal_data_with_full_seq.csv')

goodman_data_raw <- goodman_data_raw %>% select(Name, CDS.seq,variable.seq)


get_gc_pos_for_seq <- function(seq){
  gc_pos = vector()
  seq = unlist(strsplit(as.character(seq),''))
  gc_pos[1] = GC1(seq)
  gc_pos[2] = GC2(seq)
  gc_pos[3] = GC3(seq)
  names(gc_pos) = c('gc_1','gc_2','gc_3')
  return(gc_pos)
}
get_gc_pos_for_seq <- Vectorize(get_gc_pos_for_seq)


get_codon_enrichment <- function(seq){
  codons <- vector(mode = 'numeric',length = 64)
  names(codons) <- names(GENETIC_CODE)
  for (i in seq(1,33,3)){
    codon <- substr(seq,i,i+2)
    codons[codon] = codons[codon] + 1
  }
  names(codons) <- paste('codon',names(GENETIC_CODE),sep = '_')
  
  return(codons)
}
get_codon_enrichment <- Vectorize(get_codon_enrichment)


cds_gc_pos <- data.frame(t(get_gc_pos_for_seq(goodman_data_raw$CDS.seq)))
names(cds_gc_pos) <- c('gc_1_cds','gc_2_cds','gc_3_cds')

cds_codons <- data.frame(t(get_codon_enrichment(goodman_data_raw$CDS.seq)))
names(cds_gc_pos) <- c('gc_1_cds','gc_2_cds','gc_3_cds')



gc_pos <- bind_cols(select(goodman_data_raw,Name),cds_gc_pos,cds_codons)


write.csv(gc_pos,'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\gc_pos_codon_data.csv',row.names = F)


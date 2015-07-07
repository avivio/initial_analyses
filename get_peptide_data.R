
require(seqinr)
require(Peptides)

get.peptide.stats <- function(sequence){
  data(aacost)
pep.properties <- list()

sequence.split <- strsplit(sequence, "")[[1]]
tot.cost = 0
for (aa in sequence.split){
  cost = aacost[which(aacost[,2] == aa),'tot']
  tot.cost = cost + tot.cost
}
stats <- AAstat(sequence.split,FALSE)

pep.properties$pep.pi <- stats[[3]]
pep.properties$pep.cost <- tot.cost
pep.properties$pep.mw <- mw(sequence)
pep.properties$pep.aindex <- aindex(sequence)
pep.properties$pep.boman <- boman(sequence)
pep.properties$pep.charge <- charge(sequence)
pep.properties$pep.hmoment<- hmoment(sequence)
pep.properties$pep.hydro <- hydrophobicity(sequence)
pep.properties$pep.instability <- instaindex(sequence)

pep.properties <- data.frame(pep.properties)
aa.freq <- as.list(stats[[1]]/11)
aa.freq <- data.frame(aa.freq)[-1]
names(aa.freq) <- paste0('pep.',aaa(names(aa.freq)))


aa.class.percent <- data.frame(stats[[2]])
names(aa.class.percent) <- paste0('pep.',tolower(names(aa.class.percent)))

pep.properties <- bind_cols(pep.properties,aa.freq)
pep.properties <- bind_cols(pep.properties,aa.class.percent)

return(pep.properties)

}
get.peptide.stats <- Vectorize(get.peptide.stats)

gene.sequences = read.csv('C:\\Users\\dell7\\Documents\\Tzachi\\analysis\\R_project\\peptide amino acid sequences.csv')


peptide.properties <- t(get.peptide.stats(as.character(gene.sequences$pep.sequence)))

gene.sequences <- cbind(gene.sequences,peptide.properties)

gene.sequences <- gene.sequences %>%
  mutate(pep.2nd.aa = substr(pep.sequence,2,2))

fitseq.data.location  <- 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\clean_data\\goodman_salis_tuller_fitseq_2_mismatch_with_0.csv'
fitseq.data = read.csv(fitseq.data.location)
fitseq.data.with.pep.data <-  left_join(fitseq.data,gene.sequences,by = 'Gene')
fitseq.data.with.pep.data <-  transform(fitseq.data.with.pep.data, Prot = as.numeric(as.character(Prot))
                       ,Bin.Pct.1 = as.numeric(as.character((Bin.Pct.1))),
                       Count.RNA = as.numeric(as.character((Count.RNA))),
                       pep.sequence = as.character(pep.sequence),
                       pep.pi=as.numeric(pep.pi),
                       pep.cost=as.numeric(pep.cost),
                       pep.mw=as.numeric(pep.mw),
                       pep.aindex=as.numeric(pep.aindex),
                       pep.boman=as.numeric(pep.boman),
                       pep.charge=as.numeric(pep.charge),
                       pep.hmoment=as.numeric(pep.hmoment),
                       pep.hydro=as.numeric(pep.hydro),
                       pep.instability=as.numeric(pep.instability),
                       pep.Ala=as.numeric(pep.Ala),
                       pep.Cys=as.numeric(pep.Cys),
                       pep.Asp=as.numeric(pep.Asp),
                       pep.Glu=as.numeric(pep.Glu),
                       pep.Phe=as.numeric(pep.Phe),
                       pep.Gly=as.numeric(pep.Gly),
                       pep.His=as.numeric(pep.His),
                       pep.Ile=as.numeric(pep.Ile),
                       pep.Lys=as.numeric(pep.Lys),
                       pep.Leu=as.numeric(pep.Leu),
                       pep.Met=as.numeric(pep.Met),
                       pep.Asn=as.numeric(pep.Asn),
                       pep.Pro=as.numeric(pep.Pro),
                       pep.Gln=as.numeric(pep.Gln),
                       pep.Arg=as.numeric(pep.Arg),
                       pep.Ser=as.numeric(pep.Ser),
                       pep.Thr=as.numeric(pep.Thr),
                       pep.Val=as.numeric(pep.Val),
                       pep.Trp=as.numeric(pep.Trp),
                       pep.Tyr=as.numeric(pep.Tyr),
                       pep.tiny=as.numeric(pep.tiny),
                       pep.small=as.numeric(pep.small),
                       pep.aromatic=as.numeric(pep.aromatic),
                       pep.polar=as.numeric(pep.polar),
                       pep.charged=as.numeric(pep.charged),
                       pep.basic=as.numeric(pep.basic),
                       pep.acidic=as.numeric(pep.acidic),
                       pep.aliphatic=as.numeric(pep.aliphatic),
                       pep.non.polar=as.numeric(pep.non.polar))
write.csv(fitseq.data.with.pep.data,
          file = 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\pep_data_goodman_salis_tuller_fitseq_2_mismatch_with_0.csv')

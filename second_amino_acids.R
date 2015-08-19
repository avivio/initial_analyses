load_packages()

require(seqinr)
aa <- levels(factor(fitseq_data_residuals_above_14$pep_2nd_aa))
result_dir = 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\results\\2nd_amino_acid\\'

for (lineage_letter in c('A','B','C','D','E','F')){
for (acid in aa){
  png(paste0(result_dir,'lineage_' , lineage_letter,'_2nd_codon_enrichment_of_amino_acid_',acid, '.png'),units="in",  width=15, height=12, res=70)
  p <- ggplot(fitseq_data_residuals_above_14 %>% filter(pep_2nd_aa== acid, lineage == lineage_letter), aes(y= log2(freq_norm_anc_1),
                                             x = log2(Prot), col = pep_2nd_codon)) +
    geom_point(size = 3) +
    coord_cartesian(ylim = c(-4,5.5)) +
    ggtitle(paste('Second amino acid', aaa(acid),'lineage', lineage_letter))+
    ylab('Log 2 of frequency in sample over frequency in ancestor\n') + 
    xlab('\nLog 2 Protein level') +
    geom_line(aes(x =log2(Prot),y= prediction, inherit.aes=FALSE),
              color = 'black',size = 1.2,show_guide  = F) +
    guides(color=guide_legend(title='Codon')) +
    theme_aviv 
  print(p)
    dev.off()
}

}
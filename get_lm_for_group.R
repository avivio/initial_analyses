for (lineage.letter in c('A','B','C','D','E','F')){
  fitseq.current.sample <-  fitseq.xy %>% filter(Promoter.Display == 'High') %>%
    filter(day == day.number, lineage == lineage.letter)
  
  print(lineage.letter)
  print('above 14')
  data.lm.h = lm(log2(freq.norm.anc.1)~log2(Prot), data =filter(fitseq.current.sample, above.14 == TRUE) )
  
  slope.h <- signif(data.lm.h$coef[[2]], 5)
  p.value.h <- summary(data.lm.h)$coef[2,4]
  rmsd.h = signif(sqrt(mean(data.lm.h$residuals^2)))
  
  lm.text.h <- paste('Slope:',slope.h,'p value:',p.value.h,'RMSD:',rmsd.h )
  print(lm.text.h)
  print('below 14')
  
  data.lm.l = lm(log2(freq.norm.anc.1)~log2(Prot), data =filter(fitseq.current.sample, above.14 == FALSE) )
  
  slope.l <- signif(data.lm.l$coef[[2]], 5)
  p.value.l <- summary(data.lm.l)$coef[2,4]
  rmsd.l = signif(sqrt(mean(data.lm.l$residuals^2)))
  
  lm.text.l <- paste('Slope:',slope.l,'p value:',p.value.l,'RMSD:',rmsd.l )
  print(lm.text.l)
}

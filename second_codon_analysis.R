require(seqinr)
load_packages()
require(Biostrings)

codons <- names(GENETIC_CODE)

codon_data <- data.frame(cbind(toupper(rownames(caitab)),caitab$ec))
names(codon_data) <- c('Codon', 'CAI')
codon_data <- codon_data %>% filter(!Codon %in% c('TAG','TAA', 'TGA'))
codon_data$Codon <- factor(codon_data$Codon,
                                       levels = unique(codon_data$Codon))
tai_nte_tab <- read.csv('C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\codon_demand_supply.csv',quote = "'")
codon_data <- left_join( codon_data,select(tai_nte_tab,Codon, tAI, nTE), by = 'Codon')
codon_data <-  transform(codon_data, CAI = as.numeric(as.character(CAI)))

#getting data all
pos_neg_codons <- fitseq_data_residuals_above_14 %>% 
  select(day,Name,n_pos_neg,pep_2nd_codon) %>%
  group_by(day,n_pos_neg,pep_2nd_codon) %>%
  summarize(codon_abundance = n()) %>% 
  mutate(all_sum = sum(codon_abundance),codon_freq = codon_abundance/all_sum) %>% 
  select(-all_sum,-codon_abundance) %>% 
  spread(n_pos_neg,codon_freq, fill = 0.00001) %>% 
  mutate(ratio = Positive/Negative) %>% 
  left_join(codon_data, by = c('pep_2nd_codon' = 'Codon')) 

pos_neg_lim <- max(c(max(pos_neg_codons$Positive),max(pos_neg_codons$Negative)))


ggplot(pos_neg_codons %>% filter(day == 12),
       aes(y=Positive,x=Negative,color = nTE)) +
  geom_point(size = 5) + 
  scale_color_gradient2(low = hue_pal()(3)[3],
                        high = hue_pal()(3)[1],mid = hue_pal()(3)[2],midpoint =0.1
                        , limit = c(0.001,0.2)
                        ) +
  
  theme_aviv + 
  geom_abline(slope = 1, intercept = 0) +
  # geom_text(aes(label=pep_2nd_codon),hjust=0, vjust=0) + 
  coord_equal(xlim =c(0,pos_neg_lim + 0.01),ylim =c(0,pos_neg_lim + 0.01) )

pos_neg_codons_filter <- pos_neg_codons %>% filter(ratio < 100
                                                   # ,nTE < 0.25
                                                   ) 

metric_cor <- paste0('\nCorrelation: ',signif( cor(pos_neg_codons_filter$ratio,pos_neg_codons_filter$tAI, method = 'pearson'),2), '    \np-value: ',
                   signif(cor.test(pos_neg_codons_filter$ratio,pos_neg_codons_filter$tAI, method = 'pearson')$p.value,4), '    ')

ggplot(pos_neg_codons_filter %>% filter(day == 12),
       aes(y=ratio,x=tAI)) +
  geom_point(size = 5) + 
  theme_aviv + 
  geom_text(aes(label=paste(' ',pep_2nd_codon)),hjust=0, vjust=0) +
  geom_smooth(method = 'lm') +
  ylab('Ration of Positive codon frequency over Negative') +
  annotate ('text', x = Inf, y = Inf,label = metric_cor, hjust = 1, vjust = 1, face = 'bold', color = 'blue'  )
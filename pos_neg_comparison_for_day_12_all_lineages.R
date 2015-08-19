


#load packages and data
load_packages()

load_fitseq_data <- function(fitseq_data_location){
  
  fitseq_data = read.csv(fitseq_data_location)
  fitseq_data = transform(fitseq_data, Prot = as.numeric(as.character(Prot))
                          ,Bin_Pct_1 = as.numeric(as.character((Bin_Pct_1))),
                          Count_RNA = as.numeric(as.character((Count_RNA))))
  fitseq_data_tidy <- fitseq_data %>% select(-Count_A_DNA,-Count_A_RNA,-Count_B_DNA,-Count_B_RNA,
                                             -Bin_1,-Bin_2,-Bin_3,-Bin_4,-Bin_5,-Bin_6,-Bin_7,
                                             -Bin_8,-Bin_9,-Bin_10,-Bin_11,-Bin_12,-RNA_A,-RNA_B,
                                             -Bin_Pct_1,-Bin_Pct_2,-Bin_Pct_3,-Bin_Pct_4,-Bin_Pct_5,
                                             -Bin_Pct_6,-Bin_Pct_7,-Bin_Pct_8,-Bin_Pct_9,-Bin_Pct_10,
                                             -Bin_Pct_11,-Bin_Pct_12,-Insuff_Prot,-Insuff_DNA,
                                             -Insuff_RNA,-Fltr_BelowRange,-Fltr_AboveRange ,
                                             -Fltr_SetGood,-full_seq,-pep_sequence,-UTR,-CDS_seq,-Promoter_seq,-RBS_seq,-Promoter,
                                             -variable_seq,-full_peptide,-preRBS,-salis_status)
  fitseq_data_tidy  <- fitseq_data_tidy %>%   gather(sample, frequency, D_12:F_12)
  fitseq_data_tidy  <- fitseq_data_tidy %>% separate(sample, c("lineage", "day"), '_',convert = T)
  fitseq_data_tidy$RBS_Display <- factor(fitseq_data_tidy$RBS_Display,
                                         levels = c("Strong", "Mid", "Weak", "WT"))
  
  
  
  
  fitseq_data_tidy <- fitseq_data_tidy %>% 
    group_by(day,lineage) %>% 
    mutate(sum_sample= sum(frequency),norm_freq = frequency/sum_sample,
           sum_anc= sum(anc_1),norm_anc = anc_1/sum_anc,freq_norm_anc = norm_freq/norm_anc,
           sum_sample_1= sum(frequency+1),norm_freq_1 = (frequency+1)/sum_sample_1,
           sum_anc_1= sum(anc_1+1),norm_anc_1 = (anc_1+1)/sum_anc_1,freq_norm_anc_1 = norm_freq_1/norm_anc_1)
  
  
  return(fitseq_data_tidy)
}



compare_pos_neg <- function(fitseq_data_prediction,variable,label,day_number,lineage_letter,data_set_name
                            ,result_dir,y_lim_all,y_lim_rbs){
  
  fitseq_current_sample <-  fitseq_data_prediction %>% 
    filter(day == day_number, lineage == lineage_letter)
  fitseq_current_sample_no_wt <-  fitseq_current_sample %>% 
    filter(RBS_Display != 'WT')
  
  
  neg  <- fitseq_current_sample %>%
    ungroup %>% 
    filter(pos_neg== 'Negative') %>%
    mutate_(new_variable = variable)  %>% select(new_variable) 
  neg <- as.numeric( unlist(neg))  
  pos  <- fitseq_current_sample %>%
    ungroup %>% 
    filter(pos_neg== 'Positive') %>%
    mutate_(new_variable = variable)  %>% select(new_variable) 
  pos <- as.numeric( unlist(pos)) 
  
  wilcox_string_all <- get_wilcox_string(pos,neg)
  
  
  fitseq_current_sample <- fitseq_current_sample %>%
    mutate(wilcox_all = wilcox_string_all)
  
  rbs_list = c('Strong','Mid','Weak')
  wilcox_string_rbs = vector(mode = 'character',length = 3)
  names(wilcox_string_rbs) <- rbs_list
  
  for (rbs in rbs_list){
    neg_rbs  <- fitseq_current_sample_no_wt %>%
      ungroup %>% 
      filter(pos_neg== 'Negative',RBS_Display == rbs) %>%
      mutate_(new_variable = variable)  %>% select(new_variable) 
    neg_rbs <- as.numeric( unlist(neg_rbs))  
    pos_rbs  <- fitseq_current_sample_no_wt %>%
      ungroup %>% 
      filter(pos_neg== 'Positive',RBS_Display == rbs) %>%
      mutate_(new_variable = variable)  %>% select(new_variable) 
    pos_rbs <- as.numeric( unlist(pos_rbs)) 
    
    wilcox_string_rbs[rbs] <- get_wilcox_string(pos_rbs,neg_rbs)
    
    
    
  }
  fitseq_current_sample_no_wt <- fitseq_current_sample_no_wt %>%
    mutate(wilcox_rbs = wilcox_string_rbs[RBS_Display])
  
  title <- paste('Comparing',label, 'between positive and negative fitness residuals',
                 '\nLineage', lineage_letter,'day', day_number, gsub('_',' ',data_set_name))
  
  png(paste0(result_dir,'all_pos_neg_hist_' ,variable,'_lineage_',lineage_letter,'_',data_set_name, '.png'),
      units="in",  width=15, height=12, res=70)
  p <- ggplot(fitseq_current_sample, aes_string(x = variable,fill = 'pos_neg')) +
    geom_density(alpha = 0.4) +
    xlab(paste0(label,'\n'))+
    expand_limits(y = y_lim_all  ) +
    ylab('\nDensity') +  
    ggtitle(title) +
    geom_text(aes( y = Inf, x = -Inf, label = wilcox_all), color="black",hjust = 0,vjust = 1, size = 6,
              parse=FALSE,fontface = 'bold',colour="red") +
    guides(fill=guide_legend(title='Fitness\nresidual')) +
    theme_aviv
  print(p)
  dev.off()
  
  png(paste0(result_dir,'promoter_rbs_pos_neg_hist_' ,variable,'_lineage_',lineage_letter,'_',data_set_name, '.png'),
      units="in",  width=15, height=12, res=70)
  p <- ggplot(fitseq_current_sample_no_wt, aes_string(x = variable,fill = 'pos_neg')) +
    geom_density(alpha = 0.4) +
    facet_grid(RBS_Display ~ .) +
    expand_limits(y = y_lim_rbs  ) +
    xlab(paste0(label,'\n')) +
    ylab('\nDensity') +  
    ggtitle(title) +
    geom_text(aes( y = Inf, x = -Inf, label = wilcox_rbs), color="black",hjust = 0,vjust = 1, size = 6,
              parse=FALSE,fontface = 'bold',colour="red") +
    guides(fill=guide_legend(title='Fitness\nresidual')) +
    theme_aviv
  print(p)
  dev.off()
  
  
  
}




get_data_for_correlation_category <- function(source_data,y_string,y_label,x_string,x_label,col_string,col_label){
  
  
  
  
  fitseq_xy <-  source_data %>%
    mutate_(y = y_string,x = x_string,col = col_string)
  
  fitseq_xy <- fitseq_xy %>%
    filter(!is.na(y) & !is.na(x) & !is.na(col))
  
  fitseq_xy_lm_all <- fitseq_xy %>% 
    ungroup() %>%
    group_by(day,lineage,col) %>% 
    do(mod = lm(y ~ x, data = .)) %>%
    mutate(slope = summary(mod)$coeff[2],
           intercept = summary(mod)$coeff[1],
           rmsd = sqrt(mean(mod$residuals^2)),
           p_value = summary(mod)$coeff[2,4],
           lm_string =paste('\n',  'slope = ', signif(slope,2),'\n',
                            'p-value = ' ,formatC(p_value,format='e',digits = 2),'\n',
                            'RMSD = ', signif(rmsd,2),'\n')) %>% 
    select(-mod)
  
  fitseq_xy <- left_join(fitseq_xy,fitseq_xy_lm_all,by = c('day','lineage','col'))
  
  fitseq_xy_lm_group <- fitseq_xy %>% 
    ungroup() %>%
    group_by(day,lineage,col,Promoter_Display,RBS_Display ) %>%
    do(mod = lm(y ~ x, data = .)) %>%
    mutate(slope = summary(mod)$coeff[2],
           intercept = summary(mod)$coeff[1],
           rmsd = sqrt(mean(mod$residuals^2)),
           p_value = summary(mod)$coeff[2,4],
           lm_string_group = paste('\n',  'slope = ', signif(slope,2),'\n',
                                   'p-value = ' ,formatC(p_value,format='e',digits = 2),'\n',
                                   'RMSD = ', signif(rmsd,2),'\n')) %>% 
    select(day,lineage,col,Promoter_Display,RBS_Display ,lm_string_group)
  
  
  fitseq_xy <- left_join(fitseq_xy,fitseq_xy_lm_group,by = c('day','lineage','col','Promoter_Display','RBS_Display'))
  
  
  fitseq_xy <-  fitseq_xy %>% 
    mutate(text_loc_y = ifelse(above_14,-Inf,Inf),
           text_loc_x= ifelse(above_14,Inf,-Inf),
           hor = ifelse(above_14,1,0),
           ver = ifelse(above_14,0,1))
  
  fitseq_xy <-  fitseq_xy %>% 
    mutate(prediction = intercept + slope*log2(Prot),
           fit_resid = log2(freq_norm_anc_1) -  prediction, 
           pos_neg = ifelse(fit_resid > 0, 'Positive', 'Negative') )
  return(fitseq_xy)
  
}

col_plots_factor <- function(fitseq_xy,day_number,lineage_letter,data_set_name,result_dir,
                             y_string,y_label,x_string,x_label,col_string,col_label){
  
  y_limits_all <- c(-5.5,4)
  y_limits_rbs_promoter <- c(-4.7,2)
  
  fitseq_current_sample <-  fitseq_xy %>% filter(day == day_number, lineage == lineage_letter)
  fitseq_current_sample_no_wt <-  fitseq_xy %>% filter(day == day_number, lineage == lineage_letter, RBS_Display != 'WT')
  
  
  filename_all  <- paste('lineage', lineage_letter,'all','day', day_number,'x',x_string ,'vs',y_string, 'col', col_string,data_set_name,sep = '_')
  filename_all  <-clean_filename(filename_all)
  title <- paste(x_label ,'vs',y_label,'\nfor day', day_number, 'lineage', lineage_letter
                 ,gsub('_',' ',data_set_name))
  
  png(paste0(result_dir,filename_all,'.png'),units="in",  width=15, height=12, res=70)
  
  
  p <- ggplot(fitseq_current_sample,aes(x= x, y = y,color  = col)) +
    geom_point(alpha = 0.3)   +
    # xlim(x_limits) +
    coord_cartesian(ylim = y_limits_all) +
    ylab(paste0(y_label,'\n')) + 
    xlab(paste0('\n',x_label)) +
    ggtitle(title) + 
    geom_line(aes(x = x,y= prediction, color = col,inherit.aes=FALSE),
              size = 1.2,show_guide  = F)+
    # geom_abline(aes(slope = slope, intercept = intercept), color = 'dodgerblue', size = 1.2) +
    #     geom_smooth(aes_string(linetype = col_string),show_guide  = F,
    #                 color = 'red', size = 1.2, method = 'lm') +
    geom_text(aes(x=text_loc_x, y=text_loc_y, face="bold", inherit.aes=FALSE,
                  parse=FALSE,label=lm_string,colour=col,hjust=hor, vjust=ver),show_guide  = F) +
    guides(colour = guide_legend(override.aes = list(alpha = 1),title = col_label)) +
    theme_aviv
  print(p)
  
  dev.off()
  
  filename_promoter_rbs  <- paste('lineage', lineage_letter,'promoter_rbs','day', day_number,'x',x_string ,'vs',y_string, 'col',data_set_name, col_string,sep = '_')
  
  filename_promoter_rbs  <-clean_filename(filename_promoter_rbs)
  
  png(paste0(result_dir,filename_promoter_rbs,'.png'),units="in",  width=15, height=12, res=70)
  p <- ggplot(fitseq_current_sample_no_wt,aes(x= x, y = y,color  = col)) +
    geom_point(alpha = 0.3)   +
    # xlim(x_limits) +
    coord_cartesian(ylim = y_limits_rbs_promoter) +
    ylab(paste0(y_label,'\n')) + 
    xlab(paste0('\n',x_label)) +
    ggtitle(title) + 
    geom_line(aes(x = x,y= prediction, color = col, inherit.aes=FALSE),
              size = 1.2,show_guide  = F)+
    #     geom_smooth(
    #       aes_string(linetype = col_string),
    #       color = 'red', size = 1_2, method = 'lm',show_guide  = F) +
    facet_grid(Promoter_Display ~ RBS_Display)  +
    geom_text(aes(x=text_loc_x, y=text_loc_y, face="bold", inherit.aes=FALSE,
                  parse=FALSE,label=lm_string_group,colour=col,hjust=hor, vjust=ver),show_guide  = F) +
    guides(colour = guide_legend(override.aes = list(alpha = 1),title = col_label)) +
    theme_aviv
  
  print(p)
  dev.off()
}

load_packages()
fiseq_data_location <- 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\fitseq_data_with_meta_day_12.csv'
fitseq_data_residuals <- load_fitseq_data(fiseq_data_location)
fitseq_data_residuals <-fitseq_data_residuals %>% 
  mutate(above_14 = log2(Prot)> 14) %>% 
  filter(log2(Prot)<17.5,    Promoter_Display=='High')

y_string <-  'log2(freq_norm_anc_1)'
y_label <- 'Log 2 of frequency in sample over frequency in ancestor'
x_string  <- 'log2(Prot)'
x_label <- 'Log 2 Protein level'
col_string  <- 'above_14'
# col_label <- 'Number of ribosomes\nper mRNA \n(calulated by TASEP)'
col_label <- 'Protein level\nabove 14'



# get correlation data into data set
fitseq_data_residuals <- get_data_for_correlation_category(fitseq_data_residuals,y_string,y_label,x_string,x_label,col_string,col_label)

fit_summary  <- fitseq_data_residuals %>% 
  group_by(day,lineage,pos_neg) %>% 
  summarise(percentile= quantile(abs(fit_resid),probs = c(0.2))) 

fitseq_data_residuals <-   left_join(fitseq_data_residuals,fit_summary,by = c('day','lineage','pos_neg'))



fitseq_data_residuals <- fitseq_data_residuals %>% 
  filter( 
    abs(fit_resid) >  percentile  )



fitseq_data_residuals_above_14 <- fitseq_data_residuals %>% filter(above_14 == TRUE)




data_set_name <- 'high_promoter_above_14_below_wall_80_percentile'



lineage_letter = 'C'
day_number = 12
base_result_dir = 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\results\\fitness_residuals\\gc_pos\\'
# base_result_dir= paste0(base_result_dir,data_set_name,'\\')
dir.create(base_result_dir)
result_dir <- paste0(base_result_dir,'lineages_histograms\\')
dir.create(result_dir)


variables <-     c('gc_1_cds', 'gc_2_cds','gc_3_cds')

labels <-     c('Coding sequnce GC 1 content', 'Coding sequnce GC 2 content ', 'Coding sequnce GC 3 content')

# variables <- 
#   c('log2(Trans)','log2(RNA)')
# labels <- 
#   c('Log 2 Translation efficiency','Log 2 RNA level')

# variables <-     c('sd_rbs', 'sd_max_score', 'sd_max_position', 'sd_mean', 'sd_sdev', 'sd_median', 'sd_count', 'sd_gfp_max_score',
#                    'sd_gfp_max_position', 'sd_gfp_mean', 'sd_gfp_sdev', 'sd_gfp_median', 'sd_gfp_count')
# 
# labels <-     c('RBS SD affinity score', 'Variable region max SD affinity score ', 'Variable region max SD affinity score position',
#                 'Variable region mean SD affinity score', 'Variable region mean SD affinity score stdev', 'Variable region median SD affinity score',
#                 'Variable region SD affinity score count','GFP region max SD affinity score ', 'GFP region max SD affinity score position',
#                 'GFP region mean SD affinity score', 'GFP region mean SD affinity score stdev', 'GFP region median SD affinity score',
#                 'GFP region SD affinity score count')
# 
# 




#   variables <- 
#     c('pep_pi','pep_cost','pep_mw','pep_aindex','pep_boman','pep_charge','pep_hmoment','pep_hydro','pep_instability','pep_Ala','pep_Cys','pep_Asp','pep_Glu',
#       'pep_Phe','pep_Gly','pep_His','pep_Ile','pep_Lys','pep_Leu','pep_Met','pep_Asn','pep_Pro','pep_Gln','pep_Arg','pep_Ser','pep_Thr','pep_Val','pep_Trp','pep_Tyr',
#       'pep_tiny','pep_small','pep_aliphatic','pep_aromatic','pep_non_polar','pep_polar','pep_charged','pep_basic','pep_acidic')
#   labels <- 
#     c('Peptide pI','Peptide cost','Peptide molecular weight','Peptide aliphatic index','Peptide Boman index','Peptide charge','Peptide hydrophyicl moment',
#       'Peptide hydrophobicity','Peptide instability index','Peptide Alanine content','Peptide Cystine content','Peptide Aspartate content',
#       'Peptide glutamine content','Peptide Phenylalanine content','Peptide Glycine content','Peptide histidine content','Peptide Isoleucine content',
#       'Peptide Lysine content','Peptide Leucine content','Peptide Methionine content','Peptide Asparginine content','Peptide Proline content',
#       'Peptide Glutamine content','Peptide Arginine content','Peptide Serine content','Peptide Threonine content','Peptide Valine content',
#       'Peptide Tryptophan content','Peptide Tyrosine content','Peptide tiny amino acid content','Peptide small amino acid content',
#       'Peptide aliphatic amino acid content','Peptide aromitc amino acid content','Peptide non polar amino acid content','Peptide polar amino acid content',
#       'Peptide charged amino acid content','Peptide basic amino acid content','Peptide acidic amino acid content')
#   
#   names(labels) <- variables


#   variables <- 
#     c('Rel_Codon_Freq','CAI','tAI','CDS_GC','GC','dG','dG_noutr','dG_unif','log2(salis_init)',
#       'TASEP_avgRiboNum',	'TASEP_density_0',	'TASEP_bottle_neck_position',	'TASEP_bottle_neck_depth')
#   labels <- 
#     c('Relative codon frequency','CAI','tAI','Coding sequence GC %',' GC %','Delta G','Delta G no UTR',
#       'Delta G cut at -5','Log 2 RBS calculator initition rate',
#       'Average ribosome number (calculated using TASEP)','Density at start codon (calculated using TASEP)',
#       'Bottle neck position (calculated using TASEP)','Bottle neck depth (calculated using TASEP)')
names(labels) <- variables


for (lineage_letter in c('A','B','C','D','E','F')){
  print(lineage_letter)
  #   for (day_number in seq(4,28,4)){
  #     print(day_number)
  
  for (i in 1:length(labels)){
    variable  <- names(labels)[i]
    label <- labels[i]
    print(label)
    fitseq_data_residuals_above_14 <- fitseq_data_residuals_above_14 %>%
      mutate_(var = variable)
    dens_all <- fitseq_data_residuals_above_14 %>%
      filter(day == day_number) %>%
      group_by(pos_neg,lineage) %>%
      summarise(group_dens_all = max(unlist(density(var)[2]))) %>%
      ungroup() %>%
      summarise(max_dens_all = max(group_dens_all))
    y_lim_all <- as.numeric( unlist(dens_all))
    dens_rbs = vector(mode = 'numeric',length = 3)
    
    dens_rbs <- fitseq_data_residuals_above_14 %>%
      filter(day == day_number,RBS_Display != 'WT') %>%
      group_by(pos_neg,lineage,RBS_Display) %>%
      summarise(group_dens_all = max(unlist(density(var)[2]))) %>%
      ungroup() %>%
      summarise(max_dens_all = max(group_dens_all))
    y_lim_rbs <- max(dens_rbs)
    
    
    compare_pos_neg(fitseq_data_residuals_above_14,variable,label,day_number,lineage_letter,data_set_name,
                    result_dir,y_lim_all,y_lim_rbs)
  }
}
# }


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



get_wilcox_string_all_lineages  <-  function(legend_label,legend_string_true,legend_string_false){
  
  greater <- wilcox_test(legend_string_true,legend_string_false,alternative = 'greater')$p.value
  less <- wilcox.test(legend_string_true,legend_string_false,alternative = 'less')$p.value
  if(greater < less){
    return(paste('\n',legend_label ,'\n True > False\n p value =',formatC(greater,format='e',digits = 2)))
  } else{
    return(paste('\n',legend_label ,'\n True < False\n p value =',formatC(less,format='e',digits = 2)))
  }
  
  
}

compare_pos_neg_all_lineages  <- function(fitseq_data_prediction,variable,label,day_number,data_set_name
                                          ,result_dir,y_lim_all,y_lim_rbsm,legend_string,legend_label){
  
  fitseq_current_sample <-  fitseq_data_prediction %>% 
    mutate_(legend_string = legend_string)%>%
    filter(day == day_number)
  fitseq_current_sample_no_wt <-  fitseq_current_sample %>% 
    filter(RBS_Display != 'WT')
  
  
  legend_string_negative  <- fitseq_current_sample %>%
    ungroup %>% 
    filter(legend_string== 'Negative') %>% mutate_(new_variable = variable)  %>% select(new_variable) 
  legend_string_negative <- as.numeric( unlist(legend_string_negative))  
  legend_string_positive  <- fitseq_current_sample %>%
    ungroup %>% 
    filter(legend_string== 'Positive') %>%
    mutate_(new_variable = variable)  %>% select(new_variable) 
  legend_string_positive <- as.numeric( unlist(legend_string_positive)) 
  
  wilcox_string_all <- get_wilcox_string(legend_string_positive,legend_string_negative)
  
  
  fitseq_current_sample <- fitseq_current_sample %>%
    mutate(wilcox_all = wilcox_string_all)
  
  rbs_list = c('Strong','Mid','Weak')
  wilcox_string_rbs = vector(mode = 'character',length = 3)
  names(wilcox_string_rbs) <- rbs_list
  
  for (rbs in rbs_list){
    legend_string_negative_rbs  <- fitseq_current_sample_no_wt %>%
      ungroup %>% 
      filter(legend_string== 'Negative',RBS_Display == rbs) %>%
      mutate_(new_variable = variable)  %>% select(new_variable)  
    legend_string_negative_rbs <- as.numeric( unlist(legend_string_negative_rbs))  
    legend_string_positive_rbs  <- fitseq_current_sample_no_wt %>%
      ungroup %>% 
      filter(legend_string== 'Positive',RBS_Display == rbs) %>%
      mutate_(new_variable = variable)  %>% select(new_variable) 
    legend_string_positive_rbs <- as.numeric( unlist(legend_string_positive_rbs)) 
    
    wilcox_string_rbs[rbs] <- get_wilcox_string(legend_string_positive_rbs,legend_string_negative_rbs)
    
    
    
  }
  fitseq_current_sample_no_wt <- fitseq_current_sample_no_wt %>%
    mutate(wilcox_rbs = wilcox_string_rbs[RBS_Display])
  
  title <- paste('Comparing',label, 'between positive and negative fitness residuals',
                 '\nDay', day_number, gsub('_',' ',data_set_name))
  
  png(paste0(result_dir,'all_pos_neg_hist_' ,variable,'_',data_set_name, '.png'),
      units="in",  width=15, height=12, res=70)
  p <- ggplot(fitseq_current_sample, aes_string(x = variable,fill =legend_string )) +
    geom_density(alpha = 0.4) +
    xlab(paste0(label,'\n'))+
    expand_limits(y = y_lim_all  ) +
    ylab('\nDensity') +  
    ggtitle(title) +
    geom_text(aes( y = Inf, x = -Inf, label = wilcox_all), color="black",hjust = 0,vjust = 1, size = 6,
              parse=FALSE,fontface = 'bold',colour="red") +
    guides(fill=guide_legend(title=legend_label)) +
    theme_aviv
  print(p)
  dev.off()
  
  png(paste0(result_dir,'promoter_rbs_pos_neg_hist_' ,variable,'_',data_set_name, '.png'),
      units="in",  width=15, height=12, res=70)
  p <- ggplot(fitseq_current_sample_no_wt, aes_string(x = variable,fill =legend_string)) +
    geom_density(alpha = 0.4) +
    facet_grid(RBS_Display ~ .) +
    expand_limits(y = y_lim_rbs  ) +
    xlab(paste0(label,'\n')) +
    ylab('\nDensity') +  
    ggtitle(title) +
    geom_text(aes( y = Inf, x = -Inf, label = wilcox_rbs), color="black",hjust = 0,vjust = 1, size = 6,
              parse=FALSE,fontface = 'bold',colour="red") +
    guides(fill=guide_legend(title=legend_label)) +
    theme_aviv
  print(p)
  dev.off()
  
  
  
}



filter_n_fitness_residuals_csv <- function(n,fitseq_data_residuals){
  
  fitseq_data_pos_neg_lineages <- 
    fitseq_data_residuals %>%
    select(Name, day, lineage, pos_neg)

  
  fitseq_data_pos_neg_lineages <-  fitseq_data_pos_neg_lineages %>% 
    spread(lineage,pos_neg)
  
  n_pos <- fitseq_data_pos_neg_lineages %>%
    mutate(A= as.numeric(A=='Positive'),B= as.numeric(B=='Positive'),C= as.numeric(C=='Positive'),
           D= as.numeric(D=='Positive'),E=as.numeric( E=='Positive'),F= as.numeric(F=='Positive')) %>%
    mutate(sum_lineages = A+B+C+D+E+F,n_pos_neg = ifelse(sum_lineages >= n,'Positive',NA) ) %>%
    select(Name,day,n_pos_neg) %>%
    filter(!is.na(n_pos_neg))
  
  n_neg <- fitseq_data_pos_neg_lineages %>%
    mutate(A= as.numeric(A=='Negative'),B= as.numeric(B=='Negative'),C= as.numeric(C=='Negative'),
           D= as.numeric(D=='Negative'),E=as.numeric( E=='Negative'),F= as.numeric(F=='Negative')) %>%
    mutate(sum_lineages = A+B+C+D+E+F,n_pos_neg = ifelse(sum_lineages >= n,'Negative',NA) ) %>%
    select(Name,day,n_pos_neg) %>%
    filter(!is.na(n_pos_neg))
  n_pos_neg <-  bind_rows(n_pos,n_neg)
  
  
  
  return(collect(inner_join(fitseq_data_residuals,n_pos_neg, copy = TRUE)))
  
  
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

fitseq_data_residuals <- filter_n_fitness_residuals_csv(5,fitseq_data_residuals)

fitseq_data_residuals <- fitseq_data_residuals %>%  
  filter(lineage =='A') %>% ungroup() %>% select(-lineage,-pos_neg,-fit_resid,-prediction,-p_value,-rmsd,-slope,-freq_norm_anc_1,-norm_anc_1,-sum_anc_1,-norm_freq_1,
                                                 -sum_sample_1,-freq_norm_anc,-norm_anc,-sum_anc,-norm_freq,-sum_sample,-frequency,-prediction)

fitseq_data_residuals_above_14 <- fitseq_data_residuals %>% filter(above_14 == TRUE)


base_result_dir = 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\results\\fitness_residuals\\gc_pos\\'
data_set_name <- 'high_promoter_above_14_below_wall_5_or_above'
day_number = 12
legend_string <- 'n_pos_neg'
legend_label <- '5 lineages\nsame sign\nfitness residual'

result_dir <- paste0(base_result_dir,'above_n_histograms\\')
dir.create(result_dir)

variables <- 
  c('CAI','tAI','CDS.GC','dG','dG.noutr',	'TASEP.bottle.neck.position',	'TASEP.bottle.neck.depth','pep.cost','pep.hydro','pep.polar',
    'pep.mw','pep.basic','pep.acidic', 'sd.max.score','sd.mean',  'sd.count')
labels <- 
  c('CAI','tAI','Coding sequence GC %','Delta G','Delta G no UTR','Bottle neck position (calculated using TASEP)',
    'Bottle neck depth (calculated using TASEP)','Peptide cost','Peptide hydrophobicity','Peptide polar amino acid content',
    'Peptide molecular weight','Peptide basic amino acid content','Peptide acidic amino acid content',
    'Variable region max SD affinity score ','Variable region mean SD affinity score','Variable region SD affinity score count')


names(labels) <- variables
for (i in 1:length(labels)){
  variable  <- names(labels)[i]
  label <- labels[i]
  print(label)
  fitseq_data_residuals_above_14 <- fitseq_data_residuals_above_14 %>%
    mutate_(var = variable)
  dens_all <- fitseq_data_residuals_above_14 %>%
    filter(day == day_number) %>%
    group_by(n_pos_neg) %>%
    summarise(group_dens_all = max(unlist(density(var)[2]))) %>%
    ungroup() %>%
    summarise(max_dens_all = max(group_dens_all))
  y_lim_all <- as.numeric( unlist(dens_all))
  dens_rbs = vector(mode = 'numeric',length = 3)
  
  dens_rbs <- fitseq_data_residuals_above_14 %>%
    filter(day == day_number,RBS_Display != 'WT') %>%
    group_by(n_pos_neg,RBS_Display) %>%
    summarise(group_dens_all = max(unlist(density(var)[2]))) %>%
    ungroup() %>%
    summarise(max_dens_all = max(group_dens_all))
  y_lim_rbs <- max(dens_rbs)
  
  
  compare_pos_neg_all_lineages(fitseq_data_residuals_above_14,variable,label,day_number,data_set_name,
                               result_dir,y_lim_all,y_lim_rbs,legend_string,legend_label)
}

# }




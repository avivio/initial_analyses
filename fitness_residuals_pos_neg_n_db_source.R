load_packages()

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

col_plots_factor_all_lineages <- function(fitseq_xy,day_number,data_set_name,result_dir,
                                          y_string,y_label,x_string,x_label,col_string,col_label){
  
  y_limits_all <- c(-5.5,4)
  y_limits_rbs_promoter <- c(-4.7,2)
  
  fitseq_current_sample <-  fitseq_xy %>% filter(day == day_number)
  fitseq_current_sample_no_wt <-  fitseq_xy %>% filter(day == day_number,RBS_Display != 'WT')
  
  
  filename_all  <- paste('all','day', day_number,'x',x_string ,'vs',y_string, 'col', col_string,data_set_name,sep = '_')
  filename_all  <-clean_filename(filename_all)
  title <- paste(x_label ,'vs',y_label,'\nfor day', day_number,gsub('_',' ',data_set_name))
  
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
  
  filename_promoter_rbs  <- paste('promoter_rbs','day', day_number,'x',x_string ,'vs',y_string, 'col',data_set_name, col_string,sep = '_')
  
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
    #       color = 'red', size = 1.2, method = 'lm',show_guide  = F) +
    facet_grid(Promoter_Display ~ RBS_Display)  +
    geom_text(aes(x=text_loc_x, y=text_loc_y, face="bold", inherit.aes=FALSE,
                  parse=FALSE,label=lm_string_group,colour=col,hjust=hor, vjust=ver),show_guide  = F) +
    guides(colour = guide_legend(override.aes = list(alpha = 1),title = col_label)) +
    theme_aviv
  
  print(p)
  dev.off()
}


filter_n_fitness_residuals_db <- function(n,fitseq_data_residuals){
  
  fitseq_data_pos_neg_lineages <- 
    fitseq_data_residuals %>%
    select(Name, day, lineage, pos_neg)
  fitseq_data_pos_neg_lineages <- collect(fitseq_data_pos_neg_lineages)
  
  
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

fitseq_data_residuals <- load_db()

fitseq_data_residuals <- fitseq_data_residuals %>% 
  filter( 
    day == 12,
    Promoter_Display=='High'
    ,above_14 == TRUE
    ,below_wall == TRUE
    # ,abs(fit_resid) >  percentile
  )

fitseq_data_residuals <- filter_n_fitness_residuals_db(5,fitseq_data_residuals)
fitseq_data_residuals$RBS_Display <- factor(fitseq_data_residuals$RBS_Display,
                                            levels = c('Strong','Mid','Weak','WT'))

fitseq_data_residuals <- fitseq_data_residuals %>%  
  filter(lineage =='A') %>% ungroup() %>% select(-lineage,-pos_neg,-fit_resid,-prediction,-p_value,-rmsd,-slope,-freq_norm_anc_1,-norm_anc_1,-sum_anc_1,-norm_freq_1,
                                                 -sum_sample_1,-freq_norm_anc,-norm_anc,-sum_anc,-norm_freq,-sum_sample,-frequency,-prediction)
fitseq_data_residuals_above_14 <- fitseq_data_residuals %>%  filter(log2(Prot) > 14)

base_result_dir = 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\results\\fitness_residuals\\sd_scores\\'
data_set_name <- 'high_promoter_above_14_below_wall_5_or_above'
day_number = 12
legend_string <- 'n_pos_neg'
legend_label <- '5 lineages\same sign\fitness residual'

result_dir <- paste0(base_result_dir,'above_n_histograms\\')
dir.create(result_dir)

variables <-     c('sd_rbs', 'sd_max_score', 'sd_max_position', 'sd_mean', 'sd_sdev', 'sd_median', 'sd_count', 'sd_gfp_max_score',
                   'sd_gfp_max_position', 'sd_gfp_mean', 'sd_gfp_sdev', 'sd_gfp_median', 'sd_gfp_count')

labels <-     c('RBS SD affinity score', 'Variable region max SD affinity score ', 'Variable region max SD affinity score position',
                'Variable region mean SD affinity score', 'Variable region mean SD affinity score stdev', 'Variable region median SD affinity score',
                'Variable region SD affinity score count','GFP region max SD affinity score ', 'GFP region max SD affinity score position',
                'GFP region mean SD affinity score', 'GFP region mean SD affinity score stdev', 'GFP region median SD affinity score',
                'GFP region SD affinity score count')


# variables <- 
#   c('log2(Trans)','log2(RNA)')
# labels <- 
#   c('Log 2 Translation efficiency','Log 2 RNA level')


# variables <- 
#   c('Rel_Codon_Freq','CAI','tAI','CDS_GC','GC','dG','dG_noutr','dG_unif','log2(salis_init)',
#     'TASEP_avgRiboNum',	'TASEP_density_0',	'TASEP_bottle_neck_position',	'TASEP_bottle_neck_depth')
# labels <- 
#   c('Relative codon frequency','CAI','tAI','Coding sequence GC %',' GC %','Delta G','Delta G no UTR',
#     'Delta G cut at -5','Log 2 RBS calculator initition rate',
#     'Average ribosome number (calculated using TASEP)','Density at start codon (calculated using TASEP)',
#     'Bottle neck position (calculated using TASEP)','Bottle neck depth (calculated using TASEP)')


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




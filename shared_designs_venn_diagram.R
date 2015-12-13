require(dplyr)
require(tidyr)
require(ggplot2)
require(colorspace)
require(scales)
require(VennDiagram)

load_fitseq_data <- function(fitseq_data_location){
  
  fitseq_data_tidy = read.csv(fitseq_data_location)
  fitseq_data_tidy = transform(fitseq_data_tidy, Prot = as.numeric(as.character(Prot)))
  
  fitseq_data_tidy  <- fitseq_data_tidy %>%   gather(sample, frequency, A_24:F_0)
  fitseq_data_tidy  <- fitseq_data_tidy %>% separate(sample, c("lineage", "day"), '_',convert = T)
#   fitseq_data_tidy$RBS_Display <- factor(fitseq_data_tidy$RBS_Display,
#                                          levels = c("Strong", "Mid", "Weak", "WT"))
#   
#   
  
  
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

fiseq_data_location <- 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\fitseq_and_prot_data_no_meta.csv'
fitseq_data_residuals <- load_fitseq_data(fiseq_data_location)
fitseq_data_residuals <-fitseq_data_residuals %>% 
  mutate(above_14 = log2(Prot)> 14) %>% 
  filter(log2(Prot)<17.5
         ,    Promoter_Display=='High'
  )

y_string <-  'log2(freq_norm_anc_1)'
y_label <- 'Log 2 of frequency in sample over frequency in ancestor'
x_string  <- 'log2(Prot)'
x_label <- 'Log 2 Protein level'
col_string  <- 'above_14'
# col_label <- 'Number of ribosomes\nper mRNA \n(calulated by TASEP)'
col_label <- 'Protein level\nabove 14'



# get correlation data into data set
fitseq_data_residuals <- get_data_for_correlation_category(fitseq_data_residuals,y_string,y_label,x_string,x_label,col_string,col_label)



fitseq_data_residuals <- filter_n_fitness_residuals_csv(6,fitseq_data_residuals)

fitseq_data_residuals <-fitseq_data_residuals %>% filter(lineage == 'A', above_14 == T)

designs_above_day_12 <- fitseq_data_residuals %>% 
  filter(n_pos_neg == 'Positive', day == 12) %>% 
  ungroup() %>% 
  select(Name)

designs_above_day_other <- fitseq_data_residuals %>% 
  filter(n_pos_neg == 'Positive', day == 16) %>% 
  ungroup() %>% 
  select(Name)

both <-  nrow(intersect(designs_above_day_12, designs_above_day_other))
all_12  <-  nrow(designs_above_day_12) 
all_other  <-  nrow(designs_above_day_other) 

png(paste0('C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\results\\days_comparison\\pos_day_12_16_venn.png'),
    type="cairo",    units="in", width=10, height=10, pointsize=12, res=500)


draw.pairwise.venn(all_12, all_other, both, c("Day 12", "Day 16"))

dev.off()
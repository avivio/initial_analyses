require(dplyr)
require(ineq)


results.dir = 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\results\\'


pre_umi_location = 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\final_match_count_24-04_changed_label.csv'
pre_umi_design_frequency = read.csv(pre_umi_location,check.names=FALSE)

post_umi_location = 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\final_result_24-04_changed_label.csv'
post_umi_design_frequency = read.csv(post_umi_location,check.names=FALSE)

by_evo_time_pre_umi_gini  <- 
  pre_umi_design_frequency %>% 
  select(-design ,anc_1, anc_2, A_4, B_4,C_4,D_4,E_4,F_4, A_8, B_8,C_8,D_8,E_8,F_8
         , A_12, B_12,C_12,D_12,E_12,F_12, A_16, B_16,C_16,D_16,E_16,F_16, A_20, B_20,C_20,D_20,E_20,F_20
         , A_24, B_24,C_24,D_24,E_24,F_24, A_28, B_28,C_28,D_28,E_28,F_28) %>% 
  summarise_each(funs(ineq(.,type="Gini")))


by_evo_time_post_umi_gini  <- 
  post_umi_design_frequency %>% 
  select(-design ,anc_1, anc_2, A_4, B_4,C_4,D_4,E_4,F_4, A_8, B_8,C_8,D_8,E_8,F_8
         , A_12, B_12,C_12,D_12,E_12,F_12, A_16, B_16,C_16,D_16,E_16,F_16, A_20, B_20,C_20,D_20,E_20,F_20
         , A_24, B_24,C_24,D_24,E_24,F_24, A_28, B_28,C_28,D_28,E_28,F_28) %>% 
  summarise_each(funs(ineq(.,type="Gini")))




png(paste0(results.dir,'gini_barplot_by_evo_time_after_relabel.png'),units="in", width=30, height=7, res=100)

sample_names = colnames(select(pre_umi_design_frequency,-design))

barplot(
  t(
    cbind(
  matrix(
    as.numeric(
      by_evo_time_pre_umi_gini
)
),

matrix(
  as.numeric(
    by_evo_time_post_umi_gini
  )
  )
)
)

,names.arg =names(by_evo_time_pre_umi_gini)  , beside = T)

dev.off()





by_line_pre_umi_gini  <- 
  pre_umi_design_frequency %>% 
  select(-design ,anc_1, anc_2, A_4,A_8,A_12,A_16,A_20,A_24,A_28,B_4,B_8,B_12,B_16,B_20,B_24,B_28,
         C_4,C_8,C_12,C_16,C_20,C_24,C_28,D_4,D_8,D_12,D_16,D_20,D_24,D_28,
         E_4,E_8,E_12,E_16,E_20,E_24,E_28,F_4,F_8,F_12,F_16,F_20,F_24,F_28) %>% 
  summarise_each(funs(ineq(.,type="Gini")))


by_line_post_umi_gini  <- 
  post_umi_design_frequency %>% 
  select(-design ,anc_1, anc_2, A_4,A_8,A_12,A_16,A_20,A_24,A_28,B_4,B_8,B_12,B_16,B_20,B_24,B_28,
         C_4,C_8,C_12,C_16,C_20,C_24,C_28,D_4,D_8,D_12,D_16,D_20,D_24,D_28,
         E_4,E_8,E_12,E_16,E_20,E_24,E_28,F_4,F_8,F_12,F_16,F_20,F_24,F_28) %>% 
  summarise_each(funs(ineq(.,type="Gini")))




png(paste0(results.dir,'gini_barplot_by_line_after_relabel.png'),units="in", width=30, height=7, res=100)


barplot(
  t(
    cbind(
      matrix(
        as.numeric(
          by_line_pre_umi_gini
        )
      ),
      
      matrix(
        as.numeric(
          by_line_post_umi_gini
      
      )
    )
  )
  )
  
  ,names.arg =names(by_line_pre_umi_gini)  , beside = T)

dev.off()

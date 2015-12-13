


require(tidyr)
require(dplyr)
bottlenecks <-  read.csv('C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\TASEP_data_only_density.csv')
s_bottlenecks <-  sample_n(bottlenecks,1000)
s_bottlenecks_tidy <- s_bottlenecks %>%  
  gather(position, depth,-Name) %>% 
  separate(position,sep = '\\.',into = c('garbage1','garbage2', 'position'), remove = T, convert = T ) %>% 
  select(-garbage1, -garbage2)

ggplot(s_bottlenecks_tidy, aes(position,Name,fill = depth)) + geom_tile()

require(dplyr)
require(tidyr)
require(ggplot2)

# y.string <-  'TASEP.bottle.neck.depth'
# y.label <- 'Bottle neck depth (calculated by TASEP)'
# x.string  <- 'TASEP.bottle.neck.position'
# x.label <- 'Bottle neck position (calculated by TASEP)'
# col.string  <- 'log2(freq.norm.anc.1)'
# col.label <- 'Log 2 of frequency\nin sample over frequency \nin ancestor'


y.string <-  'log2(freq.norm.anc.1)'
y.label <- 'Log 2 of frequency in sample over frequency in ancestor'
x.string  <- 'log2(Trans)'
x.label <- 'Log 2 Translation efficiency'
col.string  <- 'TASEP.avgRiboNum'
col.label <- 'Number of ribosomes\nper mRNA \n(calulated by TASEP)'


result.dir = 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\results\\fitness_vs_x_and_color\\'
new.dir <- clean.filename(paste(x.string,col.string,sep = '_'))
result.dir <- paste0(result.dir,new.dir,'\\')
dir.create(result.dir)


fitseq.xy <-  fitseq.data.tidy %>%
  mutate_(y = y.string,x = x.string)



fitseq.xy <-  fitseq.xy %>% mutate_(col = col.string)
fitseq.xy <-  fitseq.xy %>% 
  ungroup() %>%
  mutate(col.mid = median(col,na.rm = T))




fitseq.xy  <- fitseq.xy %>% 
  ungroup() %>%
  group_by(day,lineage) %>%
  mutate(all.r = cor.test(x,y)$estimate,all.p = cor.test(x,y)$p.value, 
         all.cor.text = paste0('r = ', round(all.r,2),'\np-value = ' ,formatC(all.p,format='e',digits = 2)))

fitseq.xy <-  fitseq.xy %>% 
  group_by(day,lineage,Promoter.Display,RBS.Display ) %>%
  mutate(group.r = cor.test(x,y)$estimate,group.p = cor.test(x,y)$p.value,
         group.cor.text = paste0('r = ', round(group.r,2),'\np-value = ' ,formatC(group.p,format='e',digits = 2)))


# y.limits <-  fitseq.xy %>%
#   ungroup() %>%
#   summarise(min.y = min(y,na.rm = T),max.y = max(y,na.rm = T))
# y.limits <- c(y.limits$min.y,y.limits$max.y)
# y.limits.all <- y.limits
# y.limits.rbs.promoter<- y.limits

y.limits.all <- c(-5.5,4)
y.limits.rbs.promoter <- c(-4.7,2)

text.loc.y.all <- y.limits.all[2] - ((y.limits.all[2]-y.limits.all[1])/20)

text.loc.y.rbs.promoter <- y.limits.rbs.promoter[2] - ((y.limits.rbs.promoter[2]-y.limits.rbs.promoter[1])/15)


x.limits <-  fitseq.xy %>%
  ungroup() %>%
  summarise(min.x = min(x,na.rm = T),max.x = max(x,na.rm = T))
x.limits <- c(x.limits$min.x,x.limits$max.x)

text.loc.x <- x.limits[1] + ((x.limits[2]-x.limits[1])/100)

fitseq.xy <-  fitseq.xy %>% 
  mutate(text.loc.x = text.loc.x,text.loc.y.all = text.loc.y.all,
         text.loc.y.rbs.promoter = text.loc.y.rbs.promoter )

col.limits <-  fitseq.xy %>%
  ungroup() %>%
  summarise(min.col = min(col,na.rm = T),max.col = max(col,na.rm = T),
            min.quant = quantile(col,probs=.1,na.rm = T),max.quant = quantile(col,probs=.8,na.rm = T))
col.limits <- c(col.limits$min.col,col.limits$max.quant)

col.mid <- mean(col.limits)

# fitseq.xy[fitseq.xy$col > col.limits[2],'col'] <-  col.limits[2]




lineage.letter = 'C'
day.number = 12

for (lineage.letter in c('A','B','C','D','E','F')){
  print(lineage.letter)
  for (day.number in seq(4,28,4)){
  print(day.number)
    col.plots(day.number,lineage.letter)
  }
}

col.plots <- function(day.number,lineage.letter){
  

fitseq.current.sample <-  fitseq.xy %>% filter(day == day.number, lineage == lineage.letter, RBS.Display != 'WT')


filename.all  <- paste('lineage', lineage.letter,'all','day', day.number,'x',x.string ,'vs',y.string, 'col', col.string,sep = '_')
filename.all  <-clean.filename(filename.all)

png(paste0(result.dir,filename.all,'.png'),units="in",  width=15, height=12, res=70)

title <- paste(x.label ,'vs\n',y.label,'\nfor day', day.number, 'lineage', lineage.letter)

p <- ggplot(fitseq.current.sample,aes(x= x, y = y,color  = col,)) +
  geom_point(alpha = 0.3)   +
  theme_minimal() + 
  theme( axis.line = element_line(colour = "black"),
         axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold"),
         strip.text.x = element_text(size=16,face="bold"), strip.text.y = element_text(size=16,face="bold"), 
         plot.title = element_text(size = 25,face = "bold")) +
  # ylim(y.limits) + xlim(x.limits) +
  xlim(x.limits) +
  coord_cartesian(ylim = y.limits.all) +
  ylab(paste0(y.label,'\n')) + 
  xlab(paste0('\n',x.label)) +
  ggtitle(title) + 
  geom_smooth(color = 'red', size = 1.2, method = 'lm') +
   # facet_grid(Promoter.Display ~ RBS.Display)  +
  scale_color_gradient(
    low = 'blue', high = 'red'
    #     low = '#3B4CC0', high = '#B4040D',name=col.label,mid = '#DDDDDD',
    #                         midpoint = col.mid,limits =  col.limits
  ,name =col.label ) +
  geom_text(aes(x=text.loc.x, y=text.loc.y.all, face="bold", inherit.aes=FALSE,
                 parse=FALSE,hjust = 0,label=all.cor.text),colour="red")
print(p)

dev.off()
                        
filename.promoter.rbs  <- paste('lineage', lineage.letter,'promoter_rbs','day', day.number,'x',x.string ,'vs',y.string, 'col', col.string,sep = '_')

filename.promoter.rbs  <-clean.filename(filename.promoter.rbs)

png(paste0(result.dir,filename.promoter.rbs,'.png'),units="in",  width=15, height=12, res=70)
p <- ggplot(fitseq.current.sample,aes(x= x, y = y,color  = col,)) +
  geom_point(alpha = 0.3)   +
  theme_minimal() + 
  theme( axis.line = element_line(colour = "black"),
         axis.text=element_text(size=14), axis.title=element_text(size=16,face="bold"),
         strip.text.x = element_text(size=16,face="bold"), strip.text.y = element_text(size=16,face="bold"), 
         plot.title = element_text(size = 25,face = "bold")) +
  # ylim(y.limits) + xlim(x.limits) +
  xlim(x.limits) +
  coord_cartesian(ylim = y.limits.rbs.promoter) +
  ylab(paste0(y.label,'\n')) + 
  xlab(paste0('\n',x.label)) +
  ggtitle(title) + 
  geom_smooth(color = 'red', size = 1.2, method = 'lm') +
  facet_grid(Promoter.Display ~ RBS.Display)  +
  scale_color_gradient(
    low = 'blue', high = 'red'
#     low = '#3B4CC0', high = '#B4040D',name=col.label,mid = '#DDDDDD',
#                         midpoint = col.mid,limits =  col.limits
,name =col.label ) +
  geom_text(aes(x=text.loc.x, y=1.6, face="bold", inherit.aes=FALSE,
                parse=FALSE,hjust = 0,label=group.cor.text),colour="red")

print(p)
dev.off()
}

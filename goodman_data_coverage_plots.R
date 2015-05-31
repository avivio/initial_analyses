#todo binned barplot method where you give the x and y values and the data
#todo scatter plot thing where you give the data, x and y values

n.terminal.raw = read.csv("C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\n_terminal_data_with_full_seq.csv")

results.dir = 'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\results\\'

generations = 9

#promoter data sets

n.terminal.tai.high.prom = n.terminal.raw[which(n.terminal.raw[,'Promoter.Display'] == 'High'),]

n.terminal.tai.low.prom = n.terminal.raw[which(n.terminal.raw[,'Promoter.Display'] == 'Low'),]


#RBS data sets

n.terminal.tai.strong.rbs = n.terminal.raw[which(n.terminal.raw[,'RBS.Display'] == 'Strong'),]

n.terminal.tai.mid.rbs = n.terminal.raw[which(n.terminal.raw[,'RBS.Display'] == 'Mid'),]

n.terminal.tai.weak.rbs = n.terminal.raw[which(n.terminal.raw[,'RBS.Display'] == 'Weak'),]

n.terminal.tai.wt.rbs = n.terminal.raw[which(n.terminal.raw[,'RBS.Display'] == 'WT'),]






#Promoters


png(paste0(results.dir,'DNA_coverage_by_promoter_unflitered.png'),units="in", width=5, height=8.5, res=100)
boxplot(n.terminal.raw[which(n.terminal.raw[,'Promoter.Display'] == 'High'),'Count.DNA'],
        n.terminal.raw[which(n.terminal.raw[,'Promoter.Display'] == 'Low'),'Count.DNA'], 
        log = 'y',names = c('High', 'Low'), main = 'DNA coverage by promoter strength'
        ,ylab = 'Log DNA coverage',xlab = 'Promoters')
dev.off()

wilcox.test(n.terminal.raw[which(n.terminal.raw[,'Promoter.Display'] == 'High'),'Count.DNA'],
            n.terminal.raw[which(n.terminal.raw[,'Promoter.Display'] == 'Low'),'Count.DNA'],
            alternative = 'less')

median(n.terminal.raw[which(n.terminal.raw[,'Promoter.Display'] == 'High'),'Count.DNA'],na.rm = TRUE)/
  median(n.terminal.raw[which(n.terminal.raw[,'Promoter.Display'] == 'Low'),'Count.DNA'],na.rm = TRUE)

get.growth.rate.difference(mean(n.terminal.tai.high.prom[,'Count.DNA'],na.rm = TRUE)/mean( n.terminal.tai.low.prom[,'Count.DNA'],na.rm = TRUE)  ,generations)


#RBSs

png(paste0(results.dir,'DNA_coverage_by_RBS_unflitered.png'),units="in", width=8, height=8.5, res=100)
boxplot(n.terminal.raw[which(n.terminal.raw[,'RBS.Display'] == 'Strong'),'Count.DNA'],
        n.terminal.raw[which(n.terminal.raw[,'RBS.Display'] == 'Mid'),'Count.DNA'],
        n.terminal.raw[which(n.terminal.raw[,'RBS.Display'] == 'Weak'),'Count.DNA'],
        n.terminal.raw[which(n.terminal.raw[,'RBS.Display'] == 'WT'),'Count.DNA'],
        log = 'y',names = c('Strong', 'Mid','Weak', 'WT'), main = 'DNA coverage by RBS strength'
        ,ylab = 'Log DNA coverage',xlab = 'RBS')
dev.off()



wilcox.test(n.terminal.raw[which(n.terminal.raw[,'RBS.Display'] == 'Strong'),'Count.DNA'],
            n.terminal.raw[which(n.terminal.raw[,'RBS.Display'] == 'Weak'),'Count.DNA'],
            alternative = 'less')

median(n.terminal.raw[which(n.terminal.raw[,'RBS.Display'] == 'Strong'),'Count.DNA'],na.rm = TRUE)/
  median(n.terminal.raw[which(n.terminal.raw[,'RBS.Display'] == 'Weak'),'Count.DNA'],na.rm = TRUE)


get.growth.rate.difference(mean(n.terminal.tai.strong.rbs[,'Count.DNA'],na.rm = TRUE)/mean( n.terminal.tai.weak.rbs[,'Count.DNA'],na.rm = TRUE)  ,generations)


wilcox.test(n.terminal.raw[which(n.terminal.raw[,'RBS.Display'] == 'Mid'),'Count.DNA'],
            n.terminal.raw[which(n.terminal.raw[,'RBS.Display'] == 'Weak'),'Count.DNA'],
            alternative = 'less')

median(n.terminal.raw[which(n.terminal.raw[,'RBS.Display'] == 'Mid'),'Count.DNA'],na.rm = TRUE)/
  median(n.terminal.raw[which(n.terminal.raw[,'RBS.Display'] == 'Weak'),'Count.DNA'],na.rm = TRUE)


get.growth.rate.difference(mean(n.terminal.tai.mid.rbs[,'Count.DNA'],na.rm = TRUE)/mean( n.terminal.tai.weak.rbs[,'Count.DNA'],na.rm = TRUE)  ,generations)


wilcox.test(n.terminal.raw[which(n.terminal.raw[,'RBS.Display'] == 'Strong'),'Count.DNA'],
            n.terminal.raw[which(n.terminal.raw[,'RBS.Display'] == 'Mid'),'Count.DNA'],
            alternative = 'less')

median(n.terminal.raw[which(n.terminal.raw[,'RBS.Display'] == 'Strong'),'Count.DNA'],na.rm = TRUE)/
  median(n.terminal.raw[which(n.terminal.raw[,'RBS.Display'] == 'Mid'),'Count.DNA'],na.rm = TRUE)


get.growth.rate.difference(mean(n.terminal.tai.weak.rbs[,'Count.DNA'],na.rm = TRUE)/mean( n.terminal.tai.mid.rbs[,'Count.DNA'],na.rm = TRUE)  ,generations)


#RBSs and promoters


png(paste0(results.dir,'DNA_coverage_by_RBS_and_promoter_unflitered.png'),units="in", width=20, height=8.5, res=100)
boxplot(
  n.terminal.tai.high.prom[which(n.terminal.raw[,'RBS.Display'] == 'Strong'),'Count.DNA'],
  n.terminal.tai.high.prom[which(n.terminal.raw[,'RBS.Display'] == 'Mid'),'Count.DNA'],
  n.terminal.tai.high.prom[which(n.terminal.raw[,'RBS.Display'] == 'Weak'),'Count.DNA'],
  n.terminal.tai.high.prom[which(n.terminal.raw[,'RBS.Display'] == 'WT'),'Count.DNA'],
  n.terminal.tai.low.prom[which(n.terminal.raw[,'RBS.Display'] == 'Strong'),'Count.DNA'],
  n.terminal.tai.low.prom[which(n.terminal.raw[,'RBS.Display'] == 'Mid'),'Count.DNA'],
  n.terminal.tai.low.prom[which(n.terminal.raw[,'RBS.Display'] == 'Weak'),'Count.DNA'],
  n.terminal.tai.low.prom[which(n.terminal.raw[,'RBS.Display'] == 'WT'),'Count.DNA'],    
  names = c('High Promoter Strong RBS', 'High Promoter Mid RBS','High Promoter Weak RBS', 'High Promoter WT RBS'
            ,'Low Promoter Strong RBS', 'Low Promoter Mid RBS','Low Promoter Weak RBS', 'Low Promoter WT RBS'),
  log = 'y',main = 'DNA coverage by promoter and RBS strength',ylab = 'Log DNA coverage',xlab = 'RBS and Promoter')
dev.off()


wilcox.test(n.terminal.tai.high.prom[which(n.terminal.tai.high.prom[,'RBS.Display'] == 'Strong'),'Count.DNA'],
            n.terminal.tai.high.prom[which(n.terminal.tai.high.prom[,'RBS.Display'] == 'Weak'),'Count.DNA'],
            alternative = 'less')

median(n.terminal.tai.high.prom[which(n.terminal.tai.high.prom[,'RBS.Display'] == 'Strong'),'Count.DNA'],na.rm = TRUE)/
  median(n.terminal.tai.high.prom[which(n.terminal.tai.high.prom[,'RBS.Display'] == 'Weak'),'Count.DNA'],na.rm = TRUE)

get.growth.rate.difference(mean(n.terminal.tai.high.prom[which(n.terminal.tai.high.prom[,'RBS.Display'] == 'Strong'),'Count.DNA'],na.rm = TRUE)/
                             mean( n.terminal.tai.high.prom[which(n.terminal.tai.high.prom[,'RBS.Display'] == 'Weak'),'Count.DNA'],na.rm = TRUE)  ,generations)


wilcox.test(n.terminal.tai.high.prom[which(n.terminal.tai.high.prom[,'RBS.Display'] == 'Strong'),'Count.DNA'],
            n.terminal.tai.high.prom[which(n.terminal.tai.high.prom[,'RBS.Display'] == 'Mid'),'Count.DNA'],
            alternative = 'less')

median(n.terminal.tai.high.prom[which(n.terminal.tai.high.prom[,'RBS.Display'] == 'Strong'),'Count.DNA'],na.rm = TRUE)/
  median(n.terminal.tai.high.prom[which(n.terminal.tai.high.prom[,'RBS.Display'] == 'Mid'),'Count.DNA'],na.rm = TRUE)

get.growth.rate.difference(mean(n.terminal.tai.high.prom[which(n.terminal.tai.high.prom[,'RBS.Display'] == 'Strong'),'Count.DNA'],na.rm = TRUE)/
                             mean( n.terminal.tai.high.prom[which(n.terminal.tai.high.prom[,'RBS.Display'] == 'Mid'),'Count.DNA'],na.rm = TRUE)  ,generations)

wilcox.test(n.terminal.tai.high.prom[which(n.terminal.tai.high.prom[,'RBS.Display'] == 'Mid'),'Count.DNA'],
            n.terminal.tai.high.prom[which(n.terminal.tai.high.prom[,'RBS.Display'] == 'Weak'),'Count.DNA'],
            alternative = 'less')

median(n.terminal.tai.high.prom[which(n.terminal.tai.high.prom[,'RBS.Display'] == 'Mid'),'Count.DNA'],na.rm = TRUE)/
  median(n.terminal.tai.high.prom[which(n.terminal.tai.high.prom[,'RBS.Display'] == 'Weak'),'Count.DNA'],na.rm = TRUE)

get.growth.rate.difference(mean(n.terminal.tai.high.prom[which(n.terminal.tai.high.prom[,'RBS.Display'] == 'Mid'),'Count.DNA'],na.rm = TRUE)/
                             mean( n.terminal.tai.high.prom[which(n.terminal.tai.high.prom[,'RBS.Display'] == 'Weak'),'Count.DNA'],na.rm = TRUE)  ,generations)


wilcox.test(n.terminal.tai.high.prom[which(n.terminal.tai.high.prom[,'RBS.Display'] == 'Strong'),'Count.DNA'],
            n.terminal.tai.low.prom[which(n.terminal.tai.low.prom[,'RBS.Display'] == 'Strong'),'Count.DNA'],
            alternative = 'less')

median(n.terminal.tai.high.prom[which(n.terminal.tai.high.prom[,'RBS.Display'] == 'Strong'),'Count.DNA'],na.rm = TRUE)/
  median(n.terminal.tai.low.prom[which(n.terminal.tai.low.prom[,'RBS.Display'] == 'Strong'),'Count.DNA'],na.rm = TRUE)

wilcox.test(n.terminal.tai.low.prom[which(n.terminal.tai.low.prom[,'RBS.Display'] == 'Strong'),'Count.DNA'],
            n.terminal.tai.low.prom[which(n.terminal.tai.low.prom[,'RBS.Display'] == 'Weak'),'Count.DNA'],
            alternative = 'less')

median(n.terminal.tai.low.prom[which(n.terminal.tai.low.prom[,'RBS.Display'] == 'Strong'),'Count.DNA'],na.rm = TRUE)/
  median(n.terminal.tai.low.prom[which(n.terminal.tai.low.prom[,'RBS.Display'] == 'Weak'),'Count.DNA'],na.rm = TRUE)


#cds types coverage
cds.types = sort(unique(as.character(n.terminal.raw[,'CDS.type'])))
dna.cds.types = vector("list", length(cds.types))

names(dna.cds.types) = gsub('גˆ†','delta',cds.types)
for (i in 1:length(cds.types)){
  dna.cds.types[[i]] = n.terminal.raw[which(n.terminal.raw[,'CDS.type'] == cds.types[i]),'Count.DNA']
}
dna.cds.types = c(dna.cds.types[-5],dna.cds.types[5])
names(dna.cds.types)[0:2] = c('Rare', 'Common')

png(paste0(results.dir,'DNA_coverage_by_cds_type_unflitered.png'),units="in", width=15, height=8.5, res=100)

boxplot(dna.cds.types, log = 'y',main = 'DNA coverage by CDS type', ylab = 'Log DNA coverage',xlab = 'CDS type')
dev.off()


wilcox.test(n.terminal.raw[which(n.terminal.raw[,'CDS.type'] == 'Max Rare'),'Count.DNA'],
            n.terminal.raw[which(n.terminal.raw[,'CDS.type'] == 'Min Rare'),'Count.DNA'],
            alternative = 'less')

1-median(n.terminal.raw[which(n.terminal.raw[,'CDS.type'] == 'Max Rare'),'Count.DNA'],na.rm = TRUE)/
  median(n.terminal.raw[which(n.terminal.raw[,'CDS.type'] == 'Min Rare'),'Count.DNA'],na.rm = TRUE)

get.growth.rate.difference(mean(n.terminal.raw[which(n.terminal.raw[,'CDS.type'] == 'Max Rare'),'Count.DNA'],na.rm = TRUE)/
                             mean(n.terminal.raw[which(n.terminal.raw[,'CDS.type'] == 'Min Rare'),'Count.DNA'],na.rm = TRUE)  ,generations)


wilcox.test(n.terminal.raw[which(n.terminal.raw[,'CDS.type'] == "גˆ†G 10"),'Count.DNA'],
            n.terminal.raw[which(n.terminal.raw[,'CDS.type'] == 'גˆ†G 1'),'Count.DNA'],
            alternative = 'less')

1-median(n.terminal.raw[which(n.terminal.raw[,'CDS.type'] == 'גˆ†G 10'),'Count.DNA'],na.rm = TRUE)/
  median(n.terminal.raw[which(n.terminal.raw[,'CDS.type'] == 'גˆ†G 1'),'Count.DNA'],na.rm = TRUE)

get.growth.rate.difference(mean(n.terminal.raw[which(n.terminal.raw[,'CDS.type'] == 'Max Rare'),'Count.DNA'],na.rm = TRUE)/
                             mean(n.terminal.raw[which(n.terminal.raw[,'CDS.type'] == 'Min Rare'),'Count.DNA'],na.rm = TRUE)  ,generations)


#tai

tai.hist = hist(n.terminal.raw[,'tAI'],xlim = c(max(n.terminal.raw[,'tAI']),min(n.terminal.raw[,'tAI']))
                , xlab = 'Average TAI content of variable region', main = 'Histogram of average TAI content of variable region')

n.term.cut.tai = cut(n.terminal.raw[,'tAI'],tai.hist$breaks, include.lowest = TRUE)
n.term.dna.sum.tai = vector(length = length(levels(n.term.cut.tai)),mode = 'numeric')

names(n.term.dna.sum.tai) = levels(n.term.cut.tai)

for (i in 1:length(n.term.cut.tai)){
 
  if (is.na(n.terminal.raw[,'Count.DNA'][i])){
    current = 0
  } else {
    current = n.terminal.raw[,'Count.DNA'][i]
  }
  
  n.term.dna.sum.tai[n.term.cut.tai[i]] = n.term.dna.sum.tai[n.term.cut.tai[i]] + current

  
}

png(paste0(results.dir,'DNA_coverage_per_tai_level_unflitered.png'),units="in", width=15, height=8.5, res=100)

barplot(n.term.dna.sum.tai/tai.hist$counts,ylim = c(0,3100), ylab = 'DNA coverage normalized by number of variants at TAI level', xlab = 'TAI',
        main = 'DNA coverage per level of TAI normalized to ammount of varaints at TAI level')
dev.off()



png(paste0(results.dir,'DNA_coverage_vs_tai_by_promoter_unflitered.png'),units="in", width=15, height=15, res=100)
par(mfrow = c(2,1))

plot(n.terminal.tai.high.prom[,'tAI'],n.terminal.tai.high.prom[,'Count.DNA'],log = 'y'
         ,main = 'High Promoter',ylab = 'Log DNA Coverage', xlab = 'Average TAI of variable region'
         ,xlim = c(min(n.terminal.raw[,'tAI']),max(n.terminal.raw[,'tAI']))
     ,ylim = c(min(n.terminal.raw[,'Count.DNA'],na.rm = TRUE),max(n.terminal.raw[,'Count.DNA'],na.rm = TRUE))
     ,col="#00000080")
    
    plot(n.terminal.tai.low.prom[,'tAI'],n.terminal.tai.low.prom[,'Count.DNA'],log = 'y'
         ,main = 'Low Promoter',ylab = 'Log DNA Coverage', xlab = 'Average TAI of variable region'
         ,xlim = c(min(n.terminal.raw[,'tAI']),max(n.terminal.raw[,'tAI']))
         ,ylim = c(min(n.terminal.raw[,'Count.DNA'],na.rm = TRUE),max(n.terminal.raw[,'Count.DNA'],na.rm = TRUE))
         ,col="#00000080")
title("Log DNA coverage vs TAI for different promoters", outer = TRUE, line = -1)

    
dev.off()


png(paste0(results.dir,'DNA_coverage_vs_tai_by_RBS_unflitered.png'),units="in", width=20, height=15, res=100)
par(mfrow = c(2,2))
plot(n.terminal.tai.strong.rbs[,'tAI'],n.terminal.tai.strong.rbs[,'Count.DNA'],log = 'y'
     ,main = 'Strong RBS',ylab = 'Log DNA Coverage', xlab = 'Average TAI of variable region'
     ,xlim = c(min(n.terminal.raw[,'tAI']),max(n.terminal.raw[,'tAI']))
     ,ylim = c(min(n.terminal.raw[,'Count.DNA'],na.rm = TRUE),max(n.terminal.raw[,'Count.DNA'],na.rm = TRUE))
     ,col="#00000080")

plot(n.terminal.tai.mid.rbs[,'tAI'],n.terminal.tai.mid.rbs[,'Count.DNA'],log = 'y'
     ,main = 'Mid RBS',ylab = 'Log DNA Coverage', xlab = 'Average TAI of variable region'
     ,xlim = c(min(n.terminal.raw[,'tAI']),max(n.terminal.raw[,'tAI']))
     ,ylim = c(min(n.terminal.raw[,'Count.DNA'],na.rm = TRUE),max(n.terminal.raw[,'Count.DNA'],na.rm = TRUE))
     ,col="#00000080")

plot(n.terminal.tai.weak.rbs[,'tAI'],n.terminal.tai.weak.rbs[,'Count.DNA'],log = 'y'
     ,main = 'Weak RBS',ylab = 'Log DNA Coverage', xlab = 'Average TAI of variable region'
     ,xlim = c(min(n.terminal.raw[,'tAI']),max(n.terminal.raw[,'tAI']))
     ,ylim = c(min(n.terminal.raw[,'Count.DNA'],na.rm = TRUE),max(n.terminal.raw[,'Count.DNA'],na.rm = TRUE))
     ,col="#00000080")

plot(n.terminal.tai.wt.rbs[,'tAI'],n.terminal.tai.wt.rbs[,'Count.DNA'],log = 'y'
     ,main = 'WT RBS',ylab = 'Log DNA Coverage', xlab = 'Average TAI of variable region'
     ,xlim = c(min(n.terminal.raw[,'tAI']),max(n.terminal.raw[,'tAI']))
     ,ylim = c(min(n.terminal.raw[,'Count.DNA'],na.rm = TRUE),max(n.terminal.raw[,'Count.DNA'],na.rm = TRUE))
     ,col="#00000080")

title("Log DNA coverage vs TAI for different RBSs", outer = TRUE, line = -1)

dev.off()



png(paste0(results.dir,'DNA_coverage_vs_tai_by_RBS_and_promoter_unflitered.png'),units="in", width=20, height=15, res=100)
par(mfrow = c(4,2))
plot(n.terminal.tai.high.prom[which(n.terminal.raw[,'RBS.Display'] == 'Strong'),'tAI'],n.terminal.tai.high.prom[which(n.terminal.raw[,'RBS.Display'] == 'Strong'),'Count.DNA'],log = 'y'
     ,main = 'Strong RBS and High Promoter',ylab = 'Log DNA Coverage', xlab = 'Average TAI of variable region'
     ,xlim = c(min(n.terminal.raw[,'tAI']),max(n.terminal.raw[,'tAI']))
     ,ylim = c(min(n.terminal.raw[,'Count.DNA'],na.rm = TRUE),max(n.terminal.raw[,'Count.DNA'],na.rm = TRUE))
     ,col="#00000080")
plot(n.terminal.tai.high.prom[which(n.terminal.raw[,'RBS.Display'] == 'Mid'),'tAI'],n.terminal.tai.high.prom[which(n.terminal.raw[,'RBS.Display'] == 'Mid'),'Count.DNA'],log = 'y'
     ,main = 'Mid RBS and High Promoter',ylab = 'Log DNA Coverage', xlab = 'Average TAI of variable region'
     ,xlim = c(min(n.terminal.raw[,'tAI']),max(n.terminal.raw[,'tAI']))
     ,ylim = c(min(n.terminal.raw[,'Count.DNA'],na.rm = TRUE),max(n.terminal.raw[,'Count.DNA'],na.rm = TRUE))
     ,col="#00000080")
plot(n.terminal.tai.high.prom[which(n.terminal.raw[,'RBS.Display'] == 'Weak'),'tAI'],n.terminal.tai.high.prom[which(n.terminal.raw[,'RBS.Display'] == 'Weak'),'Count.DNA'],log = 'y'
     ,main = 'Weak RBS and High Promoter',ylab = 'Log DNA Coverage', xlab = 'Average TAI of variable region'
     ,xlim = c(min(n.terminal.raw[,'tAI']),max(n.terminal.raw[,'tAI']))
     ,ylim = c(min(n.terminal.raw[,'Count.DNA'],na.rm = TRUE),max(n.terminal.raw[,'Count.DNA'],na.rm = TRUE))
     ,col="#00000080")
plot(n.terminal.tai.high.prom[which(n.terminal.raw[,'RBS.Display'] == 'WT'),'tAI'],n.terminal.tai.high.prom[which(n.terminal.raw[,'RBS.Display'] == 'WT'),'Count.DNA'],log = 'y'
     ,main = 'WT RBS and High Promoter',ylab = 'Log DNA Coverage', xlab = 'Average TAI of variable region'
     ,xlim = c(min(n.terminal.raw[,'tAI']),max(n.terminal.raw[,'tAI']))
     ,ylim = c(min(n.terminal.raw[,'Count.DNA'],na.rm = TRUE),max(n.terminal.raw[,'Count.DNA'],na.rm = TRUE))
     ,col="#00000080")
plot(n.terminal.tai.low.prom[which(n.terminal.raw[,'RBS.Display'] == 'Strong'),'tAI'],n.terminal.tai.low.prom[which(n.terminal.raw[,'RBS.Display'] == 'Strong'),'Count.DNA'],log = 'y'
     ,main = 'Strong RBS and Low Promoter',ylab = 'Log DNA Coverage', xlab = 'Average TAI of variable region'
     ,xlim = c(min(n.terminal.raw[,'tAI']),max(n.terminal.raw[,'tAI']))
     ,ylim = c(min(n.terminal.raw[,'Count.DNA'],na.rm = TRUE),max(n.terminal.raw[,'Count.DNA'],na.rm = TRUE))
     ,col="#00000080")
plot(n.terminal.tai.low.prom[which(n.terminal.raw[,'RBS.Display'] == 'Mid'),'tAI'],n.terminal.tai.low.prom[which(n.terminal.raw[,'RBS.Display'] == 'Mid'),'Count.DNA'],log = 'y'
     ,main = 'Mid RBS and Low Promoter',ylab = 'Log DNA Coverage', xlab = 'Average TAI of variable region'
     ,xlim = c(min(n.terminal.raw[,'tAI']),max(n.terminal.raw[,'tAI']))
     ,ylim = c(min(n.terminal.raw[,'Count.DNA'],na.rm = TRUE),max(n.terminal.raw[,'Count.DNA'],na.rm = TRUE))
     ,col="#00000080")
plot(n.terminal.tai.low.prom[which(n.terminal.raw[,'RBS.Display'] == 'Weak'),'tAI'],n.terminal.tai.low.prom[which(n.terminal.raw[,'RBS.Display'] == 'Weak'),'Count.DNA'],log = 'y'
     ,main = 'Weak RBS and Low Promoter',ylab = 'Log DNA Coverage', xlab = 'Average TAI of variable region'
     ,xlim = c(min(n.terminal.raw[,'tAI']),max(n.terminal.raw[,'tAI']))
     ,ylim = c(min(n.terminal.raw[,'Count.DNA'],na.rm = TRUE),max(n.terminal.raw[,'Count.DNA'],na.rm = TRUE))
     ,col="#00000080")
plot(n.terminal.tai.low.prom[which(n.terminal.raw[,'RBS.Display'] == 'WT'),'tAI'],n.terminal.tai.low.prom[which(n.terminal.raw[,'RBS.Display'] == 'WT'),'Count.DNA'],log = 'y'
     ,main = 'WT RBS and Low Promoter',ylab = 'Log DNA Coverage', xlab = 'Average TAI of variable region'
     ,xlim = c(min(n.terminal.raw[,'tAI']),max(n.terminal.raw[,'tAI']))
     ,ylim = c(min(n.terminal.raw[,'Count.DNA'],na.rm = TRUE),max(n.terminal.raw[,'Count.DNA'],na.rm = TRUE))
     ,col="#00000080")
title("Log DNA coverage vs TAI for different RBS promoter pairs", outer = TRUE, line = -1)

dev.off()


png(paste0(results.dir,'DNA_coverage_vs_tai_unflitered.png'),units="in", width=15, height=8.5, res=100)


plot(n.terminal.raw[,'tAI'],n.terminal.raw[,'Count.DNA'],log = 'y'
     ,main = 'Log DNA coverage vs TAI',ylab = 'LogDNA Coverage', xlab = 'Average TAI of variable region'
     ,xlim = c(min(n.terminal.raw[,'tAI']),max(n.terminal.raw[,'tAI']))
     ,ylim = c(min(n.terminal.raw[,'Count.DNA'],na.rm = TRUE),max(n.terminal.raw[,'Count.DNA'],na.rm = TRUE))
     ,col="#00000080")

dev.off()


#gc

gc.hist = hist(n.terminal.raw[,'GC']
                , xlab = 'Average GC content of variable region', main = 'Histogram of average GC content of variable region', right = TRUE)

n.term.cut.gc = cut(n.terminal.raw[,'GC'],gc.hist$breaks, include.lowest = TRUE)
n.term.dna.sum.gc = vector(length = length(levels(n.term.cut.gc)),mode = 'numeric')

names(n.term.dna.sum.gc) = levels(n.term.cut.gc)

for (i in 1:length(n.term.cut.gc)){
  
  if (is.na(n.terminal.raw[,'Count.DNA'][i])){
    current = 0
  } else {
    current = n.terminal.raw[,'Count.DNA'][i]
  }
  
  n.term.dna.sum.gc[n.term.cut.gc[i]] = n.term.dna.sum.gc[n.term.cut.gc[i]] + current
  
  
}

png(paste0(results.dir,'DNA_coverage_per_GC_level_unflitered.png'),units="in", width=15, height=8.5, res=100)

barplot(n.term.dna.sum.gc/gc.hist$counts, ylab = 'DNA coverage normalized by number of variants at GC content level', xlab = 'GC content',
        main = 'DNA coverage per level of GC content normalized to ammount of varaints at GC content level')
dev.off()

cor(n.term.dna.sum.gc[1:18]/gc.hist$counts[1:18],c(1:18), method ='spearman')
cor(n.term.dna.sum.gc[18:length(n.term.dna.sum.gc)]/gc.hist$counts[18:length(n.term.dna.sum.gc)],c(18:length(n.term.dna.sum.gc)), method ='spearman')



png(paste0(results.dir,'DNA_coverage_vs_GC_by_promoter_unflitered.png'),units="in", width=15, height=15, res=100)
par(mfrow = c(2,1))
plot(n.terminal.tai.high.prom[,'GC'],n.terminal.tai.high.prom[,'Count.DNA'],log = 'y'
     ,main = 'High Promoter',ylab = 'Log DNA Coverage', xlab = 'GC content of variable region'
     ,xlim = c(min(n.terminal.raw[,'GC']),max(n.terminal.raw[,'GC']))
     ,ylim = c(min(n.terminal.raw[,'Count.DNA'],na.rm = TRUE),max(n.terminal.raw[,'Count.DNA'],na.rm = TRUE))
     ,col="#00000080")

plot(n.terminal.tai.low.prom[,'GC'],n.terminal.tai.low.prom[,'Count.DNA'],log = 'y'
     ,main = 'Low Promoter',ylab = 'Log DNA Coverage', xlab = 'GC content of variable region'
     ,xlim = c(min(n.terminal.raw[,'GC']),max(n.terminal.raw[,'GC']))
     ,ylim = c(min(n.terminal.raw[,'Count.DNA'],na.rm = TRUE),max(n.terminal.raw[,'Count.DNA'],na.rm = TRUE))
     ,col="#00000080")
title("Log DNA coverage vs GC content for different promoters", outer = TRUE, line = -1)

dev.off()

cor(n.terminal.tai.high.prom[,'GC'],n.terminal.tai.high.prom[,'Count.DNA'], method ='spearman',
    use = 'pairwise.complete.obs')
cor(n.terminal.tai.low.prom[,'GC'],n.terminal.tai.low.prom[,'Count.DNA'], method ='spearman',
    use = 'pairwise.complete.obs')



png(paste0(results.dir,'DNA_coverage_vs_GC_by_RBS_unflitered.png'),units="in", width=20, height=15, res=100)
par(mfrow = c(2,2))
plot(n.terminal.tai.strong.rbs[,'GC'],n.terminal.tai.strong.rbs[,'Count.DNA'],log = 'y'
     ,main = 'Strong RBS',ylab = 'Log DNA Coverage', xlab = 'Gc content of variable region'
     ,xlim = c(min(n.terminal.raw[,'GC']),max(n.terminal.raw[,'GC']))
     ,ylim = c(min(n.terminal.raw[,'Count.DNA'],na.rm = TRUE),max(n.terminal.raw[,'Count.DNA'],na.rm = TRUE))
     ,col="#00000080")

plot(n.terminal.tai.mid.rbs[,'GC'],n.terminal.tai.mid.rbs[,'Count.DNA'],log = 'y'
     ,main = 'Mid RBS',ylab = 'Log DNA Coverage', xlab = 'GC content of variable region'
     ,xlim = c(min(n.terminal.raw[,'GC']),max(n.terminal.raw[,'GC']))
     ,ylim = c(min(n.terminal.raw[,'Count.DNA'],na.rm = TRUE),max(n.terminal.raw[,'Count.DNA'],na.rm = TRUE))
     ,col="#00000080")

plot(n.terminal.tai.weak.rbs[,'GC'],n.terminal.tai.weak.rbs[,'Count.DNA'],log = 'y'
     ,main = 'Weak RBS',ylab = 'Log DNA Coverage', xlab = 'GC content of variable region'
     ,xlim = c(min(n.terminal.raw[,'GC']),max(n.terminal.raw[,'GC']))
     ,ylim = c(min(n.terminal.raw[,'Count.DNA'],na.rm = TRUE),max(n.terminal.raw[,'Count.DNA'],na.rm = TRUE))
     ,col="#00000080")

plot(n.terminal.tai.wt.rbs[,'GC'],n.terminal.tai.wt.rbs[,'Count.DNA'],log = 'y'
     ,main = 'WT RBS',ylab = 'Log DNA Coverage', xlab = 'GC content of variable region'
     ,xlim = c(min(n.terminal.raw[,'GC']),max(n.terminal.raw[,'GC']))
     ,ylim = c(min(n.terminal.raw[,'Count.DNA'],na.rm = TRUE),max(n.terminal.raw[,'Count.DNA'],na.rm = TRUE))
     ,col="#00000080")

title("Log DNA coverage vs GC content for different promoters", outer = TRUE, line = -1)

dev.off()

cor(n.terminal.tai.strong.rbs[,'GC'],n.terminal.tai.strong.rbs[,'Count.DNA'], method ='spearman',
    use = 'pairwise.complete.obs')
cor(n.terminal.tai.mid.rbs[,'GC'],n.terminal.tai.mid.rbs[,'Count.DNA'], method ='spearman',
    use = 'pairwise.complete.obs')
cor(n.terminal.tai.weak.rbs[,'GC'],n.terminal.tai.weak.rbs[,'Count.DNA'], method ='spearman',
    use = 'pairwise.complete.obs')
cor(n.terminal.tai.wt.rbs[,'GC'],n.terminal.tai.wt.rbs[,'Count.DNA'], method ='spearman',
    use = 'pairwise.complete.obs')


png(paste0(results.dir,'DNA_coverage_vs_gc_by_RBS_and_promoter_unflitered.png'),units="in", width=20, height=15, res=100)
par(mfrow = c(4,2))


plot(n.terminal.tai.high.prom[which(n.terminal.tai.high.prom[,'RBS.Display'] == 'Strong'),'GC'],n.terminal.tai.high.prom[which(n.terminal.tai.high.prom[,'RBS.Display'] == 'Strong'),'Count.DNA'],log = 'y'
     ,main = 'Strong RBS and High Promoter',ylab = 'Log DNA Coverage', xlab = 'Average GC content of variable region'
     ,xlim = c(min(n.terminal.raw[,'GC']),max(n.terminal.raw[,'GC']))
     ,ylim = c(min(n.terminal.raw[,'Count.DNA'],na.rm = TRUE),max(n.terminal.raw[,'Count.DNA'],na.rm = TRUE))
     ,col="#00000080")
plot(n.terminal.tai.high.prom[which(n.terminal.tai.high.prom[,'RBS.Display'] == 'Mid'),'GC'],n.terminal.tai.high.prom[which(n.terminal.tai.high.prom[,'RBS.Display'] == 'Mid'),'Count.DNA'],log = 'y'
     ,main = 'Mid RBS and High Promoter',ylab = 'Log DNA Coverage', xlab = 'Average GC content of variable region'
     ,xlim = c(min(n.terminal.raw[,'GC']),max(n.terminal.raw[,'GC']))
     ,ylim = c(min(n.terminal.raw[,'Count.DNA'],na.rm = TRUE),max(n.terminal.raw[,'Count.DNA'],na.rm = TRUE))
     ,col="#00000080")
plot(n.terminal.tai.high.prom[which(n.terminal.tai.high.prom[,'RBS.Display'] == 'Weak'),'GC'],n.terminal.tai.high.prom[which(n.terminal.tai.high.prom[,'RBS.Display'] == 'Weak'),'Count.DNA'],log = 'y'
     ,main = 'Weak RBS and High Promoter',ylab = 'Log DNA Coverage', xlab = 'Average GC content of variable region'
     ,xlim = c(min(n.terminal.raw[,'GC']),max(n.terminal.raw[,'GC']))
     ,ylim = c(min(n.terminal.raw[,'Count.DNA'],na.rm = TRUE),max(n.terminal.raw[,'Count.DNA'],na.rm = TRUE))
     ,col="#00000080")
plot(n.terminal.tai.high.prom[which(n.terminal.tai.high.prom[,'RBS.Display'] == 'WT'),'GC'],n.terminal.tai.high.prom[which(n.terminal.tai.high.prom[,'RBS.Display'] == 'WT'),'Count.DNA'],log = 'y'
     ,main = 'WT RBS and High Promoter',ylab = 'Log DNA Coverage', xlab = 'Average GC content of variable region'
     ,xlim = c(min(n.terminal.raw[,'GC']),max(n.terminal.raw[,'GC']))
     ,ylim = c(min(n.terminal.raw[,'Count.DNA'],na.rm = TRUE),max(n.terminal.raw[,'Count.DNA'],na.rm = TRUE))
     ,col="#00000080")
plot(n.terminal.tai.low.prom[which(n.terminal.tai.low.prom[,'RBS.Display'] == 'Strong'),'GC'],n.terminal.tai.low.prom[which(n.terminal.tai.low.prom[,'RBS.Display'] == 'Strong'),'Count.DNA'],log = 'y'
     ,main = 'Strong RBS and Low Promoter',ylab = 'Log DNA Coverage', xlab = 'Average GC content of variable region'
     ,xlim = c(min(n.terminal.raw[,'GC']),max(n.terminal.raw[,'GC']))
     ,ylim = c(min(n.terminal.raw[,'Count.DNA'],na.rm = TRUE),max(n.terminal.raw[,'Count.DNA'],na.rm = TRUE))
     ,col="#00000080")
plot(n.terminal.tai.low.prom[which(n.terminal.tai.low.prom[,'RBS.Display'] == 'Mid'),'GC'],n.terminal.tai.low.prom[which(n.terminal.tai.low.prom[,'RBS.Display'] == 'Mid'),'Count.DNA'],log = 'y'
     ,main = 'Mid RBS and Low Promoter',ylab = 'Log DNA Coverage', xlab = 'Average GC content of variable region'
     ,xlim = c(min(n.terminal.raw[,'GC']),max(n.terminal.raw[,'GC']))
     ,ylim = c(min(n.terminal.raw[,'Count.DNA'],na.rm = TRUE),max(n.terminal.raw[,'Count.DNA'],na.rm = TRUE))
     ,col="#00000080")
plot(n.terminal.tai.low.prom[which(n.terminal.tai.low.prom[,'RBS.Display'] == 'Weak'),'GC'],n.terminal.tai.low.prom[which(n.terminal.tai.low.prom[,'RBS.Display'] == 'Weak'),'Count.DNA'],log = 'y'
     ,main = 'Weak RBS and Low Promoter',ylab = 'Log DNA Coverage', xlab = 'Average GC content of variable region'
     ,xlim = c(min(n.terminal.raw[,'GC']),max(n.terminal.raw[,'GC']))
     ,ylim = c(min(n.terminal.raw[,'Count.DNA'],na.rm = TRUE),max(n.terminal.raw[,'Count.DNA'],na.rm = TRUE))
     ,col="#00000080")
plot(n.terminal.tai.low.prom[which(n.terminal.tai.low.prom[,'RBS.Display'] == 'WT'),'GC'],n.terminal.tai.low.prom[which(n.terminal.tai.low.prom[,'RBS.Display'] == 'WT'),'Count.DNA'],log = 'y'
     ,main = 'WT RBS and Low Promoter',ylab = 'Log DNA Coverage', xlab = 'Average GC content of variable region'
     ,xlim = c(min(n.terminal.raw[,'GC']),max(n.terminal.raw[,'GC']))
     ,ylim = c(min(n.terminal.raw[,'Count.DNA'],na.rm = TRUE),max(n.terminal.raw[,'Count.DNA'],na.rm = TRUE))
     ,col="#00000080")
title("Log DNA coverage vs GC content for different RBS promoter pairs", outer = TRUE, line = -1)

dev.off()



cor(n.terminal.tai.high.prom[which(n.terminal.tai.high.prom[,'RBS.Display'] == 'Strong'),'GC'],
    n.terminal.tai.high.prom[which(n.terminal.tai.high.prom[,'RBS.Display'] == 'Strong'),'Count.DNA'], 
    method ='spearman',
    use = 'pairwise.complete.obs')
cor(n.terminal.tai.high.prom[which(n.terminal.tai.high.prom[,'RBS.Display'] == 'Mid'),'GC'],
    n.terminal.tai.high.prom[which(n.terminal.tai.high.prom[,'RBS.Display'] == 'Mid'),'Count.DNA'], 
    method ='spearman',
    use = 'pairwise.complete.obs')
cor(n.terminal.tai.high.prom[which(n.terminal.tai.high.prom[,'RBS.Display'] == 'Weak'),'GC'],
    n.terminal.tai.high.prom[which(n.terminal.tai.high.prom[,'RBS.Display'] == 'Weak'),'Count.DNA'], 
    method ='spearman',
    use = 'pairwise.complete.obs')
cor(n.terminal.tai.high.prom[which(n.terminal.tai.high.prom[,'RBS.Display'] == 'WT'),'GC'],
    n.terminal.tai.high.prom[which(n.terminal.tai.high.prom[,'RBS.Display'] == 'WT'),'Count.DNA'], 
    method ='spearman',
    use = 'pairwise.complete.obs')
cor(n.terminal.tai.low.prom[which(n.terminal.tai.low.prom[,'RBS.Display'] == 'Strong'),'GC'],
    n.terminal.tai.low.prom[which(n.terminal.tai.low.prom[,'RBS.Display'] == 'Strong'),'Count.DNA'], 
    method ='spearman',
    use = 'pairwise.complete.obs')
cor(n.terminal.tai.low.prom[which(n.terminal.tai.low.prom[,'RBS.Display'] == 'Mid'),'GC'],
    n.terminal.tai.low.prom[which(n.terminal.tai.low.prom[,'RBS.Display'] == 'Mid'),'Count.DNA'], 
    method ='spearman',
    use = 'pairwise.complete.obs')
cor(n.terminal.tai.low.prom[which(n.terminal.tai.low.prom[,'RBS.Display'] == 'Weak'),'GC'],
    n.terminal.tai.low.prom[which(n.terminal.tai.low.prom[,'RBS.Display'] == 'Weak'),'Count.DNA'], 
    method ='spearman',
    use = 'pairwise.complete.obs')
cor(n.terminal.tai.low.prom[which(n.terminal.tai.low.prom[,'RBS.Display'] == 'WT'),'GC'],
    n.terminal.tai.low.prom[which(n.terminal.tai.low.prom[,'RBS.Display'] == 'WT'),'Count.DNA'], 
    method ='spearman',
    use = 'pairwise.complete.obs')


png(paste0(results.dir,'DNA_coverage_vs_gc__unflitered.png'),units="in", width=11, height=8.5, res=100)
plot(n.terminal.raw[,'GC'],n.terminal.raw[,'Count.DNA'],log = 'y'
     ,main = 'Log DNA coverage vs GC content',ylab = 'Log DNA Coverage',  xlab = 'Average GC content of variable region'    ,
     ,xlim = c(min(n.terminal.raw[,'GC']),max(n.terminal.raw[,'GC']))
     ,ylim = c(min(n.terminal.raw[,'Count.DNA'],na.rm = TRUE),max(n.terminal.raw[,'Count.DNA'],na.rm = TRUE))
     ,col="#00000080")
dev.off()

#deltag

deltag.hist = hist(n.terminal.raw[,'dG']
               , xlab = 'Average delta G content of variable region', main = 'Histogram of average delta G content of variable region', right = TRUE)

n.term.cut.deltag = cut(n.terminal.raw[,'dG'],deltag.hist$breaks, include.lowest = TRUE)
n.term.dna.sum.deltag = vector(length = length(levels(n.term.cut.deltag)),mode = 'numeric')

names(n.term.dna.sum.deltag) = levels(n.term.cut.deltag)

for (i in 1:length(n.term.cut.deltag)){
  
  if (is.na(n.terminal.raw[,'Count.DNA'][i])){
    current = 0
  } else {
    
    current = n.terminal.raw[,'Count.DNA'][i]
  }
  
  n.term.dna.sum.deltag[n.term.cut.deltag[i]] = n.term.dna.sum.deltag[n.term.cut.deltag[i]] + current
  
  
}

png(paste0(results.dir,'DNA_coverage_per_delta_G_level_unflitered.png'),units="in", width=15, height=8.5, res=100)

barplot(n.term.dna.sum.deltag/deltag.hist$counts, ylab = 'DNA coverage normalized by number of variants at delta G content level', xlab = 'delta G content',
        main = 'DNA coverage per level of delta G content normalized to ammount of varaints at delta G content level')
dev.off()

cor(n.term.dna.sum.deltag[12:length(n.term.dna.sum.deltag)]/deltag.hist$counts[12:length(n.term.dna.sum.deltag)]
    ,c(12:length(n.term.dna.sum.deltag)), method ='spearman',use = 'pairwise.complete.obs')


par(mfrow = c(3,1))
hist(log10(n.terminal.tai.high.prom[,'Count.DNA']),ylim = c(0,5000))
hist(log10(n.terminal.tai.low.prom[,'Count.DNA']),ylim = c(0,5000))
hist(log10(n.terminal.raw[,'Count.DNA']),ylim = c(0,5000))




png(paste0(results.dir,'DNA_coverage_vs_deta_G_by_promoter_unflitered.png'),units="in", width=15, height=15, res=100)
par(mfrow = c(2,1))
plot(n.terminal.tai.high.prom[,'dG'],log10(n.terminal.tai.high.prom[,'Count.DNA'])
     ,main = 'High Promoter',ylab = 'Log DNA Coverage', xlab = 'Delta G of variable region'
     ,xlim = c(min(n.terminal.raw[,'dG']),max(n.terminal.raw[,'dG']))
     ,ylim = c(log10(min(n.terminal.raw[,'Count.DNA'],na.rm = TRUE)),log10(max(n.terminal.raw[,'Count.DNA'],na.rm = TRUE)))
     ,col="#00000080")
abline(lm(log10(n.terminal.tai.high.prom[,'Count.DNA'])~n.terminal.tai.high.prom[,'dG']), col="red") # regression line (y~x) 
plot(n.terminal.tai.low.prom[,'dG'],log10(n.terminal.tai.low.prom[,'Count.DNA'])
     ,main = 'Low Promoter',ylab = 'Log DNA Coverage', xlab = 'Delta G of variable region'
     ,xlim = c(min(n.terminal.raw[,'dG']),max(n.terminal.raw[,'dG']))
     ,ylim = c(log10(min(n.terminal.raw[,'Count.DNA'],na.rm = TRUE)),log10(max(n.terminal.raw[,'Count.DNA'],na.rm = TRUE)))
     ,col="#00000080")
title("Log DNA coverage vs delta G content for different promoters", outer = TRUE, line = -1)
abline(lm(log10(n.terminal.tai.low.prom[,'Count.DNA'])~n.terminal.tai.low.prom[,'dG']), col="red") # regression line (y~x) 
dev.off()


png(paste0(results.dir,'DNA_coverage_vs_Delta_G_by_RBS_unflitered.png'),units="in", width=20, height=15, res=100)
par(mfrow = c(2,2))
plot(n.terminal.tai.strong.rbs[,'dG'],log10(n.terminal.tai.strong.rbs[,'Count.DNA'])
     ,main = 'Strong RBS',ylab = 'Log DNA Coverage', xlab = 'Delta G of variable region'
     ,xlim = c(min(n.terminal.raw[,'dG']),max(n.terminal.raw[,'dG']))
     ,ylim = c(log10(min(n.terminal.raw[,'Count.DNA'],na.rm = TRUE)),log10(max(n.terminal.raw[,'Count.DNA'],na.rm = TRUE)))
     ,col="#00000080")
abline(lm(log10(n.terminal.tai.strong.rbs[,'Count.DNA'])~n.terminal.tai.strong.rbs[,'dG']), col="red") # regression line (y~x) 

plot(n.terminal.tai.mid.rbs[,'dG'],log10(n.terminal.tai.mid.rbs[,'Count.DNA'])
     ,main = 'Mid RBS',ylab = 'Log DNA Coverage', xlab = 'Delta G of variable region'
     ,xlim = c(min(n.terminal.raw[,'dG']),max(n.terminal.raw[,'dG']))
     ,ylim = c(log10(min(n.terminal.raw[,'Count.DNA'],na.rm = TRUE)),log10(max(n.terminal.raw[,'Count.DNA'],na.rm = TRUE)))
     ,col="#00000080")
abline(lm(log10(n.terminal.tai.mid.rbs[,'Count.DNA'])~n.terminal.tai.mid.rbs[,'dG']), col="red") # regression line (y~x) 

plot(n.terminal.tai.weak.rbs[,'dG'],log10(n.terminal.tai.weak.rbs[,'Count.DNA'])
     ,main = 'Weak RBS',ylab = 'Log DNA Coverage', xlab = 'Delta G of variable region'
     ,xlim = c(min(n.terminal.raw[,'dG']),max(n.terminal.raw[,'dG']))
     ,ylim = c(log10(min(n.terminal.raw[,'Count.DNA'],na.rm = TRUE)),log10(max(n.terminal.raw[,'Count.DNA'],na.rm = TRUE)))
     ,col="#00000080")
abline(lm(log10(n.terminal.tai.weak.rbs[,'Count.DNA'])~n.terminal.tai.weak.rbs[,'dG']), col="red") # regression line (y~x) 

plot(n.terminal.tai.wt.rbs[,'dG'],log10(n.terminal.tai.wt.rbs[,'Count.DNA'])
     ,main = 'WT RBS',ylab = 'Log DNA Coverage', xlab = 'Delta G of variable region'
     ,xlim = c(min(n.terminal.raw[,'dG']),max(n.terminal.raw[,'dG']))
     ,ylim = c(log10(min(n.terminal.raw[,'Count.DNA'],na.rm = TRUE)),log10(max(n.terminal.raw[,'Count.DNA'],na.rm = TRUE)))
     ,col="#00000080")
abline(lm(log10(n.terminal.tai.wt.rbs[,'Count.DNA'])~n.terminal.tai.wt.rbs[,'dG']), col="red") # regression line (y~x) 

title("Log DNA coverage vs Delta G for different promoters", outer = TRUE, line = -1)

dev.off()





png(paste0(results.dir,'DNA_coverage_vs_delta_G_by_RBS_and_promoter_unflitered.png'),units="in", width=20, height=15, res=100)
par(mfrow = c(4,2))
plot(n.terminal.tai.high.prom[which(n.terminal.raw[,'RBS.Display'] == 'Strong'),'dG'],n.terminal.tai.high.prom[which(n.terminal.raw[,'RBS.Display'] == 'Strong'),'Count.DNA']
     ,main = 'Strong RBS and High Promoter',ylab = 'Log DNA Coverage', xlab = 'Average delta G of variable region'
     ,xlim = c(min(n.terminal.raw[,'dG']),max(n.terminal.raw[,'dG']))
     ,ylim = c(min(n.terminal.raw[,'Count.DNA'],na.rm = TRUE),max(n.terminal.raw[,'Count.DNA'],na.rm = TRUE))
     ,col="#00000080", log = 'y')
abline(lm(log10(n.terminal.tai.high.prom[which(n.terminal.raw[,'RBS.Display'] == 'Strong'),'Count.DNA'])~
            n.terminal.tai.high.prom[which(n.terminal.raw[,'RBS.Display'] == 'Strong'),'dG'])
                ,col="red") # regression line (y~x) 
text(-55,1000,paste('R^2 = ',
                   round(summary(
                     lm(log10(n.terminal.tai.high.prom[which(n.terminal.raw[,'RBS.Display'] == 'Strong'),'Count.DNA'])~
                                               n.terminal.tai.high.prom[which(n.terminal.raw[,'RBS.Display'] == 'Strong'),'dG']))$adj.r.squared
                         ,2), '\nCorrelation = ', 
                   round(cor(log10(n.terminal.tai.high.prom[which(n.terminal.raw[,'RBS.Display'] == 'Strong'),'Count.DNA']),
                       n.terminal.tai.high.prom[which(n.terminal.raw[,'RBS.Display'] == 'Strong'),'dG']
                       ,use = 'na.or.complete',method = 'spearman'),2))
     , col = 'red')

plot(n.terminal.tai.high.prom[which(n.terminal.raw[,'RBS.Display'] == 'Mid'),'dG'],
     n.terminal.tai.high.prom[which(n.terminal.raw[,'RBS.Display'] == 'Mid'),'Count.DNA']
     ,main = 'Mid RBS and High Promoter',ylab = 'Log DNA Coverage', xlab = 'Average delta G of variable region'
     ,xlim = c(min(n.terminal.raw[,'dG']),max(n.terminal.raw[,'dG']))
     ,ylim = c(min(n.terminal.raw[,'Count.DNA'],na.rm = TRUE),max(n.terminal.raw[,'Count.DNA'],na.rm = TRUE))
     ,col="#00000080", log = 'y')
abline(lm(log10(n.terminal.tai.high.prom[which(n.terminal.raw[,'RBS.Display'] == 'Mid'),'Count.DNA'])~
            n.terminal.tai.high.prom[which(n.terminal.raw[,'RBS.Display'] == 'Mid'),'dG'])
       ,col="red") # regression line (y~x) 

text(-55,1000,paste('R^2 = ',
                   round(summary(
                     lm(log10(n.terminal.tai.high.prom[which(n.terminal.raw[,'RBS.Display'] == 'Mid'),'Count.DNA'])~
                          n.terminal.tai.high.prom[which(n.terminal.raw[,'RBS.Display'] == 'Mid'),'dG']))$adj.r.squared
                     ,2), '\nCorrelation = ', 
                   round(cor(log10(n.terminal.tai.high.prom[which(n.terminal.raw[,'RBS.Display'] == 'Mid'),'Count.DNA']),
                             n.terminal.tai.high.prom[which(n.terminal.raw[,'RBS.Display'] == 'Mid'),'dG']
                             ,use = 'na.or.complete',method = 'spearman'),2))
     , col = 'red')


plot(n.terminal.tai.high.prom[which(n.terminal.raw[,'RBS.Display'] == 'Weak'),'dG']
     ,n.terminal.tai.high.prom[which(n.terminal.raw[,'RBS.Display'] == 'Weak'),'Count.DNA']
     ,main = 'Weak RBS and High Promoter',ylab = 'Log DNA Coverage', xlab = 'Average delta G of variable region'
     ,xlim = c(min(n.terminal.raw[,'dG']),max(n.terminal.raw[,'dG']))
     ,ylim = c(min(n.terminal.raw[,'Count.DNA'],na.rm = TRUE),max(n.terminal.raw[,'Count.DNA'],na.rm = TRUE))
     ,col="#00000080", log = 'y')
abline(lm(log10(n.terminal.tai.high.prom[which(n.terminal.raw[,'RBS.Display'] == 'Weak'),'Count.DNA'])~
            n.terminal.tai.high.prom[which(n.terminal.raw[,'RBS.Display'] == 'Weak'),'dG'])
       ,col="red") # regression line (y~x) 

text(-55,1000,paste('R^2 = ',
                   round(summary(
                     lm(log10(n.terminal.tai.high.prom[which(n.terminal.raw[,'RBS.Display'] == 'Weak'),'Count.DNA'])~
                          n.terminal.tai.high.prom[which(n.terminal.raw[,'RBS.Display'] == 'Weak'),'dG']))$adj.r.squared
                     ,2), '\nCorrelation = ', 
                   round(cor(log10(n.terminal.tai.high.prom[which(n.terminal.raw[,'RBS.Display'] == 'Weak'),'Count.DNA']),
                             n.terminal.tai.high.prom[which(n.terminal.raw[,'RBS.Display'] == 'Weak'),'dG']
                             ,use = 'na.or.complete',method = 'spearman'),2))
     , col = 'red')



plot(n.terminal.tai.high.prom[which(n.terminal.raw[,'RBS.Display'] == 'WT'),'dG']
     ,n.terminal.tai.high.prom[which(n.terminal.raw[,'RBS.Display'] == 'WT'),'Count.DNA']
     ,main = 'WT RBS and High Promoter',ylab = 'Log DNA Coverage', xlab = 'Average delta G of variable region'
     ,xlim = c(min(n.terminal.raw[,'dG']),max(n.terminal.raw[,'dG']))
     ,ylim = c(min(n.terminal.raw[,'Count.DNA'],na.rm = TRUE),max(n.terminal.raw[,'Count.DNA'],na.rm = TRUE))
     ,col="#00000080",log = 'y')
abline(lm(log10(n.terminal.tai.high.prom[which(n.terminal.raw[,'RBS.Display'] == 'WT'),'Count.DNA'])~
            n.terminal.tai.high.prom[which(n.terminal.raw[,'RBS.Display'] == 'WT'),'dG'])
       ,col="red") # regression line (y~x) 


text(-55,1000,paste('R^2 = ',
                   round(summary(
                     lm(log10(n.terminal.tai.high.prom[which(n.terminal.raw[,'RBS.Display'] == 'WT'),'Count.DNA'])~
                          n.terminal.tai.high.prom[which(n.terminal.raw[,'RBS.Display'] == 'WT'),'dG']))$adj.r.squared
                     ,2), '\nCorrelation = ', 
                   round(cor(log10(n.terminal.tai.high.prom[which(n.terminal.raw[,'RBS.Display'] == 'WT'),'Count.DNA']),
                             n.terminal.tai.high.prom[which(n.terminal.raw[,'RBS.Display'] == 'WT'),'dG']
                             ,use = 'na.or.complete',method = 'spearman'),2))
     , col = 'red')


plot(n.terminal.tai.low.prom[which(n.terminal.raw[,'RBS.Display'] == 'Strong'),'dG'],
     n.terminal.tai.low.prom[which(n.terminal.raw[,'RBS.Display'] == 'Strong'),'Count.DNA']
     ,main = 'Strong RBS and Low Promoter',ylab = 'Log DNA Coverage', xlab = 'Average delta G of variable region'
     ,xlim = c(min(n.terminal.raw[,'dG']),max(n.terminal.raw[,'dG']))
     ,ylim = c(min(n.terminal.raw[,'Count.DNA'],na.rm = TRUE),max(n.terminal.raw[,'Count.DNA'],na.rm = TRUE))
     ,col="#00000080", log= 'y')
abline(lm(log10(n.terminal.tai.low.prom[which(n.terminal.raw[,'RBS.Display'] == 'Strong'),'Count.DNA'])~
            n.terminal.tai.low.prom[which(n.terminal.raw[,'RBS.Display'] == 'Strong'),'dG'])
       ,col="red") # regression line (y~x) 

text(-55,1000,paste('R^2 = ',
                   round(summary(
                     lm(log10(n.terminal.tai.low.prom[which(n.terminal.raw[,'RBS.Display'] == 'Strong'),'Count.DNA'])~
                          n.terminal.tai.low.prom[which(n.terminal.raw[,'RBS.Display'] == 'Strong'),'dG']))$adj.r.squared
                     ,2), '\nCorrelation = ', 
                   round(cor(log10(n.terminal.tai.low.prom[which(n.terminal.raw[,'RBS.Display'] == 'Strong'),'Count.DNA']),
                             n.terminal.tai.low.prom[which(n.terminal.raw[,'RBS.Display'] == 'Strong'),'dG']
                             ,use = 'na.or.complete',method = 'spearman'),2))
     , col = 'red')

plot(n.terminal.tai.low.prom[which(n.terminal.raw[,'RBS.Display'] == 'Mid'),'dG'],
     n.terminal.tai.low.prom[which(n.terminal.raw[,'RBS.Display'] == 'Mid'),'Count.DNA']
     ,main = 'Mid RBS and Low Promoter',ylab = 'Log DNA Coverage', xlab = 'Average delta G of variable region'
     ,xlim = c(min(n.terminal.raw[,'dG']),max(n.terminal.raw[,'dG']))
     ,ylim = c(min(n.terminal.raw[,'Count.DNA'],na.rm = TRUE),max(n.terminal.raw[,'Count.DNA'],na.rm = TRUE))
     ,col="#00000080", log = 'y')
abline(lm(log10(n.terminal.tai.low.prom[which(n.terminal.raw[,'RBS.Display'] == 'Mid'),'Count.DNA'])~
            n.terminal.tai.low.prom[which(n.terminal.raw[,'RBS.Display'] == 'Mid'),'dG'])
       ,col="red") # regression line (y~x) 
text(-55,1000,paste('R^2 = ',
                   round(summary(
                     lm(log10(n.terminal.tai.low.prom[which(n.terminal.raw[,'RBS.Display'] == 'Mid'),'Count.DNA'])~
                          n.terminal.tai.low.prom[which(n.terminal.raw[,'RBS.Display'] == 'Mid'),'dG']))$adj.r.squared
                     ,2), '\nCorrelation = ', 
                   round(cor(log10(n.terminal.tai.low.prom[which(n.terminal.raw[,'RBS.Display'] == 'Mid'),'Count.DNA']),
                             n.terminal.tai.low.prom[which(n.terminal.raw[,'RBS.Display'] == 'Mid'),'dG']
                             ,use = 'na.or.complete',method = 'spearman'),2))
     , col = 'red')


plot(n.terminal.tai.low.prom[which(n.terminal.raw[,'RBS.Display'] == 'Weak'),'dG'],
     n.terminal.tai.low.prom[which(n.terminal.raw[,'RBS.Display'] == 'Weak'),'Count.DNA']
     ,main = 'Weak RBS and Low Promoter',ylab = 'Log DNA Coverage', xlab = 'Average delta G of variable region'
     ,xlim = c(min(n.terminal.raw[,'dG']),max(n.terminal.raw[,'dG']))
     ,ylim = c(min(n.terminal.raw[,'Count.DNA'],na.rm = TRUE),max(n.terminal.raw[,'Count.DNA'],na.rm = TRUE))
     ,col="#00000080", log = 'y')
abline(lm(log10(n.terminal.tai.low.prom[which(n.terminal.raw[,'RBS.Display'] == 'Weak'),'Count.DNA'])~
            n.terminal.tai.low.prom[which(n.terminal.raw[,'RBS.Display'] == 'Weak'),'dG'])
       ,col="red") # regression line (y~x) 

text(-55,1000,paste('R^2 = ',
                   round(summary(
                     lm(log10(n.terminal.tai.low.prom[which(n.terminal.raw[,'RBS.Display'] == 'Weak'),'Count.DNA'])~
                          n.terminal.tai.low.prom[which(n.terminal.raw[,'RBS.Display'] == 'Weak'),'dG']))$adj.r.squared
                     ,2), '\nCorrelation = ', 
                   round(cor(log10(n.terminal.tai.low.prom[which(n.terminal.raw[,'RBS.Display'] == 'Weak'),'Count.DNA']),
                             n.terminal.tai.low.prom[which(n.terminal.raw[,'RBS.Display'] == 'Weak'),'dG']
                             ,use = 'na.or.complete',method = 'spearman'),2))
     , col = 'red')


plot(n.terminal.tai.low.prom[which(n.terminal.raw[,'RBS.Display'] == 'WT'),'dG'],
     n.terminal.tai.low.prom[which(n.terminal.raw[,'RBS.Display'] == 'WT'),'Count.DNA']
     ,main = 'WT RBS and Low Promoter',ylab = 'Log DNA Coverage', xlab = 'Average delta G of variable region'
     ,xlim = c(min(n.terminal.raw[,'dG']),max(n.terminal.raw[,'dG']))
     ,ylim = c(min(n.terminal.raw[,'Count.DNA'],na.rm = TRUE),max(n.terminal.raw[,'Count.DNA'],na.rm = TRUE))
     ,col="#00000080", log = 'y')
abline(lm(log10(n.terminal.tai.low.prom[which(n.terminal.raw[,'RBS.Display'] == 'WT'),'Count.DNA'])~
            n.terminal.tai.low.prom[which(n.terminal.raw[,'RBS.Display'] == 'WT'),'dG'])
       ,col="red") # regression line (y~x) 
text(-55,1000,paste('R^2 = ',
                   round(summary(
                     lm(log10(n.terminal.tai.low.prom[which(n.terminal.raw[,'RBS.Display'] == 'WT'),'Count.DNA'])~
                          n.terminal.tai.low.prom[which(n.terminal.raw[,'RBS.Display'] == 'WT'),'dG']))$adj.r.squared
                     ,2), '\nCorrelation = ', 
                   round(cor(log10(n.terminal.tai.low.prom[which(n.terminal.raw[,'RBS.Display'] == 'WT'),'Count.DNA']),
                             n.terminal.tai.low.prom[which(n.terminal.raw[,'RBS.Display'] == 'WT'),'dG']
                             ,use = 'na.or.complete',method = 'spearman'),2))
     , col = 'red')
title("Log DNA coverage vs delta G for different RBS promoter pairs", outer = TRUE, line = -1)

dev.off()


png(paste0(results.dir,'DNA_coverage_vs_delta_G_unflitered.png'),units="in", width=11, height=8.5, res=100)
plot(n.terminal.raw[,'dG'],n.terminal.raw[,'Count.DNA'],log = 'y'
     ,main = 'Log DNA coverage vs dG content',ylab = 'Log DNA Coverage',
     xlab = 'Average dG content of variable region',col="#00000080")
dev.off()



#C 2 R plots
#ADD ABline diagonal
c.2.r.matrix = read.csv('C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\c_2_r_summary_prot_dna.csv')


c.2.r.matrix.strong   <- c.2.r.matrix[which(c.2.r.matrix[,'RBS'] == 'Strong'),]
c.2.r.matrix.mid   <- c.2.r.matrix[which(c.2.r.matrix[,'RBS'] == 'Mid'),]
c.2.r.matrix.weak   <- c.2.r.matrix[which(c.2.r.matrix[,'RBS'] == 'Weak'),]
c.2.r.matrix.wt   <- c.2.r.matrix[which(c.2.r.matrix[,'RBS'] == 'WT'),]

c.2.r.matrix.high   <- c.2.r.matrix[which(c.2.r.matrix[,'Promoter'] == 'High'),]
c.2.r.matrix.low  <- c.2.r.matrix[which(c.2.r.matrix[,'Promoter'] == 'Low'),]



png(paste0(results.dir,"common_vs_rare_DNA_coverage_for_all_RBS.png"),units="in", width=10, height=10, res=100)

par(mfcol = c(2,2),oma = c(2,2,2,2))
plot(c.2.r.matrix.strong[,'Rare_DNA'],c.2.r.matrix.strong[,'Common_DNA'],
     main = 'Strong RBS', ylab = 'Common variant DNA coverage',
     xlab = 'Rare variant DNA coverage', ylim = c(0,8000), xlim = c(0,8000))
abline(0,1)
plot(c.2.r.matrix.weak[,'Rare_DNA'],c.2.r.matrix.weak[,'Common_DNA'],
     main = 'Weak RBS', ylab = 'Common variantDNA coverage',
     xlab = 'Rare variant DNA coverage', ylim = c(0,8000), xlim = c(0,8000))
abline(0,1)
plot(c.2.r.matrix.mid[,'Rare_DNA'],c.2.r.matrix.mid[,'Common_DNA'],
     main = 'Mid RBS', ylab = 'Common variant DNA coverage',
     xlab = 'Rare variant DNA coverage', ylim = c(0,8000), xlim = c(0,8000))
abline(0,1)
plot(c.2.r.matrix.wt[,'Rare_DNA'],c.2.r.matrix.wt[,'Common_DNA'],
     main = 'WT RBS', ylab = 'Common variant DNA coverage',
     xlab = 'Rare variant DNA coverage', ylim = c(0,8000), xlim = c(0,8000))
abline(0,1)
mtext("DNA coverage of rare vs common codon varaints of all RBSs", side = 3, outer = TRUE)

dev.off()


png(paste0(results.dir,"common_vs_rare_DNA_coverage_for_all_Promoters.png"),units="in", width=10, height=10, res=100)

par(mfcol = c(2,1),oma = c(2,2,2,2))
plot(c.2.r.matrix.high[,'Rare_DNA'],c.2.r.matrix.high[,'Common_DNA'],
     main = 'High Promoter', ylab = 'Common variant DNA coverage',
     xlab = 'Rare variant DNA coverage', ylim = c(0,8000), xlim = c(0,8000))
abline(0,1)
plot(c.2.r.matrix.low[,'Rare_DNA'],c.2.r.matrix.low[,'Common_DNA'],
     main = 'Low Promoter', ylab = 'Common variant DNA coverage',
     xlab = 'Rare variant DNA coverage', ylim = c(0,8000), xlim = c(0,8000))
abline(0,1)
mtext("DNA coverage of rare vs common codon varaints of all promoters", side = 3, outer = TRUE)

dev.off()


png(paste0(results.dir,"common_vs_rare_DNA_coverage_for_all_RBS_and_promoters_w.png"),units="in", width=12, height=6, res=100)

par(mfcol = c(2,4),oma = c(2,2,2,2))
plot(c.2.r.matrix.high[which(c.2.r.matrix.high[,'RBS'] == 'Strong'),'Rare_DNA'],
     c.2.r.matrix.high[which(c.2.r.matrix.high[,'RBS'] == 'Strong'),'Common_DNA'],
     main = 'High Promoter Strong RBS', ylab = 'Common variant DNA coverage',
     xlab = 'Rare variant DNA coverage', ylim = c(0,8000), xlim = c(0,8000))
abline(0,1)
plot(c.2.r.matrix.low[which(c.2.r.matrix.low[,'RBS'] == 'Strong'),'Rare_DNA'],
     c.2.r.matrix.low[which(c.2.r.matrix.low[,'RBS'] == 'Strong'),'Common_DNA'],
     main = 'Low Promoter Strong RBS', ylab = 'Common variant DNA coverage',
     xlab = 'Rare variant DNA coverage', ylim = c(0,8000), xlim = c(0,8000))
abline(0,1)
plot(c.2.r.matrix.high[which(c.2.r.matrix.high[,'RBS'] == 'Mid'),'Rare_DNA'],
     c.2.r.matrix.high[which(c.2.r.matrix.high[,'RBS'] == 'Mid'),'Common_DNA'],
     main = 'High Promoter Mid RBS', ylab = 'Common variant DNA coverage',
     xlab = 'Rare variant DNA coverage', ylim = c(0,8000), xlim = c(0,8000))
abline(0,1)
plot(c.2.r.matrix.low[which(c.2.r.matrix.low[,'RBS'] == 'Mid'),'Rare_DNA'],
     c.2.r.matrix.low[which(c.2.r.matrix.low[,'RBS'] == 'Mid'),'Common_DNA'],
     main = 'Low Promoter Mid RBS', ylab = 'Common variant DNA coverage',
     xlab = 'Rare variant DNA coverage', ylim = c(0,8000), xlim = c(0,8000))
abline(0,1)
plot(c.2.r.matrix.high[which(c.2.r.matrix.high[,'RBS'] == 'Weak'),'Rare_DNA'],
     c.2.r.matrix.high[which(c.2.r.matrix.high[,'RBS'] == 'Weak'),'Common_DNA'],
     main = 'High Promoter Weak RBS', ylab = 'Common variant DNA coverage',
     xlab = 'Rare variant DNA coverage', ylim = c(0,8000), xlim = c(0,8000))
abline(0,1)
plot(c.2.r.matrix.low[which(c.2.r.matrix.low[,'RBS'] == 'Weak'),'Rare_DNA'],
     c.2.r.matrix.low[which(c.2.r.matrix.low[,'RBS'] == 'Weak'),'Common_DNA'],
     main = 'Low Promoter Weak RBS', ylab = 'Common variant DNA coverage',
     xlab = 'Rare variant DNA coverage', ylim = c(0,8000), xlim = c(0,8000))
abline(0,1)
plot(c.2.r.matrix.high[which(c.2.r.matrix.high[,'RBS'] == 'WT'),'Rare_DNA'],
     c.2.r.matrix.high[which(c.2.r.matrix.high[,'RBS'] == 'WT'),'Common_DNA'],
     main = 'High Promoter WT RBS', ylab = 'Common variant DNA coverage',
     xlab = 'Rare variant DNA coverage', ylim = c(0,8000), xlim = c(0,8000))
abline(0,1)
plot(c.2.r.matrix.low[which(c.2.r.matrix.low[,'RBS'] == 'WT'),'Rare_DNA'],
     c.2.r.matrix.low[which(c.2.r.matrix.low[,'RBS'] == 'WT'),'Common_DNA'],
     main = 'Low Promoter WT RBS', ylab = 'Common variant DNA coverage',
     xlab = 'Rare variant DNA coverage', ylim = c(0,8000), xlim = c(0,8000))
abline(0,1)
mtext("DNA coverage of rare vs common codon varaints of all RBSs and promoters", side = 3, outer = TRUE)

dev.off()

cds.type.boxplot = function(data.set){

  cds.types = sort(unique(as.character(n.terminal.raw[,'CDS.type'])))
  dna.cds.types.split = vector("list", length(cds.types))
  
  names(dna.cds.types.split) = gsub('גˆ†','delta',cds.types)
  for (i in 1:length(cds.types)){
    dna.cds.types.split[[i]] = data.set[which(data.set[,'CDS.type'] == cds.types[i]),'Count.DNA']
  }
  dna.cds.types.split = c(dna.cds.types.split[-5],dna.cds.types.split[5])
  names(dna.cds.types.split)[0:2] = c('Rare', 'Common')
  
  return(dna.cds.types.split)  
  
  
}

boxplot(cds.type.boxplot(n.terminal.tai.weak.rbs), log = 'y',main = 'DNA coverage by CDS type', ylab = 'Log DNA coverage',xlab = 'CDS type')

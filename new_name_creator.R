require(gdata)
n.term = 
  read.csv('C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\n_terminal_data_with_full_seq.csv')

new.name = vector()

for (i in 1:nrow(n.term)){
  promoter.letter  = substr(n.term[i,'Promoter.Display'],1,1)
  RBS.letter  = toupper(substr(n.term[i,'RBS.Display'],1,2))
  gene.number = as.numeric(n.term[i,'Gene'])
  if (startsWith(n.term[i,'CDS.type'],'Max')){
    type = 'R'
  } else if (startsWith(n.term[i,'CDS.type'],'Min')){
    type = 'C'
  
} else if (startsWith(n.term[i,'CDS.type'],'WT')){
  type = 'W'
} else if (startsWith(n.term[i,'CDS.type'],'גˆ†G')){
  type = substr(unlist(strsplit(as.character(n.term[i,'CDS.type']),split = ' '))[2],1,2)
}
new.name <- append(new.name,paste(promoter.letter,RBS.letter,gene.number,type,sep = '_'))

}

write.csv(cbind(as.character(n.term$Name), new.name),'C:\\Users\\dell7\\Documents\\Tzachi\\workspace\\data\\new_name.csv',quote = F)

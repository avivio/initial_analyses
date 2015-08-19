


scatter  <- ggplot(cor_set,aes(x = TASEP.bottle.neck.depth, y = dG.noutr, color = n.pos.neg)) +
  geom_point(alpha = 0.5)  + 
  ylab('Delta G no UTR') + 
  xlab('Bottleneck depth calculated by TASEP') + 
  ggtitle('Delta G no UTR vs Bottleneck depth\n day 12 no wall high promoter above 14') +  
  guides(color=guide_legend(title='5 lineages\nfitness residual')) + 
  theme_aviv +  
  theme(legend.position = c(1, 1), legend.justification = c(1, 1)) 


p.top <- p  <- ggplot(cor_set,aes(x = TASEP.bottle.neck.depth, fill = n.pos.neg)) +
  geom_density(alpha = 0.5)  + 
  xlab('Bottleneck depth calculated by TASEP') + 
  theme(legend.position = "none") + 
  theme_aviv + 
  theme(legend.position = "none") 



p.right <- p  <- ggplot(cor_set,aes(x = dG.noutr, fill = n.pos.neg)) +
  geom_density(alpha = 0.5)  + 
  xlab('Delta G no UTR') + 
  coord_flip() + 
  theme_aviv + 
  theme(legend.position = "none") 
  

empty <- ggplot() + geom_point(aes(1, 1), colour = "white") + theme(plot.background = element_blank(), 
                                                                    panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                                                    panel.border = element_blank(), panel.background = element_blank(), axis.title.x = element_blank(), 
                                                                    axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), 
                                                                    axis.ticks = element_blank())
grid.arrange(p.top, empty, scatter, p.right, ncol = 2, nrow = 2, widths = c(4, 
                                                                                  1), heights = c(1, 4))


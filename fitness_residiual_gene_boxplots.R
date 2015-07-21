sorted.genes  <- fitseq.data.residuals.above.14%>%
  filter(day == 12) %>% 
  group_by(Gene) %>% 
  summarize(med.fit = median(fit.resid)) %>%
  arrange(med.fit)

fitseq.data.residuals.above.14$Gene <- factor(fitseq.data.residuals.above.14$Gene,
                                       levels = sorted.genes$Gene)

ggplot(fitseq.data.residuals.above.14,
       aes(x = Gene, y = fit.resid)) +
  geom_boxplot() +
  theme_aviv +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ylab('Fitness residual') +
  ggtitle('Fitness residuals by gene')
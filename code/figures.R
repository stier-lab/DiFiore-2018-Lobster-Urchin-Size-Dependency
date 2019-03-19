###############################################################################
## Basic plot of number consumed by predator and prey size at max density
##############################################################################

p1 <- df %>%
  filter(density > 24) %>%
ggplot(aes(x = class, y = num_consumed))+
  geom_boxplot(outlier.shape = 4, fill = NA)+
  geom_jitter(aes(color = class), width = 0.3, show.legend = F)+
  scale_color_manual(values = c('#1b9e77','#d95f02','#7570b3','#e7298a'))+
  labs(y = "Number of prey eaten", x = "Predator size category")+
  theme_pubclean()


p2 <- df %>%
  filter(density > 24) %>%
  ggplot(aes(x = treatment, y = num_consumed))+
  geom_boxplot(outlier.shape = 4, fill = NA)+
  geom_jitter(aes(color = treatment), width = 0.3, show.legend = F)+
  scale_color_manual(values = c('#d95f02','#7570b3','#e7298a'))+
  labs(y = "Number of prey eaten", x = "Prey size category")+
  scale_x_discrete(labels = c('large', 'medium', 'small'))+
  theme_pubclean()

fig1 <- cowplot::plot_grid(p1, p2 + labs(y = NULL))


png("figures/fig1_190312.png", width = 623*1.2, height = 369*1.2)
fig1
dev.off()

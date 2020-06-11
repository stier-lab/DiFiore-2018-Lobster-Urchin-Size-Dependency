library(here)
source(here("code", "1_setup.R"))

#--------------------------------------------------------------------------------
## Figure on max consumption
#--------------------------------------------------------------------------------

df <- read.table(here("data/cleaned","loburc_cleaned.csv"), header = T, sep = ",") %>% 
  arrange(treatment, id) %>%
  group_by(id) %>%
  filter(initial == max(initial)) %>%
  mutate(mort = killed/initial) %>%
  filter(id != "N07") %>%
  ungroup() %>%
  mutate(treatment.clean = recode(treatment, urc_small = "Small\nurchins", urc_medium = "Medium\nurchins", urc_large = "Large\nurchins") )


res <- cbind(df$killed, df$initial-df$killed)

p1 <- ggplot(df, aes(x = treatment.clean, y = mort))+
  geom_boxplot(aes(fill = treatment.clean), outlier.shape = NA)+
  scale_fill_manual(values = colors)+
  geom_jitter(aes(shape = treatment.clean, fill = treatment.clean), color = "white", size = 2, width = 0.25)+
  scale_shape_manual(values = c(21,23,24))+
  labs(x = "", y = "Proportional mortality")+
  annotate("text", x = c(1,2,3.1), y = c(0.25, 0.35, 0.9), label = c("a", "b", "b"))+
  theme(legend.position = "none")

glmer1 <- glmer(res ~ 0 + treatment + (1|id), data = df, family = "binomial")
summary(glmer1)
summary(multcomp::glht(glmer1, linfct = multcomp::mcp(treatment = "Tukey")))
inverse.logit(fixef(glmer1))
hist(residuals(glmer1))


glmer2 <- glmer(res ~ scale(size) * treatment + (1|id), data = df, family = "binomial")
summary(glmer2)
summary(multcomp::glht(glmer2, linfct = multcomp::mcp(treatment = "Tukey")))

newdat <- expand.grid(treatment = unique(df$treatment), size = seq(min(df$size), max(df$size), length.out = 1000))

newdat$predicted <- predict(glmer2, newdata = newdat, re.form = NA, type = "response")

newdat <- newdat %>%
  mutate(treatment.clean = recode(treatment, urc_small = "Small\nurchins", urc_medium = "Medium\nurchins", urc_large = "Large\nurchins") )


p2 <- ggplot(df, aes(x = size, y = mort))+
  geom_line(data = newdat, aes(y = predicted, color = treatment.clean, linetype = treatment.clean), lwd = 1.1)+
  geom_point(aes(fill = treatment.clean, shape = treatment.clean), color = "white", size = 2)+
  scale_shape_manual(values = c(21,23,24))+
  scale_fill_manual(values = colors)+
  scale_color_manual(values = colors)+
  scale_linetype_manual(values = c(1,2,3))+
  labs(x = "\nCarpace length (mm)", y = NULL, color = NULL, fill = NULL, linetype = NULL, shape = NULL)+
  theme(legend.key.size = unit(1.5, "cm"))


fig2 <- cowplot::plot_grid(p1,p2, nrow = 1, rel_widths = c(1, 1.75))

ggsave(here("figures", "fig2-mortality.png"), fig2, device = "png", width = 8.5, height = 5 )


















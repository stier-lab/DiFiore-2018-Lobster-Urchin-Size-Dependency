library(here)
source(here("code", "1_setup.R"))
source(here("code", "Base_functions/functions.R"))


#####################################
## Get data
#####################################

df <- read.table(here("data/cleaned","loburc_cleaned.csv"), header = T, sep = ",") %>% 
  arrange(treatment, id) %>%
  group_by(id) %>%
  # filter(initial == max(initial)) %>%
  mutate(m.int = initial*mr, 
         m.kill = killed*mr,
         prop.m.kill = m.kill/m.int) %>%
  filter(id != "N07") %>%
  ungroup() %>%
  mutate(treatment.clean = recode(treatment, urc_small = "Small\nurchins", urc_medium = "Medium\nurchins", urc_large = "Large\nurchins") )


mf <- df %>% mutate(m.intr = round(m.int, 0), 
              m.killr = round(m.kill, 0))

res <- cbind(mf$m.killr, mf$m.intr-mf$m.killr)

mod <- glm(res ~ mc*treatment, mf, family = "binomial")
summary(mod)

mod2 <- glm(res ~ mc, mf, family = "binomial")
summary(mod2)

AIC(mod, mod2)
anova(mod2, mod)


pframe <- expand.grid(treatment = unique(df$treatment), mc = seq(min(mf$mc, na.rm =T), max(mf$mc, na.rm =T), length.out = 100))
pp <- predict(mod, newdata = pframe, se.fit = T)


linkinv <- family(mod)$linkinv ## inverse-link function
#Put together the prediction, lower and upper bounds, and back-transform to the probability scale:
  
pframe$pred0 <- pp$fit
pframe$prop.m.kill <- linkinv(pp$fit)
alpha <- 0.95
sc <- abs(qnorm((1-alpha)/2))  ## Normal approx. to likelihood
alpha2 <- 0.5
sc2 <- abs(qnorm((1-alpha2)/2))  ## Normal approx. to likelihood
pframe <- transform(pframe,
                    lwr=linkinv(pred0-sc*pp$se.fit),
                    upr=linkinv(pred0+sc*pp$se.fit),
                    lwr2=linkinv(pred0-sc2*pp$se.fit),
                    upr2=linkinv(pred0+sc2*pp$se.fit))


ps1 <- ggplot(df, aes(x = mc, y = prop.m.kill))+
  geom_point(aes(color = treatment), pch = 21)+
  geom_ribbon(data = pframe, aes(x = mc, ymin = lwr, ymax = upr, group = treatment), alpha = 0.25, color = "gray")+
  geom_line(data = pframe, aes(x = mc, y = prop.m.kill, color = treatment), lwd = 1.2)+
  scale_color_manual(values = c('#AF8DC3','#C3AF8D','#8DC3AF'))+
  labs(x = "Predator body size (g)", y = "Grams of prey eaten per gram offered", color = "")+
  theme(legend.position = c(0.6, 0.8))

ggsave(here::here("figures/", "biomass-maxconsumption.png"), ps1, width = 8, height = 6)


#---------------------------------------------------------------------------------
## Bayesian version using RStanArm
#---------------------------------------------------------------------------------

stan1 <- rstanarm::stan_glm(res ~ mc * treatment, data = mf, family = binomial(link = "logit"))
plot(stan1)
launch_shinystan(stan1)


p1 <- ggplot(mf, aes(x = treatment, y = prop.m.kill))+
  geom_boxplot(aes(fill = treatment), outlier.shape = NA)+
  scale_fill_manual(values = c('#AF8DC3','#C3AF8D','#8DC3AF') )+
  geom_jitter(aes(shape = treatment, fill = treatment), color = "white", size = 2, width = 0.25)+
  scale_shape_manual(values = c(21,23,24))+
  labs(x = "Urchin size", y = expression(paste("Proportion of biomass eaten (g g"^"-1", ")")), fill = "Urchin size", shape = "Urchin size")+
  scale_x_discrete(labels = c("Large", "Medium", "Small"))+
  theme(legend.position = "none")

p2 <- mf %>%
  group_by(treatment) %>%
  data_grid(mc = seq_range(mc, n = 51)) %>%
  add_fitted_draws(stan1) %>%
  ggplot(aes(x = mc, y = prop.m.kill, color = treatment, fill = treatment)) +
  stat_lineribbon(aes(y = .value), .width = 0.95, fill = "grey90") +
  geom_point(data = df, aes(shape = treatment)) +
  scale_shape_manual(values = c(21,23,24), labels = c("Large", "Medium", "Small"))+
  scale_color_manual(values = c('#AF8DC3','#C3AF8D','#8DC3AF'), labels = c("Large", "Medium", "Small"))+
  scale_fill_manual(values = c('#AF8DC3','#C3AF8D','#8DC3AF'), labels = c("Large", "Medium", "Small"))+
  labs(x = "Predator body size (g)", y = expression(paste("Proportion of biomass eaten (g g"^"-1", ")")), color = "Urchin size", fill = "Urchin size", shape = "Urchin size")+
  theme(legend.position = c(0.7, 0.8))

p2_bayes <- cowplot::plot_grid(p1, p2+labs(y = ""), align = "h", nrow = 1, rel_widths = c(0.5, 1))  

ggsave(filename = here::here("figures/", "p2_bayesian.png"), plot = p2_bayes, device = "png", width = 8.5, height = 8.5/2)










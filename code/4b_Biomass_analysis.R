library(here)
source(here("code", "1_setup.R"))
source(here("code", "Base_functions/functions.R"))
source(here("code", "theme.R"))


#####################################
## Get data
#####################################

df <- read.csv(here("data/cleaned","loburc_cleaned.csv"), header = T) %>% 
  arrange(treatment, id) %>%
  group_by(id) %>%
  filter(initial == max(initial)) %>%
  mutate(m.int = initial*mr, 
         m.kill = killed*mr,
         prop.m.kill = m.kill/m.int) %>%
  filter(id != "N07") %>%
  ungroup() %>%
  mutate(treatment.clean = recode(treatment, urc_small = "Small\nurchins", urc_medium = "Medium\nurchins", urc_large = "Large\nurchins") )

#-----------------------------------------
## Summary stats for P2 of results
#-----------------------------------------


# Across all densities only the largest lobsters regularly consumed the largest urchins, and 0.0X% of large urchins were consumed by lobster less than XXX g.

read.csv(here("data/cleaned","loburc_cleaned.csv"), header = T) %>% 
  arrange(treatment, id) %>%
  group_by(id) %>%
  mutate(m.int = initial*mr, 
         m.kill = killed*mr,
         prop.m.kill = m.kill/m.int) %>%
  filter(id != "N07") %>%
  ungroup() %>%
  mutate(treatment.clean = recode(treatment, urc_small = "Small\nurchins", urc_medium = "Medium\nurchins", urc_large = "Large\nurchins") ) %>%
  filter(treatment == "urc_large") %>%
  filter(mc <= median(mc)) %>% 
  summarize(initial = sum(initial), 
            killed = sum(killed))


#---------------------------
## Frequentist version
#---------------------------

mf <- df %>% mutate(m.intr = round(m.int, 0), 
              m.killr = round(m.kill, 0)) %>%
  group_by(id)

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

stan1 <- rstanarm::stan_glm(cbind(m.killr, m.intr-m.killr) ~ mc + treatment, data = mf, family = binomial(link = "logit"))
plot(stan1)
launch_shinystan(stan1)

get_variables(stan1)
prior_summary(stan1)

mod.predict <- ggeffects::ggpredict(stan1, terms = ~treatment)

mod.predict <- mod.predict %>% as.data.frame() %>%
  rename(treatment = x)

p1 <- ggplot(mf, aes(x = treatment, y = prop.m.kill))+
  scale_fill_manual(values = c('#AF8DC3','#C3AF8D','#8DC3AF') )+
  scale_color_manual(values = c('#AF8DC3','#C3AF8D','#8DC3AF') )+
  geom_jitter(aes(shape = treatment, fill = treatment, color = treatment), size = 2, alpha = 0.5)+
  scale_shape_manual(values = c(21,23,24))+
  geom_pointinterval(data = mod.predict, aes(x = treatment, y = predicted, ymin = conf.low, ymax = conf.high, shape = treatment, fill = treatment), size = 7.5)+
  labs(x = "Urchin size", y = expression(paste("Proportion of biomass eaten (g g"^"-1", ")")), fill = "Urchin size", shape = "Urchin size")+
  scale_x_discrete(labels = c("Large", "Medium", "Small"))+
  theme_bd()+
  theme(legend.position = "none")


mod.predict <- ggeffects::ggpredict(stan1, terms = ~mc + treatment)
plot(mod.predict)

mod.predict <- mod.predict %>% as.data.frame() %>%
  rename(mc = x,
         treatment = group)

p2 <- ggplot(mf, aes(x = mc, y = prop.m.kill, color = treatment))+
  geom_point(data = mf, aes(shape = treatment, fill = treatment), size = 2, alpha = 0.5) +
  geom_ribbon(data = mod.predict, aes(x = mc, y = predicted, ymin = conf.low, ymax = conf.high, group = treatment), fill = "black", alpha = 0.25, color = "transparent")+
  geom_line(data = mod.predict, aes(x = mc, y = predicted, color = treatment), lwd = 1.5)+
  scale_shape_manual(values = c(21,23,24), labels = c("Large", "Medium", "Small"))+
  scale_color_manual(values = c('#AF8DC3','#C3AF8D','#8DC3AF'), labels = c("Large", "Medium", "Small"))+
  scale_fill_manual(values = c('#AF8DC3','#C3AF8D','#8DC3AF'), labels = c("Large", "Medium", "Small"))+
  labs(x = "Predator body size (g)", y = expression(paste("Proportion of biomass eaten (g g"^"-1", ")")), color = "Urchin size", fill = "Urchin size", shape = "Urchin size")+
  theme_bd()


p2_bayes <- cowplot::plot_grid(p1, p2+labs(y = ""), align = "vh", nrow = 1, rel_widths = c(0.5, 1), axis = "bl") 

p2_bayes
ggsave(filename = here::here("figures/", "p2_bayesian.png"), plot = p2_bayes, device = "png", width = 10, height = 10/2)


#----------------------------------------------------------
## Alternative versions of the figure
#----------------------------------------------------------
temp <- mf %>%
  group_by(id, size, treatment) %>%
  filter(m.int == max(m.int))

mod1 <- lm(m.kill ~ mc * treatment, data = temp)
summary(mod1)
hist(residuals(mod1), breaks = 30)

mod2 <- lm(log(m.kill) ~ log(mc) * treatment, data = temp[temp$m.kill > 0, ])
summary(mod2)
hist(residuals(mod1), breaks = 30)

pred <- ggeffects::ggpredict(mod1, terms = ~mc*treatment) %>%
  as.data.frame() %>%
  rename(mc = x, 
         treatment = group)


p1 <- mf %>%
  group_by(id, size, treatment) %>%
  filter(initial == max(initial)) %>%
  ggplot(aes(x = treatment, y = m.kill, color = treatment))+
  geom_jitter(aes(shape = treatment, fill = treatment), size = 2, alpha = 0.5)+
  stat_summary(fun = "mean", geom = "point", size = 5, aes(shape = treatment, fill = treatment))+
  scale_shape_manual(values = c(21,23,24), labels = c("Large", "Medium", "Small"))+
  scale_color_manual(values = c('#AF8DC3','#C3AF8D','#8DC3AF'), labels = c("Large", "Medium", "Small"))+
  scale_fill_manual(values = c('#AF8DC3','#C3AF8D','#8DC3AF'), labels = c("Large", "Medium", "Small"))+
  labs(x = "Predator body size (g)", y = "Biomass of urchins consumed", color = "Urchin size", fill = "Urchin size", shape = "Urchin size")+
  theme_bd()

ggsave("figures/p1_forAS.png", plot = p1, width = 6, height = 4)

p2 <- mf %>%
  group_by(id, size, treatment) %>%
  filter(m.int == max(m.int)) %>%
  ggplot(aes(x = mc, y = m.kill, color = treatment))+
  geom_point(aes(shape = treatment, fill = treatment), size = 2) +
  geom_ribbon(data = pred, aes(x = mc, y = predicted, ymin = conf.low, ymax = conf.high, group = treatment), fill = "black", alpha = 0.1, color = "transparent")+
  geom_line(data = pred, aes(x = mc, y = predicted, color = treatment), lwd = 1.5)+
  scale_shape_manual(values = c(21,23,24), labels = c("Large", "Medium", "Small"))+
  scale_color_manual(values = c('#AF8DC3','#C3AF8D','#8DC3AF'), labels = c("Large", "Medium", "Small"))+
  scale_fill_manual(values = c('#AF8DC3','#C3AF8D','#8DC3AF'), labels = c("Large", "Medium", "Small"))+
  labs(x = "Predator body size (g)", y = "Biomass of urchins consumed", color = "Urchin size", fill = "Urchin size", shape = "Urchin size")+
  theme_bd()

ggsave("figures/p2_forAS.png", plot = p2, width = 6, height = 4)


pred2 <- ggeffects::ggpredict(mod2, terms = ~mc*treatment, back.transform = F) %>%
  as.data.frame() %>%
  rename(mc = x, 
         treatment = group)

p3 <- mf %>%
  group_by(id, size, treatment) %>%
  filter(m.int == max(m.int)) %>%
  filter(m.kill > 0) %>%
  ggplot(aes(x = log(mc), y = log(m.kill), color = treatment))+
  geom_point(aes(shape = treatment, fill = treatment), size = 2) +
  geom_ribbon(data = pred2, aes(x = log(mc), y = predicted, ymin = conf.low, ymax = conf.high, group = treatment), fill = "black", alpha = 0.1, color = "transparent")+
  geom_line(data = pred2, aes(x = log(mc), y = predicted, color = treatment), lwd = 1.5)+
  scale_shape_manual(values = c(21,23,24), labels = c("Large", "Medium", "Small"))+
  scale_color_manual(values = c('#AF8DC3','#C3AF8D','#8DC3AF'), labels = c("Large", "Medium", "Small"))+
  scale_fill_manual(values = c('#AF8DC3','#C3AF8D','#8DC3AF'), labels = c("Large", "Medium", "Small"))+
  labs(x = "Predator body size (g)", y = "Biomass of urchins consumed", color = "Urchin size", fill = "Urchin size", shape = "Urchin size")+
  theme_bd()

ggsave("figures/p3_forAS.png", plot = p3, width = 6, height = 4)

cowplot::plot_grid(p1+theme(legend.position = "none"), p2, rel_widths = c(0.5, 1))


# This analysis doesn't make sense in the context of our experiment because we did not control biomass, for example we would have needed to have a much higher density of small urchins in order to have the same total biomass as a high density of large urchins. Large urchin biomass in the high density trials ~ 1800 g, while the small urchin biomass in the high density trials was only ~80 g. Obviously, then we can't compare without controlling for the amount of biomass offered. In a different experiment (if we had fixed biomass and allowed density to vary), I would expect to see an effect of prey size on the biomass consumed. For example small predators might consume more small prey while large predators may consume more prey per unit mass of predator. In other words the slope and intercept of the relationship between predator size and biomass consumed would depende on prey size class. However, in our case we don't have a good estimate for how much biomass a large predator 


mf %>%
  group_by(id, size, treatment) %>%
  # filter(m.int == max(m.int)) %>%
  mutate(m.kill_div_mc = m.kill/mc) %>%
  ggplot(aes(x = mc, y = m.kill_div_mc, color = treatment))+
  geom_point(aes(shape = treatment, fill = treatment), size = 2, alpha = 0.5) +
  geom_smooth(method = "lm", aes(group = treatment))+
  scale_shape_manual(values = c(21,23,24), labels = c("Large", "Medium", "Small"))+
  scale_color_manual(values = c('#AF8DC3','#C3AF8D','#8DC3AF'), labels = c("Large", "Medium", "Small"))+
  scale_fill_manual(values = c('#AF8DC3','#C3AF8D','#8DC3AF'), labels = c("Large", "Medium", "Small"))+
  labs(x = "Predator body size (g)", y = "Biomass of urchins consumed / Biomass of lobster", color = "Urchin size", fill = "Urchin size", shape = "Urchin size")+
  theme_bd()




mf %>%
  group_by(id, size, treatment) %>%
  filter(initial == max(initial)) %>%
  mutate(m.kill_div_mc = m.kill/mc) %>%
  ggplot(aes(x = treatment, y = m.kill_div_mc, color = treatment))+
  geom_jitter(aes(shape = treatment, fill = treatment), size = 2, alpha = 0.5)+
  stat_summary(fun = "mean", geom = "point")































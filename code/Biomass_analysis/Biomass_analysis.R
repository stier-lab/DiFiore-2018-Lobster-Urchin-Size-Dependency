library(here)
source(here("code", "1_setup.R"))
source(here("code", "Base_functions/functions.R"))


#####################################
## Get data
#####################################

df <- read.table(here("data/cleaned","loburc_cleaned.csv"), header = T, sep = ",") %>% 
  arrange(treatment, id) %>%
  group_by(id) %>%
  filter(initial == max(initial)) %>%
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


ggplot(df, aes(x = mc, y = prop.m.kill))+
  geom_point(aes(color = treatment), pch = 21)+
  geom_ribbon(data = pframe, aes(x = mc, ymin = lwr, ymax = upr, group = treatment), alpha = 0.5, color = "gray")+
  geom_line(data = newdat, aes(x = mc, y = predicted, color = treatment), lwd = 1.2)+
  labs(x = "Lobster size (g)", y = "Grams of urchins eaten per gram offered", color = "")+
  theme(legend.position = c(0.6, 0.8))







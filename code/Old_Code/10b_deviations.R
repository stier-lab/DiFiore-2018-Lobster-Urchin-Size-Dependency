#-------------------------------------------
## Setup
#-------------------------------------------

library(here)
source(here("code", "Base_functions/Functions.R"))
source(here("code", "1_setup.R"))

#--------------------------------------------
## Get data
#--------------------------------------------

df <- read.csv(here::here("data/cleaned/posteriors", "observational_predictions.csv")) %>% as_tibble()

# Full error propogations
d <- df %>% group_by(year, site, estimate) %>%
  mutate(id.vec = 1:n()) %>%
  group_by(year, site, id.vec) %>%
  pivot_wider(names_from = estimate, values_from = prediction) %>%
  group_by(year, site) %>%
  fill(Rall, MTE.fixed) %>%
  # mutate(null = (null - mu)/(mu + null), 
  #        full = (full - mu)/(mu+ full), 
  #        MTE.fixed = (MTE.fixed - mu)/(mu + MTE.fixed), 
  #        Rall = (Rall - mu)/(mu + Rall)) %>%
  # 
  mutate_at(.vars = c("null", "full", "MTE.fixed", "Rall"), funs((. - mu))) %>%
  select(-mu, -id.vec) %>%
  pivot_longer(values_to = "prediction", names_to = "estimate", -c(year, site))

points <- d %>% group_by(year, site, estimate) %>% mutate(.draw = 1:n()) %>%
  sample_draws(500) %>%
  filter(estimate != "MTE")

# W/ parameter uncertaintly
d  %>%
  ungroup() %>%
  group_by(estimate) %>%
  mutate(.draw = 1:n()) %>%
  mean_qi(prediction, .width = c(0.75, 0.95)) %>%
  # summarize(mean = mean(prediction), 
  #           ci.lower = PI(prediction)[1], 
  #           ci.upper = PI(prediction)[2]) %>%
  filter(estimate != "MTE") %>%
  ggplot(aes(x = estimate, y = prediction))+
  geom_jitter(data = points, aes( x = estimate , y = prediction, color = estimate), pch = 21, alpha = 0.25)+
  scale_color_manual(values = c('#AF8DC3','#C3AF8D', '#c3958d', '#8DC3AF', '#c38db5'))+
  geom_pointinterval(alpha = 0.5)+
  geom_hline(yintercept = 0, lty = "dashed")+
  coord_cartesian(ylim = c(-0.002, 0.005) )
  # ylim(c(-10, 10))
  ylim(c(-0.002, 0.005))
  ylim(c(-0.002, 0.0085))


# Averaged across parameter uncertaintly

e <- df %>% group_by(year, site, estimate) %>% 
  summarize(prediction = mean(prediction)) %>%
  pivot_wider(names_from = estimate, values_from = prediction) %>%
  group_by(year, site) %>%
  mutate_at(.vars = c("null", "full", "MTE.fixed", "Rall"), funs((. - mu))) %>%
  select(-mu) %>%
  pivot_longer(values_to = "prediction", names_to = "estimate", -c(year, site)) %>%
  filter(estimate != "MTE")

e %>%
  group_by(estimate) %>%
  summarize(mean = mean(prediction), 
            ci.lower = PI(prediction)[1], 
            ci.upper = PI(prediction)[2]) %>%
  ggplot(aes(x = estimate, y = mean))+
  geom_point()+
  geom_jitter(data = e, aes(x= estimate, y = prediction))+
  geom_errorbar(aes(ymin = ci.lower, ymax = ci.upper))+
  geom_hline(yintercept = 0, lty = "dashed")+
  coord_cartesian(ylim = c(-0.001, 0.003))


# Combination



# 



points <- d %>% group_by(year, site, estimate) %>% mutate(.draw = 1:n()) %>%
  sample_draws(500) %>%
  filter(estimate != "MTE", estimate %in% c("full", "null")) %>%
  bind_rows(filter(e, estimate %in% c("MTE.fixed", "Rall"))) %>%
  ungroup() %>%
  mutate(estimate = forcats::fct_relevel(as.factor(estimate), "full", "null", "MTE.fixed", "Rall"), 
         estimate = recode(estimate, full = "Fully-integrated\nsize-structured", null = "Unstructured", MTE.fixed = "First principles\nexpectations", Rall = "Metanalysis\nexpectations"))

e %>%
  ungroup() %>%
  group_by(estimate) %>%
  filter(estimate %in% c("MTE.fixed", "Rall")) %>%
  mutate(.draw = 1:n()) %>%
  mean_qi(prediction, .width = c(0.75, 0.95)) -> half

panelc <- d  %>%
  ungroup() %>%
  group_by(estimate) %>%
  mutate(.draw = 1:n()) %>%
  mean_qi(prediction, .width = c(0.75, 0.95)) %>%
  filter(estimate %in% c("full", "null")) %>%
  bind_rows(half) %>%
  ungroup() %>%
  mutate(estimate = forcats::fct_relevel(as.factor(estimate), "full", "null", "MTE.fixed", "Rall"), 
         estimate = recode(estimate, full = "Fully-integrated\nsize-structured", null = "Unstructured", MTE.fixed = "First principles\nexpectations", Rall = "Metanalysis\nexpectations")) %>%
  ggplot(aes(x = estimate, y = prediction))+
  geom_jitter(data = points, aes( x = estimate , y = prediction, fill = estimate), pch = 21, alpha = 0.25, show.legend = F)+
  scale_color_manual(values = c('#AF8DC3','#C3AF8D', '#c3958d', '#8DC3AF', '#c38db5'))+
  geom_pointinterval()+
  geom_hline(yintercept = 0, lty = "dashed")+
  coord_cartesian(ylim = c(-0.002, 0.004) )+
  labs(x = "", y = expression(paste("Deviations in predicted consumption (ind. m"^-2,"h"^-1,")", sep = "")))#+
  # annotate("text", x = 4.5, y = 0.001, label = "Mean size-\nstructured model")


ggsave(here::here("figures/", "panelc_version10k.png"), panelc)



p6 <- cowplot::plot_grid(timeseries, panelc, align = "horizontal")
ggsave(here::here("figures/", "fig6panelbc_alt.png"), p6, width = 8.25*2, height = 8.25/3*2)










#--------------------------------------------------------
## scrap
#--------------------------------------------------------

ggplot(p, aes(x= prediction, y = estimate))+
  geom_pointintervalh(aes(color = as.factor(year)), position = position_dodge())+
  facet_wrap(~site)

ggplot(p, aes(x = prediction, y = estimate))+
  geom_point(aes(color = estimate))+
  geom_vline(xintercept = 0, linetype = "dashed")
  
ggplot(p2, aes(x = diff, y = estimate))+
  geom_pointintervalh()+
  geom_vline(xintercept = 0, linetype = "dashed")

ggplot(p2, aes(x = diff, y = estimate))+
  geom_point()+
  geom_vline(xintercept = 0, linetype = "dashed")



test <- data.frame(m1 = rnorm(100), m2 = rnorm(100, 1), m3 = rnorm(100, -1)) %>%
  gather(model, prediction) %>%
  group_by(model) %>%
  mean_qi(prediction) %>%
  mutate(diff = prediction - prediction[1], 
         .lower = .lower - .lower[1], 
         .upper = .upper - .upper[1])


ggplot(test, aes(x = diff, y = model))+
  geom_pointintervalh()+
  geom_vline(xintercept = 0, linetype = "dashed")



require(tidyverse)
df1 <- data.frame(m1 = rnorm(100, 5), m2 = rnorm(100, 7), m3 = rnorm(100, 3)) %>%
  gather(model, prediction) %>%
  group_by(model) %>% 
  summarize(mean = mean(prediction)) %>%
  mutate(diff = mean - mean[1])

df1 %>% ggplot(aes(x = diff, y = model))+
  geom_point()+
  geom_vline(xintercept = 0, linetype = "dashed")

ggplot(aes(x = diff, y = model))+
  geom_point()+
  geom_vline(xintercept = 0, linetype = "dashed")

















































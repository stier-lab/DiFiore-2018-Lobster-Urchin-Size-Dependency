library(here)
library(tidyverse)
library(rstan)



gauss_model <- rstan::stan_model('code/STAN_models/allometric.stan')
#gauss_model <- rstan::stan_model('code/STAN_models/allometric_poisson.stan')


df <- read.table(here("data/cleaned","loburc_cleaned.csv"), header = T, sep = ",") %>%
  arrange(id, treatment) %>%
  drop_na(mc)

meta <- distinct(df, id, mc, mr) %>% 
  mutate(id = as.numeric(as.factor(id)))

stan_data = list("initial"= df$initial,
                 "killed" = df$killed,
                 "N" = length(df$initial), 
                 "mc" = meta$mc, 
                 "mr" = meta$mr,
                 "Nind" = length(unique(df$id)), 
                 "id" = as.numeric(as.factor(df$id))
) # named list

stanfit_gauss <- sampling(gauss_model, data = stan_data, chains = 3, thin = 3,
                          iter = 10000,
                          seed = 2131231, control = list(adapt_delta = 0.95))

shinystan::launch_shinystan(stanfit_gauss)

stanfit_gauss %>%
tidybayes::recover_types(df) %>%
  tidybayes::spread_draws(mualphaa, mualphah, beta1a, beta2a, beta1h, beta2h) %>%
  #tidybayes::sample_draws(n = 10000) %>%
  pivot_longer(cols = c(mualphaa:beta2h), names_to = "parameter", values_to = "estimate") %>%
  ggplot(aes(x=.iteration, y=estimate, color=as.factor(.chain))) +
  geom_line(alpha=0.5) +
  facet_grid(parameter~.chain, scale="free_y") +
  geom_smooth(method="loess") + labs(color="chain")

stanfit_gauss %>%
  tidybayes::recover_types(df) %>%
  tidybayes::spread_draws(mualphaa, mualphah, beta1a, beta2a, beta1h, beta2h) %>%
  #tidybayes::sample_draws(n = 10000) %>%
  pivot_longer(cols = c(mualphaa:beta2h), names_to = "parameter", values_to = "estimate") %>%
  group_by(parameter) %>%
  tidybayes::median_qi()

pairs(stanfit_gauss, pars = c("beta1a", "beta1h"), las = 1)

pairs(stanfit_gauss, pars = c("beta1h", "beta2h", "alphah"), las = 1)





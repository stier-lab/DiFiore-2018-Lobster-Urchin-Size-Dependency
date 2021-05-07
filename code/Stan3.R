library(here)
library(tidyverse)
library(rstan)

sink(here::here("code/STAN_models/FR_forloop.stan"))
cat("
    data {
    int<lower=1> Nind;
    int<lower=0> Ndensities;
    int killed[Nind, Ndensities];
    int initial[Nind, Ndensities];
    real<lower=0> mc[Nind];
    real<lower=0> mr[Nind];
    }
    
    parameters {
    real<lower=0> alphaa;
    real<lower=0> varalphaa;
    real beta1a; 
    real<lower=0> varbeta1a;
    real beta2a;
    real<lower=0> varbeta2a;  
    real<lower=0> alphah;
    real<lower=0> varalphah;
    real beta1h; 
    real<lower=0> varbeta1h;
    real beta2h;
    real<lower=0> varbeta2h;
    
    }
    
    transformed parameters{
    // likelihooh

    real<lower=0, upper=1> prob[Nind, Ndensities];
    
    for(i in 1:Nind){
    for(j in 1:Ndensities){
    prob[i,j] = 1/(1/(alphaa*mc[i]^beta1a*mr[i]^beta2a) + alphah*mc[i]^beta1h*mr[i]^beta2h*initial[i,j]);
    }
    
    }
    
    }
    
    
    model {
    
    //priors
    
    // Population level prior on attack rate intercept
    alphaa ~ normal(0, varalphaa);
    //log(alphaa) ~ normal(0, varalphaa);
    //target += -log(alphaa);
    varalphaa ~ exponential(1);
    
    // Allometric scaling exponents for attack rate
    beta1a ~ normal(0, varbeta1a);
    varbeta1a ~ gamma(1,3);
    
    beta2a ~ normal(0, varbeta2a);
    varbeta2a ~ gamma(1,3);
    
    //--------------------------------------------------------------------------
    
    //Population level prior on handling time intercept
    alphah ~ normal(0, varalphah);
    //log(alphah) ~ normal(0, varalphah);
    //target += -log(alphah);
    varalphah ~ exponential(1);
    
    // Allometric scaling exponents for handling time
    beta1h ~ normal(0, varbeta1h);
    varbeta1h ~ gamma(1,3);
    
    beta2h ~ normal(0, varbeta2h);
    varbeta2h ~ gamma(1,3);
    
    // likelihood
    
    for(i in 1:Nind){
    for(j in 1:Ndensities){
    killed[i,j] ~ binomial(initial[i,j], prob[i,j]);
    }
    }
    
    }",fill = TRUE)
sink()


gauss_model <- rstan::stan_model('code/STAN_models/FR_forloop.stan')

df <- read.table(here("data/cleaned","loburc_cleaned.csv"), header = T, sep = ",") %>%
  arrange(id, treatment) %>%
  drop_na(mc)

mat.initial <- df %>%
  arrange(id, initial) %>%
  group_by(id) %>%
  mutate(trial = 1:n()) %>%
  select(trial, id, initial) %>%
  pivot_wider(names_from = trial, values_from = initial) %>%
  drop_na()

mat.initial <- as.matrix(mat.initial[, 2:7])

mat.killed <- df %>%
  arrange(id, initial) %>%
  group_by(id) %>%
  mutate(trial = 1:n()) %>%
  select(trial, id, killed) %>%
  pivot_wider(names_from = trial, values_from = killed) %>% 
  drop_na()

meta <- distinct(df, id, mc, mr) %>% 
  filter(id %in% mat.killed$id) %>%
  mutate(id = as.numeric(as.factor(id)))

mat.killed <- as.matrix(mat.killed[, 2:7])

mc <- meta$mc
mr <- meta$mr

stan_data = list("initial"= mat.initial,
                 "killed" = mat.killed,
                 "mc" = mc, 
                 "mr" = mr,
                 "Nind" = dim(mat.initial)[1],
                 "Ndensities" = dim(mat.initial)[2]
) # named list

stanfit_gauss <- sampling(gauss_model, data = stan_data, chains = 4, thin = 10,
                          iter = 1*10^5, seed = 2131231, control = list(adapt_delta = 0.85))

temp <- stanfit_gauss %>%
  tidybayes::recover_types(df) %>%
  tidybayes::spread_draws(alphaa, alphah, beta1a, beta2a, beta1h, beta2h)
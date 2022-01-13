library(here)
library(tidyverse)
library(rstan)

sink(here::here("code/STAN_models/allometric.stan"))
cat("
    data {
    int N;
    int<lower=0> killed[N];
    int<lower=0> initial[N];
    int Nind;
    int<lower=1> id[N]; // this will need to be character???
    real<lower=0> mc[Nind];
    real<lower=0> mr[Nind];
    }
    
    parameters {
    real mualphaa;
    real<lower=0> varmualphaa;
    real alphaa[Nind];
    real<lower=0> varaind;
    real beta1a; 
    real<lower=0> varbeta1a;
    real beta2a;
    real<lower=0> varbeta2a;  
    real mualphah;
    real<lower=0> varmualphah;
    real alphah[Nind];
    real<lower=0> varhind;
    real beta1h; 
    real<lower=0> varbeta1h;
    real beta2h;
    real<lower=0> varbeta2h;
    
    }
    
    transformed parameters{
    // likelihood
    
    real<lower=0> a[Nind];
    real<lower=0> h[Nind]; 
    real loga[Nind];
    real logh[Nind];
    real<lower=0, upper=1> prob[N];
    
    for(i in 1:Nind){
    loga[i] = alphaa[i] + beta1a*log(mc[i]) + beta2a*log(mr[i]);
    a[i] = exp(loga[i]);
    logh[i] = alphah[i] + beta1h*log(mc[i]) + beta2h*log(mr[i]);
    h[i] = exp(logh[i]);
    }
    
    for(i in 1:N){
    //prob[i] = 1/(1/a[id[i]] + h[id[i]]*initial[i]);
    prob[i] = a[id[i]] / (1 + a[id[i]]*h[id[i]]*initial[i]);
    }
    
    }
    
    
    
    model {
    
    //priors
    
    // Population level prior on attack rate intercept
    mualphaa ~ normal(0, varmualphaa);
    //varmualphaa ~ uniform(0,10);
    
    // Individual level variation in alphaa
    for(i in 1:Nind){
    alphaa[i] ~ normal(mualphaa, varaind);
    }
    varaind ~ uniform(0,10);
    
    // Allometric scaling exponents for attack rate
    beta1a ~ normal(0.75, varbeta1a);
    varbeta1a ~ gamma(2,1);
    
    beta2a ~ normal(0.5, varbeta2a);
    varbeta2a ~ gamma(2,1);
    
    //--------------------------------------------------------------------------
    
    //Population level prior on handling time intercept
    mualphah ~ normal(0, varmualphah);
    varmualphah ~ uniform(0,10);
    
    //Individual level variation in handling time interpect 
    for(i in 1:Nind){
    alphah[i] ~ normal(mualphah, varhind);
    }
    varhind ~ uniform(0,10);
    
    // Allometric scaling exponents for handling time
    beta1h ~ normal(-0.75, varbeta1h);
    varbeta1h ~ gamma(2,1);
    
    beta2h ~ normal(0.5, varbeta2h);
    varbeta2h ~ gamma(2,1);
    
    // likelihood
    
    for(i in 1:N){
    killed[i] ~ binomial(initial[i], prob[i]);
    }
    
    
    }",fill = TRUE)
sink()


gauss_model <- rstan::stan_model('code/STAN_models/allometric.stan')

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
  tidybayes::sample_draws(n = 10000) %>%
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





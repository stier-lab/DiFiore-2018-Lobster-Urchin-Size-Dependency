library(tidybayes)
# This code I adapted from this posting here: https://discourse.mc-stan.org/t/difficulties-with-non-linear-hierarchical-model/20519. The poster gave me the idea to move to a poisson distribution with a non-centered parameterization.


sink("code/STAN_models/allometric_poissonNC.stan")
cat("
    data {
    int N;
    int<lower=0> killed[N];
    int<lower=0> initial[N];
    int Nind;
    int id[N]; 
    real<lower=0> mc[N];
    real<lower=0> mr[N];
    }
    
    parameters {

    // These are all the hyperparameters estimated across all individuals
    real alphaa;
    real<lower=0> varalphaa;
    real beta1a; 
    real<lower=0> varbeta1a;
    real beta2a;
    real<lower=0> varbeta2a;  
    real alphah;
    real<lower=0> varalphah;
    real beta1h; 
    real<lower=0> varbeta1h;
    real beta2h;
    real<lower=0> varbeta2h;
    
    // These are parameters associated with variation in attack rates and handling times between individuals. In other words they are the random effect of individual. 
    real<lower=0> sigmaa;
    vector[Nind] zza; 
    real<lower=0> sigmah; 
    vector[Nind] zzh;
    }
    
    transformed parameters{
    
    vector[N] a;
    vector[N] h; 
    vector[N] lambda; // This is the estimated number of prey eaten in any trial i. Lambda is the parameters for a poisson distribution.
    
    for(i in 1:N){
    a[i] = exp(alphaa + beta1a*log(mc[i]) + beta2a*log(mr[i]) + zza[id[i]]*sigmaa);
    h[i] = exp(alphah + beta1h*log(mc[i]) + beta2h*log(mr[i]) + zzh[id[i]]*sigmah);
    lambda[i] = a[i]*initial[i] / (1 + a[i]*h[i]*initial[i]);
    }
    

    }
    
    
    model {
    
    //priors
    
    // individual variation
    sigmaa ~ gamma(2,1);
    for(i in 1:Nind){
    zza[i] ~ normal(0, sigmaa);
    }
    
    sigmah ~ gamma(2,1);
    for(i in 1:Nind){
    zzh[i] ~ normal(0, sigmah);
    }
    
    // Population level prior on attack rate intercept
    alphaa ~ normal(0, varalphaa);
    varalphaa ~ gamma(1,1);    
    
    // Allometric scaling exponents for attack rate
    beta1a ~ normal(0.75, varbeta1a);
    varbeta1a ~ gamma(1,1);
    
    beta2a ~ normal(0.5, varbeta2a);
    varbeta2a ~ gamma(1,1);
    
    //--------------------------------------------------------------------------
    
    //Population level prior on handling time intercept
    alphah ~ normal(0, varalphah);
    varalphah ~ gamma(1,1);    
    
    // Allometric scaling exponents for handling time
    beta1h ~ normal(-0.75, varbeta1h);
    varbeta1h ~ gamma(1,1);
    
    beta2h ~ normal(0.5, varbeta2h);
    varbeta2h ~ gamma(1,1);
    
    // likelihood
    
    killed ~ poisson(lambda);
    
    }",fill = TRUE)
sink()


poissonNC <- rstan::stan_model('code/STAN_models/allometric_poissonNC.stan')

df <- read.table(here("data/cleaned","loburc_cleaned.csv"), header = T, sep = ",") %>%
  arrange(id, treatment) %>%
  drop_na(mc)

stan_data = list("initial"= df$initial,
                 "killed" = df$killed,
                 "N" = length(df$initial), 
                 "mc" = df$mc, 
                 "mr" = df$mr,
                 "Nind" = length(unique(df$id)), 
                 "id" = as.numeric(as.factor(df$id))
) # named list

stanfit <- sampling(poissonNC, data = stan_data, chains = 4, thin = 3,
                          iter = 20000,
                          seed = 123123, control = list(adapt_delta = 0.99))

formerge <- df %>% mutate(i = 1:n()) %>%
  select(i, id)

df.ind <- stanfit %>%
  recover_types(df) %>%
  gather_draws(a[i], h[i]) %>%
  left_join(formerge) %>%
  ungroup() %>%
  select(-i) %>%
  distinct()

df.pop <- stanfit %>%
  recover_types(df) %>%
  gather_draws(alphaa, beta1a, beta2a, alphah, beta1h, beta2h, varalphaa, varbeta1a, varbeta2a, varalphah, varbeta1h, varbeta2h)

write.csv(df.ind, here::here("data/cleaned/posteriors", "allometric_individualSTAN.csv"), row.names = F)
write.csv(df.pop, here::here("data/cleaned/posteriors", "allometric_populationSTAN.csv"), row.names = F)

#----------------------------------------------------
## Model evaluation
#----------------------------------------------------


shinystan::launch_shinystan(stanfit_gauss)

stanfit %>%
  tidybayes::recover_types(df) %>%
  tidybayes::spread_draws(alphaa, alphah, beta1a, beta2a, beta1h, beta2h) %>%
  #tidybayes::sample_draws(n = 10000) %>%
  pivot_longer(cols = c(alphaa:beta2h), names_to = "parameter", values_to = "estimate") %>%
  ggplot(aes(x=.iteration, y=estimate, color=as.factor(.chain))) +
  geom_line(alpha=0.5) +
  facet_grid(parameter~.chain, scale="free_y") +
  geom_smooth(method="loess") + labs(color="chain")


  
  
  
  
  
stanfit_gauss %>%
  tidybayes::recover_types(df) %>%
  tidybayes::spread_draws(h[i], a[i]) %>%
  tidybayes::median_qi()%>%
  select(h,a) %>%
  bind_cols(df) %>%
  ggplot(aes(x = initial, y = killed, group = id))+
  geom_point()+
  geom_line(aes(x = initial, y = holling2(N = initial, a = a, h =h, P = 1, T = 1)))+
  facet_wrap(~id)


stanfit_gauss %>%
  tidybayes::recover_types(df) %>%
  tidybayes::spread_draws(h[i], a[i]) %>%
  tidybayes::median_qi()%>%
  select(h,a) %>%
  bind_cols(df) %>%
  pivot_longer(cols = c(a,h)) %>%
  ggplot(aes(x = mc, y = value, group = name))+
  geom_point(aes(color = as.factor(mr)))+
  facet_wrap(~name)+
  scale_x_log10()+
  scale_y_log10()
  
  
  
  
  
  
  
  
  
  
  
  
  
  


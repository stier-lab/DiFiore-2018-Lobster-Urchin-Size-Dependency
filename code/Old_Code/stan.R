library(rstan)

sink(here::here("code", "FR.stan"))
cat("
    
    data {
    int<lower=1> N;
    int<lower=0> killed[N];
    int<lower=0> initial[N];
    int<lower=0> Ntreats;
    int<lower=1> Nind;
    int<lower=0> tind[Nind]; 
    int<lower=1> id[N];
    }
    
    parameters {
    real logit_a[Nind];
    real log_h[Nind];
    real t_logit_a[Ntreats];
    real t_log_h[Ntreats];
    real mu_logit_a;
    real mu_log_h;
    real<lower=0> sigma_t_a;
    real<lower=0> sigma_t_h;
    real<lower=0> sigma_ind_a;
    real<lower=0> sigma_ind_h;
    
    }

    transformed parameters{

    real<lower=0, upper=1> mu_a;
    real<lower=0> mu_h;
    real<lower=0, upper=1> t_a[Ntreats];
    real<lower=0> t_h[Ntreats];
    real<lower=0, upper=1> a[Nind];
    real<lower=0> h[Nind];
    vector[N] prob;
    

    //hyperpriors transformed

    mu_a = exp(mu_logit_a)/(1+exp(mu_logit_a));
    mu_h = exp(mu_log_h);

    //Treatment level priors transformed

    for(i in 1:Ntreats){
    t_a[i] = exp(t_logit_a[i])/(1+exp(t_logit_a[i]));
    t_h[i] = exp(t_log_h[i]);
    }

    //Individual level priors transformed

    for(i in 1:Nind){
    a[i] = exp(logit_a[i])/(1+exp(logit_a[i]));
    h[i] = exp(log_h[i]);
    }

    // likelihood

    for(i in 1:N){
    prob[i] = 1/(1/a[id[i]] + h[id[i]]*initial[i]);
    }

    
    }

    
    model {

    //priors
    
    //hyperpriors
    
    mu_logit_a ~ normal(0, 10);
    mu_log_h ~ normal(0, 10);
    
    // treatment level priors
    
    for(i in 1:Ntreats){
    t_logit_a[i] ~ normal(mu_logit_a, sigma_t_a);
    t_log_h[i] ~ normal(mu_log_h, sigma_t_h);
    }
    
    // individual level priors
    
    for(i in 1:Nind){
    logit_a[i] ~ normal(t_logit_a[tind[i]], sigma_ind_a);
    log_h[i] ~ normal(t_log_h[tind[i]], sigma_ind_h);
    }
    
    // Variance priors
    
    sigma_ind_a ~uniform(0,10);
    sigma_ind_h ~uniform(0,10);
    
    sigma_t_a ~ uniform(0,10);
    sigma_t_h ~ uniform(0,10);
    
    //likelihood
    
    
    for(i in 1:N){

    killed[i] ~ binomial(initial[i], prob[i]);
    }
    
    }
    
    
    
    ",fill = TRUE)
sink()

fr_model <- stan_model(here::here('code', 'FR.stan'))


df <- read.table(here("data/cleaned","loburc_cleaned.csv"), header = T, sep = ",") %>%
  arrange(id, treatment)

tind <- distinct(df, treatment, id) %>%
  arrange(id, treatment)

tind <- as.numeric(as.factor(as.vector(tind$treatment)))

id <- as.numeric(as.factor(as.vector(df$id)))


stan_data = list("initial"= df$initial,
                 "killed" = df$killed,
                 "id" = id ,
                 "Nind" = length(unique(df$id)),
                 "N" = length(df$initial), 
                 "tind" = tind, 
                 "Ntreats" = length(unique(df$treatment))
) # named list


fr_stan <- sampling(fr_model, data = stan_data, chains = 3,
                           iter = 10000, seed = 2131231)

temp <- as.data.frame(fr_stan)




# Ok, so I need to simulate a series of data sets that build in complexity. Specifically, I need a basic dataset of number of prey eaten as a function of number of prey offered, where it is simulated based on known attack rates and handling times. Then I need a structured data set, where at the very least individual predators forage at each prey density. 

# Simple data set 

fr <- function(a,h,N){
  a*N / (1+ a*h*N)
}

fr.prob <- function(a,h,N){
  (a*N / (1+ a*h*N))/N
}

a = 1
h = 0.1
N = seq(1,50, by = 1)

df <- data.frame(N = N, actual = fr(a = a, h = h, N = N))
df$predicted <- rbinom(n = length(N), size = N, prob = fr.prob(a = a, h = h, N = N))
plot(predicted ~ N, df, ylim = c(0, 14))
lines(actual~N, df, col = "darkred")
abline(a = 0, b = 1)
abline(h = 0)
abline(v = 1)


# Structured data set

# so this data set has 10 predators, each forage at six different densities of prey, attack rates, and handling times for each predator are determined by their body size. 

N = c(2, 5, 10, 20, 30, 50)
beta.a = 0.75
beta.h = -0.75
mass = seq(200, 2000, length.out = 10)
a0 = rlnorm(10, meanlog = log(1/max(mass^beta.a)/2), sdlog = 0.5) # so this is a little assbackwards but basically it ensures that most attack rates are less than 1, and all are positive. 
h0 = rlnorm(10, meanlog = 1, sdlog = 0.1)

a0 = 1/max(mass^beta.a)
h0 = 10

a <- vector()
h <- vector()

for(i in 1:10){
  a[i] <- a0[i]*mass[i]^beta.a
  h[i] <- h0[i]*mass[i]^beta.h
}

for(i in 1:10){
  a[i] <- a0*mass[i]^beta.a
  h[i] <- h0*mass[i]^beta.h
}

d <- expand.grid(N = N, ind = as.character(1:10))
d$a <- rep(a, each = 6)
d$h <- rep(h, each = 6)
d$mass <- rep(mass, each = 6)
d$killed <- fr(a = d$a, h = d$h, N = d$N)

ggplot(d, aes(x = N, y = killed))+
  geom_line(aes(color = as.factor(ind)))



out <- matrix(nrow=6, ncol = 10)
for(i in 1:10){
  out[,i] <- rbinom(n = length(N), size = N, prob = fr.prob(a = a[i], h = h[i], N = N))
}

df <- as.data.frame(out) %>% 
  mutate(N = N) %>%
  pivot_longer(V1:V10) %>%
  arrange(name, N) %>%
  rename(ind = name, data.sim = value) %>%
  separate(ind, into = c("junk", "ind"), sep = "(?<=[A-Za-z])(?=[0-9])") %>%
  right_join(d)

ggplot(df, aes(x = N))+
  geom_jitter(aes(x = N, y = data.sim, color = ind))+
  geom_line(aes(x = N, y = killed, color = ind))

df %>% pivot_longer(cols = c(a, h), names_to = "parameters", values_to = "estimate") %>% 
  ggplot(aes(x = mass, y = estimate))+
  geom_point(aes(color = ind))+
  geom_smooth(method = "lm", se = F)+
  facet_wrap(~parameters, scales = "free")+
  scale_x_log10()+
  scale_y_log10()

summary(lm(log(a) ~ log(mass), df))
summary(lm(log(h) ~ log(mass), df))



sink(here::here("code/STAN_models/FR_nopreymass.stan"))
cat("
    data {
    int<lower=1> Nind;
    int<lower=0> Ndensities;
    int killed[Nind, Ndensities];
    int initial[Nind, Ndensities];
    real<lower=0> mc[Nind];
    }
    
    parameters {
    real<lower=0> alphaa;
    real<lower=0> varalphaa;
    real beta1a; 
    real<lower=0> varbeta1a;
    real<lower=0> alphah;
    real<lower=0> varalphah;
    real beta1h; 
    real<lower=0> varbeta1h;
    
    }
    
    transformed parameters{
    // likelihood
    
    vector[Nind] loga;
    vector[Nind] logh;
    vector[Nind] a;
    vector[Nind] h;
    real<lower=0, upper=1> prob[Nind, Ndensities];
    
    for(i in 1:Nind){
    loga[i] = alphaa + beta1a*log(mc[i]);
    a[i] = exp(loga[i]);
    logh[i] = alphah + beta1h*log(mc[i]);
    h[i] = exp(logh[i]);
    
    for(j in 1:Ndensities){
    prob[i,j] = 1/(1/a[i] + h[i]*initial[i,j]);
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
    beta1a ~ normal(0.75, varbeta1a);
    varbeta1a ~ gamma(2,1);
    
    //--------------------------------------------------------------------------
    
    //Population level prior on handling time intercept
    alphah ~ normal(0, varalphah);
    //log(alphah) ~ normal(0, varalphah);
    //target += -log(alphah);
    varalphah ~ exponential(1);
    
    // Allometric scaling exponents for handling time
    beta1h ~ normal(-0.75, varbeta1h);
    varbeta1h ~ gamma(2,1);
    
    // likelihood
    
    for(i in 1:Nind){
    for(j in 1:Ndensities){
    killed[i,j] ~ binomial(initial[i,j], prob[i,j]);
    }
    }
    
    }",fill = TRUE)
sink()


gauss_model <- rstan::stan_model('code/STAN_models/FR_nopreymass.stan')







killed <- t(out)
initial <- t(matrix(rep(unique(df$N), 10), nrow = 6, ncol = 10))


stan_data = list("initial"= initial,
                 "killed" = killed,
                 "mc" = mass, 
                 "Nind" = dim(initial)[1],
                 "Ndensities" = dim(initial)[2]
) # named list

stanfit_gauss <- sampling(gauss_model, data = stan_data, chains = 4, thin = 10,
                          iter = 50000, seed = 2131231, control = list(adapt_delta = 0.85))



pairs(stanfit_gauss, pars = c("beta1a", "beta1h", "alphaa", "alphah"), las = 1)







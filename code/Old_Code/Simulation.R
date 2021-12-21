# Structured data set - no underlying size structure

# so this data set has 10 predators, each forage at six different densities of prey, and on one size class of prey. Attack rates, and handling times vary randomly between predators within a treatment. The goal is to determine if the model can accurately recover the attack rates and handling times for the population, treatments, and individuals. 

fr <- function(a,h,N){
  a*N / (1+ a*h*N)
}

fr.prob <- function(a,h,N){
  (a*N / (1+ a*h*N))/N
}

set.seed(10001)

N = c(2, 5, 10, 20, 30, 50)
N.ind = 10 # within each treatment, so 30 total lobsters
N.treats = 3
mass = seq(200, 2000, length.out = N.ind)

# The population level parameters

mu.logit.a <- rnorm(n = 1)
mu.a <- exp(mu.logit.a)/(1+exp(mu.logit.a))

mu.log.h <- rnorm(n = 1)
mu.h <- exp(max(min(mu.log.h, 20), -20))

t.a <- vector()
t.h <- vector()
t.logit.a <- vector()
t.log.h <- vector()

for(i in 1:N.treats){
  
  # a
  t.logit.a[i] <- rnorm(n = 1, mean = mu.logit.a, sd = 1)
  t.a[i] <- exp(t.logit.a[i])/(1+exp(t.logit.a[i]))
  
  # h
  t.log.h[i] <- rnorm(n = 1, mean = mu.log.h, sd = 1)
  t.h[i] <- exp(max(min(t.log.h[i],20),-20))
}


logit.a <- matrix(nrow = 10, ncol = 3)
log.h <- matrix(nrow = 10, ncol = 3)
a <- matrix(nrow = 10, ncol = 3)
h <- matrix(nrow = 10, ncol = 3)

for(j in 1:N.treats){
  for(i in 1:N.ind){
    # a
    logit.a[i, j] <- rnorm(n = 1, mean = t.logit.a[j], sd = 1 )
    a[i, j] <-  exp(logit.a[i, j])/(1+exp(logit.a[i, j]))
    
    # h
    log.h[i, j] <- rnorm(n = 1, mean = t.log.h[j], sd = 1)
    h[i, j] <- exp(max(min(log.h[i, j],20),-20))
  }
}

a <- as.data.frame(a)
names(a) <- c("t1", "t2", "t3")
a$ind <- as.character(1:10)
a <- pivot_longer(a, cols = t1:t3, names_to = "treatment", values_to = "estimate") %>% mutate(parameter = "a")

h <- as.data.frame(h)
names(h) <- c("t1", "t2", "t3")
h$ind <- as.character(1:10)
h <- pivot_longer(h, cols = t1:t3, names_to = "treatment", values_to = "estimate") %>% mutate(parameter = "h")

d <- bind_rows(a, h) %>% mutate(ind = as.character(ind), 
                                id = as.factor(paste(treatment, ind, sep = "_")))

df <- expand.grid(N = N, treatment = c("t1", "t2", "t3"), ind = as.character(1:N.ind)) %>% 
  arrange(treatment, ind, N) %>% 
  left_join(d) %>% 
  pivot_wider(names_from = parameter, values_from = estimate) %>%
  mutate(actual = fr(a = a, h = h, N = N), 
         prob = fr.prob(a = a, h = h, N = N))
df$sim <- rbinom(n = length(df$N), size = df$N, prob = df$prob)


df %>% ggplot(aes(x = N, y = sim))+
  geom_jitter(aes(group = id), show.legend= F)+
  geom_line(aes(y = actual, group = id), show.legend = F)+
  facet_wrap(~treatment)

# Fit the data using the hierarchical model

jagsscript = cat("

                 model{
                 # PRIORS
                 
                 # hyperprior
                 
                 mu.logit.a ~ dnorm(0, mu.tau.a) 
                 mu.a <- exp(mu.logit.a)/(1+exp(mu.logit.a))
                 
                 mu.log.h ~ dnorm(0, mu.tau.h)
                 mu.h <- exp(max(min(mu.log.h, 20), -20))
                 
                 
                 
                 # Treatment level piors
                 for(i in 1:Ntreats){
                 #a
                 t.logit.a[i] ~ dnorm(mu.logit.a, t.tau.a)
                 t.a[i] <- exp(t.logit.a[i])/(1+exp(t.logit.a[i]))

                 #h
                 t.log.h[i] ~ dnorm(mu.log.h, t.tau.h)
                 t.h[i] <- exp(max(min(t.log.h[i],20),-20))

                 }
                 
                 
                 # Individual level piors
                 for(i in 1:num.ind){
                 
                 # a
                 loga[i] ~ dnorm(t.logit.a[tind[i]], tau_int.a)
                 a[i] <- exp(loga[i])/(1+exp(loga[i]))
                 
                 # h
                 logh[i] ~ dnorm(t.log.h[tind[i]], tau_int.h)
                 h[i] <- exp(max(min(logh[i],10),-20))
                 }
                 
                 # functional response likelihood
                 
                  for(i in 1:n){
                 
                  prob[i] <- max(0.0001,min(0.9999,T/(1/a[id[i]] + h[id[i]]*initial[i])))
                 
                  killed[i] ~ dbin(prob[i],initial[i])
                 
                 }
                 
                 # Piors on variances for all levels
                 sigma_int.a ~dunif(0,10)
                 tau_int.a <- 1/(sigma_int.a*sigma_int.a)
                 sigma_int.h ~dunif(0,10)
                 tau_int.h <- 1/(sigma_int.h*sigma_int.h)
                 
                 sigma_t.a ~ dunif(0,10)
                 t.tau.a <- 1/(sigma_t.a*sigma_t.a)
                 sigma_t.h ~ dunif(0,10)
                 t.tau.h <- 1/(sigma_t.h*sigma_t.h)
                 
                 sigma_mu.a ~ dunif(0,10)
                 mu.tau.a <- 1/(sigma_mu.a*sigma_mu.a)
                 #mu.tau.a ~ dgamma(20, 40)
                 
                 sigma_mu.h ~ dunif(0,10)
                 mu.tau.h <- 1/(sigma_mu.h*sigma_mu.h)
                 #mu.tau.h ~ dgamma(20,40)
                 
                 }
                 ", file = here("code", "JAGs_models/heirarchical_jags.txt"))

model.loc=here("code","JAGs_models/heirarchical_jags.txt")
jags.params=c("a", "h", "t.a", "t.h", "mu.a", "mu.h", 
              "sigma_int.a", "sigma_int.h", "sigma_t.a", "sigma_t.h", 
              "loga", "logh", "t.logit.a", "t.log.h", "sigma_mu.a", "sigma_mu.h", "mu.tau.a", "mu.tau.h"
)


tind <- distinct(df, treatment, id) %>%
  arrange(id)

tind <- as.factor(as.vector(tind$treatment))

jags.data = list("initial"= df$N,
                 "killed" = df$sim,
                 "id" = df$id,
                 "num.ind" = length(unique(df$id)),
                 "n" = length(df$N), 
                 "tind" = tind, 
                 "T" = 1, 
                 "Ntreats" = length(unique(df$treatment))
) # named list

n.chains = 3
n.burnin = 250000/4
n.thin = 10
n.iter = 500000/4
model = R2jags::jags(jags.data,parameters.to.save=jags.params,inits=NULL,
                     model.file=model.loc, n.chains = n.chains, n.burnin=n.burnin,n.thin=n.thin, n.iter=n.iter, DIC=TRUE)


View(distinct(df, id, a, h))


df.treat <- coda::as.mcmc(model) %>%
  recover_types(df) %>%
  spread_draws(t.a[treatment], t.h[treatment]) 

df.treat %>%
  group_by(treatment) %>%
  mean_qi()

df.ind <- coda::as.mcmc(model) %>%
  recover_types(df) %>%
  spread_draws(a[id], h[id]) %>%
  pivot_longer(cols = c(a,h), names_to = "parameter", values_to = "prediction")

df.ind %>%
  group_by(id) %>%
  median_qi() %>% 
  View()

a %>% bind_rows(h) %>% 
  mutate(id = paste(treatment, ind, sep = "_")) %>% 
  left_join(df.ind) %>%
  rename(actual = estimate) %>%
  group_by(id, parameter, actual) %>%
  median_qi(prediction) %>% 
  ggplot(aes(x = id))+
  geom_pointinterval(aes(y = prediction, ymin = .lower, ymax = .upper))+
  geom_point(aes(y = actual), col = "red")+
  facet_wrap(~parameter, scales = "free")+
  coord_flip()


s1 = 1
s2 = 4

temp <- rbeta(n = 10000, shape1 = s1, shape2 = s2)
hist(temp)
mean(temp)

s1/(s1+s2)





















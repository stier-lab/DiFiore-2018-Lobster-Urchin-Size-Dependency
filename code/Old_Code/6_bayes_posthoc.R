#-------------------------------------------
## Setup
#-------------------------------------------

library(here)
source(here("code", "Base_functions/Functions.R"))
source(here("code", "1_setup.R"))

#-------------------------------------------
## Get data
#-------------------------------------------

meta <- read.csv(here::here("data/", "lob-metadata.csv"))
df <- read.csv(here::here("data/cleaned/posteriors", "posteriors_individuals.csv")) %>% 
  left_join(meta) %>%
  filter(id != "N07")



ltow <- function(udiam){0.000592598*udiam^2.872636198*1.01}
ltow(c(10,30,50,70))


#---------------------------------------------------------------------------
## Construct a model that incorporates uncertainty in the h estimates
#---------------------------------------------------------------------------

d2 <- read.csv(here::here("data/cleaned/posteriors", "posteriors_individuals.csv")) %>% 
  left_join(meta) %>%
  group_by(id) %>%
  sample_draws(500) %>%
  filter(id != "N07")

jagsscript = cat("
                 
                 model{
                 
                 for(i in 1:n){
                 
                 y[i] ~ dnorm(mu[i], tauy)
                 mu[i] <- alpha + beta1*mc[i] + beta2*mr[i]
                 
                 }
                 
                 #priors
                 alpha ~ dnorm(0, 0.01)
                 beta1 ~ dnorm(-0.75, tau.beta1)
                 beta2 ~ dnorm(0.5, tau.beta2)
                 
                 # priors on variance
                 sigmay ~ dunif(0,10)
                 tauy <- 1/(sigmay*sigmay)
                 
                 sigma.beta1 ~ dgamma(2,2)
                 tau.beta1 <- 1/(sigma.beta1*sigma.beta1)
                 
                 sigma.beta2 ~ dgamma(2,2)
                 tau.beta2 <- 1/(sigma.beta2*sigma.beta2)
                 
                 }
                 ", file = here("code", "post-hocH_wparameteruncertainty.txt"))



model.loc=here("code","post-hocH_wparameteruncertainty.txt")
jags.params=c("alpha", "beta1", "beta2")


jags.data = list("y" = log(d2$h),
                 "mc" = log(d2$mc), 
                 "mr" = log(d2$mr),
                 "n" = length(d2$id)
) # named list

n.chains = 3
n.burnin = 10000
n.thin = 3
n.iter = 30000
model_wpuncert = R2jags::jags(jags.data,parameters.to.save=jags.params,inits=NULL,
                              model.file=model.loc, n.chains = n.chains, n.burnin=n.burnin,n.thin=n.thin, n.iter=n.iter, DIC=TRUE)

df.hwpuncert <- mcmcplots::as.mcmc.rjags(model_wpuncert) %>%
  spread_draws(alpha, beta1, beta2) %>%
  mutate(model = "model_wpuncert")


write.csv(df.hwpuncert, here::here("data/cleaned/posteriors", "posteriors_posthoc_h.csv"), row.names = F)





#---------------------------------------------------------------------------
## Construct a model that incorporates uncertainty in the a estimates
#---------------------------------------------------------------------------

beta1.theory <- mean(c(0.58,0.92))
beta2.theory <- mean(c(0.33,0.66))

  jagsscript = cat("
                 
                 model{
                 
                 for(i in 1:n){
                 
                 y[i] ~ dnorm(mu[i], tauy)
                 mu[i] <- alpha + beta1*mc[i] + beta2*mr[i]
                 
                 }
                 
                 #priors
                 alpha ~ dnorm(0, 0.01)
                 beta1 ~ dnorm(beta1.theory, tau.beta1)
                 beta2 ~ dnorm(beta2.theory, tau.beta2)
                 
                 # priors on variance
                 sigmay ~ dunif(0,10)
                 tauy <- 1/(sigmay*sigmay)
                 
                 sigma.beta1 ~ dgamma(2,2)
                 tau.beta1 <- 1/(sigma.beta1*sigma.beta1)
                 
                 sigma.beta2 ~ dgamma(2,2)
                 tau.beta2 <- 1/(sigma.beta2*sigma.beta2)
                 }
                 
                 ", file = here("code", "post-hocA_wparameteruncertainty.txt"))

model.loc=here("code","post-hocA_wparameteruncertainty.txt")
jags.params=c("alpha", "beta1", "beta2")


jags.data = list("y" = log(d2$a),
                 "mc" = log(d2$mc), 
                 "mr" = log(d2$mr),
                 "n" = length(d2$id), 
                 "beta1.theory" = beta1.theory, 
                 "beta2.theory" = beta2.theory
) # named list

n.chains = 3
n.burnin = 10000
n.thin = 3
n.iter = 30000
model_wpuncert.a = R2jags::jags(jags.data,parameters.to.save=jags.params,inits=NULL,
                              model.file=model.loc, n.chains = n.chains, n.burnin=n.burnin,n.thin=n.thin, n.iter=n.iter, DIC=TRUE)


df.awpuncert <- mcmcplots::as.mcmc.rjags(model_wpuncert.a) %>%
  spread_draws(alpha, beta1, beta2) %>%
  mutate(model = "model_wpuncert")


write.csv(df.awpuncert, here::here("data/cleaned/posteriors", "posteriors_posthoc_a.csv"), row.names = F)

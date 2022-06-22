library(here)
source(here("code", "1_setup.R"))
source(here("code", "Base_functions/functions.R"))


#####################################
## Get data
#####################################


df <- read.table(here("data/cleaned","loburc_cleaned.csv"), header = T, sep = ",") %>%
  arrange(id, treatment)


dunn.h <- 0.741 * 24 #days
log.dh <- log(dunn.h)

# tank area in Dunn 
area.dunn <- pi*(2.18/2)^2 + 0.91*2.18
dunn.a <- 0.194 / area.dunn / 
  24 / #convert to hours
  tsize #convert to units of our tanks ~ 2 m2.
logit.da <- rstanarm::logit(dunn.a)

##################################################################
## the model
##################################################################

jagsscript = cat("

model{
# PRIORS

  # hyperprior
  
      mu.logit.a ~ dnorm(logit.da, mu.tau.a) # Dunn prior on attack rate
      mu.a <- exp(mu.logit.a)/(1+exp(mu.logit.a))
  
      mu.log.h ~ dnorm(log.dh, mu.tau.h) # Dunn prior on h
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

          #sigma_mu.a ~ dunif(0,10)
          #mu.tau.a <- 1/(sigma_mu.a*sigma_mu.a)
          mu.tau.a ~ dgamma(20, 40)

          #sigma_mu.h ~ dunif(0,10)
          #mu.tau.h <- 1/(sigma_mu.h*sigma_mu.h)
          mu.tau.h ~ dgamma(20,40)

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

jags.data = list("initial"= df$initial,
                 "killed" = df$killed,
                 "id" = as.factor(df$id),
                 "num.ind" = length(unique(df$id)),
                 "n" = length(df$initial), 
                 "tind" = tind, 
                 "T" = 48, 
                 "Ntreats" = length(unique(df$treatment)), 
                 "log.dh" = log.dh, 
                 "logit.da" = logit.da
) # named list

n.chains = 3
n.burnin = 250000
n.thin = 10
n.iter = 500000
model = R2jags::jags(jags.data,parameters.to.save=jags.params,inits=NULL,
                     model.file=model.loc, n.chains = n.chains, n.burnin=n.burnin,n.thin=n.thin, n.iter=n.iter, DIC=TRUE)


df.treat <- as.mcmc(model) %>%
  recover_types(df) %>%
  spread_draws(t.a[treatment], t.h[treatment])

df.ind <- as.mcmc(model) %>%
  recover_types(df) %>%
  spread_draws(a[id], h[id])

df.pop <- as.mcmc(model) %>%
  recover_types(df) %>%
  spread_draws(mu.a, mu.h, mu.tau.a, mu.tau.h, sigma_int.a, sigma_int.h, sigma_t.a, sigma_t.h)

write.csv(df.pop, here::here("data/cleaned/posteriors", "posteriors_population.csv"), row.names = F)
write.csv(df.treat, here::here("data/cleaned/posteriors", "posteriors_treatments.csv"), row.names = F)
write.csv(df.ind, here::here("data/cleaned/posteriors", "posteriors_individuals.csv"), row.names = F)
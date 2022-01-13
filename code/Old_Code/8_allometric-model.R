jagsscript = cat("
model{
                 
# Piors

#-----------------------------------------------------------------------------

  # Population level prior on attack rate intercept
    mu.alpha.a ~ dnorm(-4.5, var.mu.alpha.a)
    #mu.alpha.a ~ dnorm(0, var.mu.alpha.a)
    var.mu.alpha.a ~ dunif(0,5)

  # Individual level variation in alpha.a
    for(i in 1:Nind){
      alpha.a[i] ~ dnorm(mu.alpha.a, var.a.ind)
      }
    var.a.ind ~ dunif(0,10)

  # Allometric scaling exponents for attack rate
    beta1.a ~ dnorm(0.75, var.beta1.a)
    var.beta1.a ~ dgamma(2,1)

    beta2.a ~ dnorm(0.5, var.beta2.a)
    var.beta2.a ~ dgamma(2,1)
    
#------------------------------------------------------------------------------

  # Population level prior on handling time intercept
    mu.alpha.h ~ dnorm(11, var.mu.alpha.h)
    #mu.alpha.h ~ dnorm(0, var.mu.alpha.h)
    var.mu.alpha.h ~ dunif(0,5)
   
  # Individual level variation in handling time interpect 
    for(i in 1:Nind){
      alpha.h[i] ~ dnorm(mu.alpha.h, var.h.ind)
      }
    var.h.ind ~ dunif(0,10)

  # Allometric scaling exponents for handling time
    beta1.h ~ dnorm(-0.75, var.beta1.h)
    var.beta1.h ~ dgamma(10,1)

    beta2.h ~ dnorm(0.5, var.beta2.h)
    var.beta2.h ~ dgamma(2,1)

    #var.beta1.h <- 10 # what happens when you crank down the variance on the prior, such that all of the probability mass is right around -0.75? The result is that the data is still able to overwhelm the prior!!! 

#------------------------------------------------------------------------------

# functional response likelihood

    for(i in 1:Nind){
      log(a[i]) <- alpha.a[i] + beta1.a*log(mc[i]) + beta2.a*log(mr[i])
      log(h[i]) <- alpha.h[i] + beta1.h*log(mc[i]) + beta2.h*log(mr[i])
      }

    for(i in 1:n){
      killed[i] ~ dbin(prob[i],initial[i])
      prob[i] <- max(0.0001,min(0.9999,T/(1/a[id[i]] + h[id[i]]*initial[i])))
      }
#----------------------------------------------------------------------------
                 
                 }
                 ", file = here("code/JAGs_models", "allometric_model.txt"))

model.loc=here("code/JAGs_models","allometric_model.txt")
jags.params=c("alpha.a", "beta1.a", "beta2.a", "alpha.h", "beta1.h", "beta2.h", "mu.alpha.h", "mu.alpha.a", "var.mu.alpha.a", "var.mu.alpha.h", "var.beta1.a", "var.beta2.a", "var.beta1.h", "var.beta2.h", "var.a.ind", "var.h.ind")

df <- read.table(here("data/cleaned","loburc_cleaned.csv"), header = T, sep = ",") %>%
  mutate(id = as.factor(id)) %>%
  arrange(id, treatment) %>%
  drop_na(mc)

df$id <- droplevels(df$id)

meta <- distinct(df, id, mc, mr)

jags.data = list("initial"= df$initial,
                 "killed" = df$killed,
                 "n" = length(df$initial), 
                 "T" = 48,
                 "mc" = meta$mc,
                 "mr" = meta$mr,
                 "Nind" = length(unique(df$id)), 
                 "id" = df$id
) # named list

n.chains = 3
n.burnin = 100 * 1000
n.thin = 10
n.iter = 250 * 1000
model = R2jags::jags(jags.data,parameters.to.save=jags.params,inits=NULL,
                     model.file=model.loc, n.chains = n.chains, n.burnin=n.burnin,n.thin=n.thin, n.iter=n.iter, DIC=TRUE)
total.size <- n.chains*n.iter/n.thin - n.chains*n.burnin/n.thin

# with uninformative priors on mu.X.0 parameters 
# mu.alpha.a      -4.034   1.256  -6.865  -4.785  -3.886  -3.154  -1.759
# mu.alpha.h       9.960   1.387   7.135   9.009  10.001  10.993  12.486 



df.ind <- as.mcmc(model) %>%
  recover_types(df) %>%
  spread_draws(alpha.a[id], alpha.h[id])

df.pop <- as.mcmc(model) %>%
  recover_types(df) %>%
  spread_draws(mu.alpha.h, mu.alpha.a, beta1.a, beta2.a, beta1.h, beta2.h, var.a.ind, var.h.ind, var.mu.alpha.a, var.mu.alpha.h, var.beta1.a, var.beta2.a, var.beta1.h, var.beta2.h)

write.csv(df.ind, here::here("data/cleaned/posteriors", "allometric_individual.csv"), row.names = F)
write.csv(df.pop, here::here("data/cleaned/posteriors", "allometric_population.csv"), row.names = F)

#----------------------------------------------------
## Model evaluation
#----------------------------------------------------

df.pop %>%
  sample_draws(n = 10000) %>%
  select(-c(var.a.ind:var.beta2.h)) %>%
  pivot_longer(cols = c(mu.alpha.h:beta2.h), names_to = "parameter", values_to = "estimate") %>%
  ggplot(aes(x=.iteration, y=estimate, color=as.factor(.chain))) +
  geom_line(alpha=0.5) +
  facet_grid(parameter~.chain, scale="free_y") +
  geom_smooth(method="loess") + labs(color="chain")





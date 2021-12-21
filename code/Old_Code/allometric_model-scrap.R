jagsscript = cat("
model{
                 
                 # Piors
                 
                 #-----------------------------------------------------------------------------
                 
                 # Population level prior on attack rate intercept
                 log.alpha.a ~ dnorm(0, 1)
                 #var.mu.alpha.a ~ dunif(0,5)
                 alpha.a <- exp(max(min(log.alpha.a, 20), -20))
                 
                 # Allometric scaling exponents for attack rate
                 beta1.a ~ dnorm(0.5, var.beta1.a)
                 var.beta1.a ~ dgamma(2,1)
                 
                 beta2.a ~ dnorm(0.5, var.beta2.a)
                 var.beta2.a ~ dgamma(2,1)
                 
                 
                 #------------------------------------------------------------------------------
                 
                 # Population level prior on handling time intercept
                 log.alpha.h ~ dnorm(0, 1)
                 #var.mu.alpha.h ~ dunif(0,5)
                 alpha.h <- exp(max(min(log.alpha.h, 20), -20))
                 
                 # Allometric scaling exponents for handling time
                 beta1.h ~ dnorm(-0.75, var.beta1.h)
                 var.beta1.h ~ dgamma(2,1)
                 

                 beta2.h ~ dnorm(0.5, var.beta2.h)
                 var.beta2.h ~ dgamma(2,1)
                 
                 #var.beta1.h <- 10 # what happens when you crank down the variance on the prior, such that all of the probability mass is right around -0.75? The result is that the data is still able to overwhelm the prior!!! 
                 
                 #------------------------------------------------------------------------------
                 
                 # functional response likelihood
                 
                 for(i in 1:Nind){
                 log(a[i]) <- log.alpha.a + beta1.a*log(mc[i]) + beta2.a*log(mr[i])
                 log(h[i]) <- log.alpha.h + beta1.h*log(mc[i]) + beta2.h*log(mr[i])
                 }
                 
                 for(i in 1:n){
                 killed[i] ~ dbin(prob[i],initial[i])
                 prob[i] <- max(0.0001,min(0.9999,T/(1/a[id[i]] + h[id[i]]*initial[i])))
                 }
                 #----------------------------------------------------------------------------
                 
                 }
                 ", file = here("code/JAGs_models", "allometric_model-scrap.txt"))

model.loc=here("code/JAGs_models","allometric_model-scrap.txt")
jags.params=c("alpha.a", "beta1.a", "beta2.a", "alpha.h", "beta1.h", "beta2.h", "mu.alpha.h", "mu.alpha.a", "var.mu.alpha.a", "var.mu.alpha.h", "var.beta1.a", "var.beta2.a", "var.beta1.h", "var.beta2.h", "var.a.ind", "var.h.ind", "log.alpha.a", "log.alpha.h")

df <- read.table(here("data/cleaned","loburc_cleaned.csv"), header = T, sep = ",") %>%
  mutate(id = as.factor(id)) %>%
  arrange(id, treatment) %>%
  drop_na(mc)

df$id <- droplevels(df$id)

meta <- distinct(df, id, mc, mr)

jags.data = list("initial"= df$initial,
                 "killed" = df$killed,
                 "n" = length(df$initial), 
                 "T" = 1,
                 "mc" = meta$mc,
                 "mr" = meta$mr,
                 "Nind" = length(unique(df$id)), 
                 "id" = df$id
) # named list

n.chains = 3
n.burnin = 100 * 1000 /4
n.thin = 10
n.iter = 250 * 1000 /4
model = R2jags::jags(jags.data,parameters.to.save=jags.params,inits=NULL,
                     model.file=model.loc, n.chains = n.chains, n.burnin=n.burnin,n.thin=n.thin, n.iter=n.iter, DIC=TRUE)
total.size <- n.chains*n.iter/n.thin - n.chains*n.burnin/n.thin

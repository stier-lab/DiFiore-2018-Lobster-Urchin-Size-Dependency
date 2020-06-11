jagsscript = cat("
model{
                 
# Piors

#-----------------------------------------------------------------------------

# Population level prior on attack rate intercept
mu.alpha.a ~ dnorm(0, var.mu.alpha.a)
var.mu.alpha.a ~ dunif(0,10)

# Individual level variation in alpha.a    
for(i in 1:Nind){
alpha.a[i] ~ dnorm(mu.alpha.a, var.a.ind)
}
var.a.ind ~ dunif(0.5,2)

# Allometric scaling exponents for attack rate
beta1.a ~ dnorm(0.58, var.beta1.a)
#var.beta1.a ~ dgamma(100,1)

beta2.a ~ dnorm(0.33, var.beta2.a)
var.beta2.a ~ dgamma(100,1)

#------------------------------------------------------------------------------

# Population level prior on handling time intercept
mu.alpha.h ~ dnorm(0, var.mu.alpha.h)
var.mu.alpha.h ~ dunif(0,10)

# Individual level variation in handling time interpect 
for(i in 1:Nind){
alpha.h[i] ~ dnorm(mu.alpha.h, var.h.ind)
}
var.h.ind ~ dunif(0,10)

# Allometric scaling exponents for handling time
beta1.h ~ dnorm(-0.75, var.beta1.h)
var.beta1.h ~ dgamma(100,1)

beta2.h ~ dnorm(0.5, var.beta2.h)
var.beta2.h ~ dgamma(100,1)

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
", file = here("code/JAGs_models", "MTE_model.txt"))

model.loc=here("code/JAGs_models","MTE_model.txt")
jags.params=c("alpha.a", "beta1.a", "beta2.a", "alpha.h", "beta1.h", "beta2.h", "mu.alpha.h", "mu.alpha.a", "var.mu.alpha.a", "var.mu.alpha.h", "var.beta1.a", "var.beta2.a", "var.beta1.h", "var.beta2.h", "var.a.ind", "var.h.ind")

df <- read.table(here("data/cleaned","loburc_cleaned.csv"), header = T, sep = ",") %>%
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
n.burnin = 100 * 1000 /4
n.thin = 10
n.iter = 250 * 1000 /4
model = R2jags::jags(jags.data,parameters.to.save=jags.params,inits=NULL,
                     model.file=model.loc, n.chains = n.chains, n.burnin=n.burnin,n.thin=n.thin, n.iter=n.iter, DIC=TRUE)



df.mte <- as.mcmc(model) %>%
  recover_types(df) %>%
  spread_draws(mu.alpha.h, mu.alpha.a, beta1.a, beta2.a, beta1.h, beta2.h, var.a.ind, var.h.ind, var.mu.alpha.a, var.mu.alpha.h, var.beta1.a, var.beta2.a, var.beta1.h, var.beta2.h)

write.csv(df.mte, here::here("data/cleaned/posteriors", "mte_population.csv"), row.names = F)



jagsscript = cat("
model{
                 
                 # Piors
                 
                 #-----------------------------------------------------------------------------
                 
                 # Population level prior on attack rate intercept
                 alpha.a ~ dnorm(-4.5, var.alpha.a)
                 var.alpha.a ~ dgamma(2,1)
                 
                 # Allometric scaling exponents for attack rate
                 beta1.a ~ dnorm(0.75, var.beta1.a)
                 var.beta1.a ~ dgamma(2,1)
                 
                 beta2.a ~ dnorm(0.5, var.beta2.a)
                 var.beta2.a ~ dgamma(2,1)
                 
                 #------------------------------------------------------------------------------
                 
                 # Population level prior on handling time intercept
                 alpha.h ~ dnorm(11, var.alpha.h)
                 var.alpha.h ~ dgamma(2,1)
                 
                 # Allometric scaling exponents for handling time
                 beta1.h ~ dnorm(-0.75, var.beta1.h)
                 var.beta1.h ~ dgamma(10,1)
                 
                 beta2.h ~ dnorm(0.5, var.beta2.h)
                 var.beta2.h ~ dgamma(2,1)

                 #------------------------------------------------------------------------------
                 for(i in 1:Nind){
                 zz.a ~ dnorm(0, sigma.a)
                 }
                 sigma.a ~ dgamma(2,1)
                 
                 for(i in 1:Nind){
                 zz.h ~ dnorm(0, sigma.h)
                 }
                 sigma.h ~ dgamma(2,1)
                 
                 
                 #------------------------------------------------------------------------------
                 
                 # functional response likelihood
                 
                 for(i in 1:Nind){
                 log(a[i]) <- alpha.a[i] + beta1.a*log(mc[i]) + beta2.a*log(mr[i]) + zz.a[i]
                 log(h[i]) <- alpha.h[i] + beta1.h*log(mc[i]) + beta2.h*log(mr[i]) + zz.h[i]
                 }
                 
                 for(i in 1:n){
                 killed[i] ~ dbin(prob[i],initial[i])
                 prob[i] <- max(0.0001,min(0.9999,T/(1/a[id[i]] + h[id[i]]*initial[i])))
                 }
                 #----------------------------------------------------------------------------
                 
                 }
                 ", file = here("code/JAGs_models", "allometric_altrandom.txt"))


model.loc=here("code/JAGs_models","allometric_model.txt")
jags.params=c("alpha.a", "beta1.a", "beta2.a", "alpha.h", "beta1.h", "beta2.h", "var.alpha.a", "var.alpha.h", "var.beta1.a", "var.beta2.a", "var.beta1.h", "var.beta2.h", "zz.a", "zz.h", "sigma.a", "sigma.h")

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
n.burnin = 100 * 1000 / 5
n.thin = 10
n.iter = 250 * 1000 / 5
model = R2jags::jags(jags.data,parameters.to.save=jags.params,inits=NULL,
                     model.file=model.loc, n.chains = n.chains, n.burnin=n.burnin,n.thin=n.thin, n.iter=n.iter, DIC=TRUE)



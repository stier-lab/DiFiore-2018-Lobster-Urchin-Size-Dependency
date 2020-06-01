source(here::here("code", "setup.R"))

post.a <- read.csv(here::here("data/cleaned", "posteriors_posthoc_a.csv"))
post.h <- read.csv(here::here("data/cleaned", "posteriors_posthoc_h.csv"))

median.a <- post.a %>% group_by(model, parameter) %>%
  summarise(median = median(estimate))

median.h <- post.h %>% group_by(model, parameter) %>%
  summarise(median = median(estimate))

jagsscript = cat("

model{
  
  # Piors

mu.alpha.a ~ dnorm(exp(-4.5), 0.01)T(0,)

for(i in 1:Nind){
  alpha.a[i] ~ dnorm(mu.alpha.a, var1)
}
  var1 ~ dunif(0, 1)
  beta1.a ~ dnorm(0, 0.01)
  beta2.a ~ dnorm(0, 0.01)

mu.alpha.h ~ dnorm(exp(11), 0.01)T(0,)

for(i in 1:Nind){
  alpha.h[i] ~ dnorm(mu.alpha.h, var2)
}
  var2 ~ dunif(0,1)
  beta1.h ~ dnorm(0, 0.01)
  beta2.h ~ dnorm(0, 0.01)


  # functional response likelihood
  
for(i in 1:n){

killed[i] ~ dbin(prob[i],initial[i])

prob[i] <- max(0.0001,min(0.9999,T/(1/(alpha.a[id[i]] * mc[i]^beta1.a * mr[i]^beta2.a) + (alpha.h[id[i]] * mc[i]^beta1.h * mr[i]^beta2.h)*initial[i])))



}
  
  
}
", file = here("code", "test2.txt"))

model.loc=here("code","test2.txt")
jags.params=c("alpha.a", "beta1.a", "beta2.a", "alpha.h", "beta1.h", "beta2.h", "var1", "var2", "mu.alpha.a", "mu.alpha.h")

df <- read.table(here("data/cleaned","loburc_cleaned.csv"), header = T, sep = ",") %>%
  arrange(id, treatment) %>%
  drop_na(mc)

df$id <- droplevels(df$id)

jags.data = list("initial"= df$initial,
                 "killed" = df$killed,
                 "n" = length(df$initial), 
                 "T" = 48, 
                 "mc" = df$mc, 
                 "mr" = df$mr, 
                 "Nind" = length(unique(df$id)), 
                 "id" = df$id
) # named list

n.chains = 3
n.burnin = 100 * 1000 /2
n.thin = 3
n.iter = 250 * 1000 /2
model = R2jags::jags(jags.data,parameters.to.save=jags.params,inits=NULL,
                     model.file=model.loc, n.chains = n.chains, n.burnin=n.burnin,n.thin=n.thin, n.iter=n.iter, DIC=TRUE)













set.seed(1001)


jagsscript = cat("

                 model{
                 
                 # Piors
                 
                 
                 alpha.a ~ dnorm(exp(-4.5), var1)
                 var1 ~ dunif(0, 1)
                 beta1.a ~ dnorm(0, var2)
                 var2 ~ dunif(0,1)
                 beta2.a ~ dnorm(0, var3)
                 var3 ~ dunif(0,1)
                 

                 alpha.h ~ dnorm(exp(11), var4)
                 var4 ~ dunif(0,1)
                 beta1.h ~ dnorm(0, var5)
                 var5 ~ dunif(0,1)
                 beta2.h ~ dnorm(0, var6)
                 var6 ~ dunif(0,1)
                 
                 
                 # functional response likelihood
                 
                 for(i in 1:n){
                 
                 killed[i] ~ dbin(prob[i],initial[i])
                 
                 prob[i] <- max(0.0001,min(0.9999,T/(1/(alpha.a * mc[i]^beta1.a * mr[i]^beta2.a) + (alpha.h * mc[i]^beta1.h * mr[i]^beta2.h)*initial[i])))
                 
                 }
                 
                 
                 }
                 ", file = here("code", "test.txt"))



model.loc=here("code","test.txt")
jags.params=c("alpha.a", "beta1.a", "beta2.a", "alpha.h", "beta1.h", "beta2.h", "var1", "var2", "var3", "var4", "var5", "var6")

df <- read.table(here("data/cleaned","loburc_cleaned.csv"), header = T, sep = ",") %>%
  arrange(id, treatment) %>%
  drop_na(mc)

df$id <- droplevels(df$id)

jags.data = list("initial"= df$initial,
                 "killed" = df$killed,
                 "n" = length(df$initial), 
                 "T" = 48, 
                 "mc" = df$mc, 
                 "mr" = df$mr
) # named list

n.chains = 3
n.burnin = 10 * 1000
n.thin = 3
n.iter = 25 * 1000
model = R2jags::jags(jags.data,parameters.to.save=jags.params,inits=NULL,
                     model.file=model.loc, n.chains = n.chains, n.burnin=n.burnin,n.thin=n.thin, n.iter=n.iter, DIC=TRUE)


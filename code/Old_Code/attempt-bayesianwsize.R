source(here::here("code", "1_setup.R"))

post.a <- read.csv(here::here("data/cleaned/posteriors", "posteriors_posthoc_a.csv"))
post.h <- read.csv(here::here("data/cleaned/posteriors", "posteriors_posthoc_h.csv"))

median.a <- post.a %>% group_by(model, parameter) %>%
  summarise(median = median(estimate))

median.h <- post.h %>% group_by(model, parameter) %>%
  summarise(median = median(estimate))

jagsscript = cat("

model{
  
  # Piors

mu.alpha.a ~ dnorm(-4.5, 1)

for(i in 1:Nind){
alpha.a[i] ~ dnorm(mu.alpha.a, var5)
}
  var5 ~ dunif(0.5,2)
  beta1.a ~ dnorm(0, 1)
  beta2.a ~ dnorm(0, 1)

mu.alpha.h ~ dnorm(11, var1)
var1 ~ dunif(0,10)

for(i in 1:Nind){
  alpha.h[i] ~ dnorm(mu.alpha.h, var2)
}
  var2 ~ dunif(0,10)
  beta1.h ~ dnorm(0, var3)
  var3 ~ dunif(0,10)
  beta2.h ~ dnorm(0, var4)
  var4 ~ dunif(0,10)

  


  # functional response likelihood
  
for(i in 1:Nind){
log(a[i]) <- alpha.a[i] + beta1.a*log(mc[i]) + beta2.a*log(mr[i]) 
log(h[i]) <- alpha.h[i] + beta1.h*log(mc[i]) + beta2.h*log(mr[i]) 
}

for(i in 1:n){
killed[i] ~ dbin(prob[i],initial[i])
prob[i] <- max(0.0001,min(0.9999,T/(1/a[id[i]] + h[id[i]]*initial[i])))
}
  
  
}
", file = here("code/JAGs_models", "test2.txt"))

model.loc=here("code/JAGs_models","test2.txt")
jags.params=c("alpha.a", "beta1.a", "beta2.a", "alpha.h", "beta1.h", "beta2.h", "mu.alpha.h", "mu.alpha.a")

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
n.burnin = 100 * 1000
n.thin = 10
n.iter = 250 * 1000
model = R2jags::jags(jags.data,parameters.to.save=jags.params,inits=NULL,
                     model.file=model.loc, n.chains = n.chains, n.burnin=n.burnin,n.thin=n.thin, n.iter=n.iter, DIC=TRUE)













set.seed(1001)


jagsscript = cat("

                 model{
                 
                 # Piors
                 
                 
                 alpha.a ~ dnorm(exp(-4.5), 1)
                 beta1.a ~ dnorm(0, 1)
                 #var2 ~ dunif(0,1)
                 beta2.a ~ dnorm(0, 1)
                 #var3 ~ dunif(0,1)
                 

                 alpha.h ~ dnorm(exp(11), var4)
                 var4 ~ dunif(0,10)
                 beta1.h ~ dnorm(0, var5)
                 var5 ~ dunif(0,10)
                 beta2.h ~ dnorm(0, var6)
                 var6 ~ dunif(0,10)
                 
                 
                 # functional response likelihood
                 
                 for(i in 1:n){
                 
                 killed[i] ~ dbin(prob[i],initial[i])
                 
                 prob[i] <- max(0.0001,min(0.9999,T/(1/(alpha.a * mc[i]^beta1.a * mr[i]^beta2.a) + (alpha.h * mc[i]^beta1.h * mr[i]^beta2.h)*initial[i])))
                 
                 }
                 
                 
                 }
                 ", file = here("code/JAGs_models", "test.txt"))



model.loc=here("code/JAGs_models","test.txt")
jags.params=c("alpha.a", "beta1.a", "beta2.a", "alpha.h", "beta1.h", "beta2.h", "var4", "var5", "var6")

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


jagscript = cat("
model {
                for (i in 1:N) {
                ar[i] <- max(0.0001,avec[block[i]]*pow(size[i]/d,gamma)*exp(1-size[i]/d))
                prob[i] <- max(0.0001,min(0.9999,1/(1/ar[i]+h*N0[i])))
                killed[i] ~ dbin(prob[i],N0[i])	
                }
                for (i in 1:nblock) {
                avec0[i] ~ dnorm(0,tau.a)
                avec[i] <- a*exp(avec0[i])
                }
                ## Specify priors
                tau.a <- pow(sd.a,-2)
                sd.a ~ dunif(0,1)
                d ~ dlnorm(0,0.01)
                gamma ~ dlnorm(0,0.01)
                h ~ dlnorm(0,0.01)
                a ~ dlnorm(0,1)
                }", file = here("code", "PR_jags.txt"))


jagscript = cat("
model {
                for (i in 1:n) {
                #ar[i] <- avec[id[i]]*pow(mc[i],beta1.a)*pow(mr[i], beta2.a)
                ht[i] <- hvec[id[i]]*pow(mc[i],beta1.h)*pow(mr[i], beta2.h)
                prob[i] <- max(0.0001,min(0.9999,1/(1/a+ht[i]*initial[i])))
                killed[i] ~ dbin(prob[i],initial[i])	
                }

                for (i in 1:Nind) {
                # avec0[i] ~ dnorm(0,tau.a)
                # avec[i] <- a*exp(avec0[i])

                hvec0[i] ~ dnorm(0,tau.h)
                hvec[i] <- h*exp(hvec0[i])
                }

                ## Specify priors
                # tau.a <- pow(sd.a,-2)
                # sd.a ~ dunif(0,1)
                tau.h <- pow(sd.h,-2)
                sd.h ~ dunif(0,1)

                h ~ dlnorm(0,0.01)
                a ~ dlnorm(0,1)

                # Allometric scaling exponents for attack rate
                beta1.a ~ dnorm(0.75, var.beta1.a)
                var.beta1.a ~ dgamma(2,1)
                
                beta2.a ~ dnorm(0.5, var.beta2.a)
                var.beta2.a ~ dgamma(2,1)
                
                # Allometric scaling exponents for handling time
                beta1.h ~ dnorm(-0.75, var.beta1.h)
                var.beta1.h ~ dgamma(10,1)
                
                beta2.h ~ dnorm(0.5, var.beta2.h)
                var.beta2.h ~ dgamma(2,1)


                }", file = here::here("code/JAGs_models", "PR_jags.txt"))


model.loc=here::here("code/JAGs_models","PR_jags.txt")
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
                 "mc" = df$mc,
                 "mr" = df$mr,
                 "Nind" = length(unique(df$id)), 
                 "id" = df$id
) # named list

n.chains = 3
n.burnin = 100 * 1000 /4
n.thin = 10
n.iter = 250 * 1000 /4
model = R2jags::jags(jags.data,parameters.to.save=jags.params,inits=NULL,
                     model.file=model.loc, n.chains = n.chains, n.burnin=n.burnin,n.thin=n.thin, n.iter=n.iter, DIC=TRUE)


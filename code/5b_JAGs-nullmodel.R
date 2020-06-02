library(here)
source(here("code", "1_setup.R"))
source(here("code", "Base_code/functions.R"))


#####################################
## Get data
#####################################


df <- read.table(here("data/cleaned","loburc_cleaned.csv"), header = T, sep = ",") %>%
  arrange(id, treatment)

##################################################################
## the model
##################################################################

jagsscript = cat("
                 
                 model{
                 # PRIORS
                 
                 logit.a ~ dnorm(0, 0.01)
                 a <- exp(logit.a)/(1+exp(logit.a))
                 
                 log.h ~ dnorm(0, 0.01)
                 h <- exp(max(min(log.h, 20), -20))

                 # sigma_a ~ dunif(0,10)
                 # t.tau.a <- 1/(sigma_a*sigma_a)
                 # sigma_h ~ dunif(0,10)
                 # t.tau.h <- 1/(sigma_h*sigma_h)
                 
                 
                 # FUNCTIONAL RESPONSE LIKELIHOOD
                 
                 for(i in 1:n){
                 
                 prob[i] <- max(0.0001,min(0.9999,T/(1/a + h*initial[i])))
                 
                 killed[i] ~ dbin(prob[i],initial[i])
                 
                 }
                 
                 
                 }
                 ", file = here("code", "JAGs_models/Null_jags.txt"))



model.loc=here("code","JAGs_models/Null_jags.txt")
jags.params=c("a", "h", 'sigma_a', 'sigma_h'
)


jags.data = list("initial"= df$initial,
                 "killed" = df$killed,
                 "n" = length(df$initial), 
                 "T" = 48
) # named list

n.chains = 3
n.burnin = 10000
n.thin = 2
n.iter = 25000
model = R2jags::jags(jags.data,parameters.to.save=jags.params,inits=NULL,
                     model.file=model.loc, n.chains = n.chains, n.burnin=n.burnin,n.thin=n.thin, n.iter=n.iter, DIC=TRUE)


df.null <- as.mcmc(model) %>%
  recover_types(df) %>%
  spread_draws(a, h)

write.csv(df.null, here::here("data/cleaned/posteriors", "posteriors_null.csv"), row.names = F)







meta <- read.csv(here::here("data/", "lob-metadata.csv"))
df <- read.csv(here::here("data/cleaned", "posteriors_individuals.csv")) %>% 
  left_join(meta) %>%
  filter(id != "N07")


ggplot(df, aes(x = log(mc), y = log(h)))+
  geom_point()+
  facet_wrap(~treatment)


ltow <- function(udiam){0.000592598*udiam^2.872636198*1.01}
ltow(c(10,30,50,70))


#------------------------------------------------------------------
## Construct a model of the medians
#-----------------------------------------------------------------

d1 <- read.csv(here::here("data/cleaned", "posteriors_individuals.csv")) %>% 
  left_join(meta) %>%
  group_by(id, mr, mc) %>%
  filter(id != "N07") %>%
  summarise(mu.a = median(a), mu.h = median(h))


jagsscript = cat("
                 
                 model{
                 
                 for(i in 1:n){
                 
                 y[i] ~ dnorm(mu[i], tauy)
                 mu[i] <- alpha + beta1*mc[i] + beta2*mr[i]
                 
                 }
                 
                 #priors
                 alpha ~ dnorm(0, 0.01)
                 beta1 ~ dnorm(0, 0.01)
                 beta2 ~ dnorm(0, 0.01)
                 sigmay ~ dunif(0,50)
                 tauy <- 1/(sigmay*sigmay)
                 
                 }
                 ", file = here("code", "post-hoc_medians.txt"))



model.loc=here("code","post-hoc_medians.txt")
jags.params=c("alpha", "beta1", "beta2")


jags.data = list("y" = log(d1$mu.h),
                 "mc" = log(d1$mc),
                 "mr" = log(d1$mr),
                 "n" = length(d1$id)
) # named list

n.chains = 3
n.burnin = 10000
n.thin = 3
n.iter = 30000
medianh = R2jags::jags(jags.data,parameters.to.save=jags.params,inits=NULL,
                       model.file=model.loc, n.chains = n.chains, n.burnin=n.burnin,n.thin=n.thin, n.iter=n.iter, DIC=TRUE)

medianh


jags.data = list("y" = log(d1$mu.a),
                 "mc" = log(d1$mc),
                 "mr" = log(d1$mr),
                 "n" = length(d1$id)
) # named list

mediana = R2jags::jags(jags.data,parameters.to.save=jags.params,inits=NULL,
                       model.file=model.loc, n.chains = n.chains, n.burnin=n.burnin,n.thin=n.thin, n.iter=n.iter, DIC=TRUE)

mediana


#---------------------------------------------------------------------------
## Construct a model that incorporates uncertainty in the a and h estimates
#---------------------------------------------------------------------------

d2 <- read.csv(here::here("data/cleaned", "posteriors_individuals.csv")) %>% 
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
                 beta1 ~ dnorm(0, 0.01)
                 beta2 ~ dnorm(0, 0.01)
                 sigmay ~ dunif(0,10)
                 tauy <- 1/(sigmay*sigmay)
                 
                 }
                 ", file = here("code", "post-hoc_wparameteruncertainty.txt"))



model.loc=here("code","post-hoc_wparameteruncertainty.txt")
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

jags.data = list("y" = log(d2$a),
                 "mc" = log(d2$mc), 
                 "mr" = log(d2$mr),
                 "n" = length(d2$id)
) # named list

n.chains = 3
n.burnin = 10000
n.thin = 3
n.iter = 30000
model_wpuncert.a = R2jags::jags(jags.data,parameters.to.save=jags.params,inits=NULL,
                              model.file=model.loc, n.chains = n.chains, n.burnin=n.burnin,n.thin=n.thin, n.iter=n.iter, DIC=TRUE)



#---------------------------------------------------------------------------------------------
## Construct a model that incorporates uncertainty in the a and h estimates and urchin size
#---------------------------------------------------------------------------------------------


# use a normal prior on the mr's with a large variance and truncate based on differences from the mean

# plot up the urchin size posterior

precision <- function(sd){1/sd^2}

temp <- unique(d2$udiam)
prior <- list()
for(i in 1:3){
  t <- rnorm(1000, temp[i], 1/0.4^2)
  prior[[i]]  <- t[t > temp[i] - 10 & t < temp[i] +10]
}
names(prior) <- temp
out <- do.call("cbind", prior) # this is some funkiness to get a list of data.frames to a single data.frame

matplot(out, type = "p")


jagsscript = cat("
                 
                 model{
                 
                 for(i in 1:n){
                 
                 x[i] ~ dnorm(udiam[i], 6.25)T(udiam[i]-10, udiam[i]+10)
                 x.m[i] <- log(0.000592598*x[i]^2.872636198*1.01)
                 
                 y[i] ~ dnorm(mu[i], tauy)
                 mu[i] <- alpha + beta1*mc[i] + beta2*x.m[i]
                 
                 }
                 
                 #priors
                 alpha ~ dnorm(0, 0.01)
                 beta1 ~ dnorm(0, 0.01)
                 beta2 ~ dnorm(0, 0.01)
                 sigmay ~ dunif(0,10)
                 tauy <- 1/(sigmay*sigmay)
                 
                 }
                 ", file = here("code", "post-hoc_wmeasurementerror.txt"))



model.loc=here("code","post-hoc_wmeasurementerror.txt")
jags.params=c("alpha", "beta1", "beta2")


jags.data = list("y" = log(d2$h),
                 "mc" = log(d2$mc),
                 "n" = length(d2$id), 
                 "udiam" = d2$udiam
) # named list

n.chains = 3
n.burnin = 10000
n.thin = 3
n.iter = 30000
model_wmeasurement = R2jags::jags(jags.data,parameters.to.save=jags.params,inits=NULL,
                                  model.file=model.loc, n.chains = n.chains, n.burnin=n.burnin,n.thin=n.thin, n.iter=n.iter, DIC=TRUE)



jags.data = list("y" = log(d2$a),
                 "mc" = log(d2$mc),
                 "n" = length(d2$id), 
                 "udiam" = d2$udiam
) # named list

n.chains = 3
n.burnin = 10000
n.thin = 3
n.iter = 30000
model_wmeasurement.a = R2jags::jags(jags.data,parameters.to.save=jags.params,inits=NULL,
                                  model.file=model.loc, n.chains = n.chains, n.burnin=n.burnin,
                                  n.thin=n.thin, n.iter=n.iter, DIC=TRUE)














df.amedian <- as.mcmc(mediana) %>%
  spread_draws(alpha, beta1, beta2) %>%
  mutate(model = "model_median")

df.awpuncert <- as.mcmc(model_wpuncert.a) %>%
  spread_draws(alpha, beta1, beta2) %>%
  mutate(model = "model_wpuncert")

df.awmeasurement <- as.mcmc(model_wmeasurement.a) %>%
  spread_draws(alpha, beta1, beta2) %>%
  mutate(model = "model_wmeasurement")

df.a <- bind_rows(df.amedian, df.awpuncert, df.awmeasurement) %>%
  gather(parameter, estimate, -c(model, .chain, .iteration, .draw))

write.csv(df.a, here::here("data/cleaned/", "posteriors_posthoc_a.csv"), row.names = F)


df.hmedian <- as.mcmc(medianh) %>%
  spread_draws(alpha, beta1, beta2) %>%
  mutate(model = "model_median")

df.hwpuncert <- as.mcmc(model_wpuncert) %>%
  spread_draws(alpha, beta1, beta2) %>%
  mutate(model = "model_wpuncert")

df.hwmeasurement <- as.mcmc(model_wmeasurement) %>%
  spread_draws(alpha, beta1, beta2) %>%
  mutate(model = "model_wmeasurement")

df.h <- bind_rows(df.hmedian, df.hwpuncert, df.hwmeasurement) %>%
  gather(parameter, estimate, -c(model, .chain, .iteration, .draw))

write.csv(df.h, here::here("data/cleaned/", "posteriors_posthoc_h.csv"), row.names = F)















#-----------------------------------
## Scrap code
#-----------------------------------

df.h %>%
  ggplot(aes(x = estimate))+
  geom_density(aes(fill = model), alpha = 0.5)+
  facet_wrap(~parameter, scales = "free")


ggplot(d2, aes(x = mc, y = h))+
  stat_summary(fun = "median", geom = "point", col = "red")+
  facet_wrap(~treatment)


calculate.CI <- function(a.parameter, h.parameter, prob = 0.95){
  require(rethinking)
  
  mu.link <- function(N, a, h){
    df.model[,a.parameter]*N/(1+df.model[,a.parameter]*df.model[,h.parameter]*N)
  } # defines a function to predict the prey killed at combination of a and h in the posteriors
  
  N.seq <- seq( from=0 , to=60 , length.out = 100 ) # define a sequence of initial densities
  mu <- sapply( N.seq , mu.link) # apply the mu.link funciton to each N in the sequence
  
  mu.median <- apply( mu , 2 , median ) # calculate the median predicted value for each N
  mu.PI <- t(apply( mu , 2 , PI , prob=prob )) # calculate the credible interval for each value of N
  
  return(data.frame(N.seq = N.seq, mu = mu.median, mu.lower = mu.PI[,1], mu.upper = mu.PI[,2])) # return a data.frame to organize output!
  
}

library(here)
source(here("code", "1_setup.R"))
source(here("code", "Base_functions/functions.R"))


#####################################
## Get data
#####################################


df <- read.table(here("data/cleaned","loburc_cleaned.csv"), header = T, sep = ",") %>%
  arrange(id, treatment)


dunn.h <- 0.741 * 24 #days
log(dunn.h)

# tank area in Dunn 
area.dunn <- pi*(2.18/2)^2 + 0.91*2.18
dunn.a <- 0.194 / area.dunn / 
  24 / #convert to hours
  0.5 #convert to units of our tanks ~ 2 m2.
logit(dunn.a)

##################################################################
## the model
##################################################################

jagsscript = cat("
                 
 model{
 # PRIORS
 
 # hyperprior
 
 mu.logit.a ~ dnorm(-5.865, mu.tau.a) # Dunn prior on attack rate
 mu.a <- exp(mu.logit.a)/(1+exp(mu.logit.a))
 
 mu.log.h ~ dnorm(2.878, mu.tau.h) # Dunn prior on h
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

# Generate posterior predictive distributions

#Population level
post.logit.mu.a ~ dnorm(mu.logit.a, mu.tau.a)
post.mu.a <- exp(post.logit.mu.a)/(1+exp(post.logit.mu.a))

post.log.mu.h ~ dnorm(mu.log.h, mu.tau.h)
post.mu.h <- exp(post.log.mu.h)

for(i in 1:n){
prob.mu[i] <- max(0.0001,min(0.9999,T/(1/post.mu.a + post.mu.h*initial[i])))
post.mu[i] ~ dbin(prob.mu[i],initial[i])
e.mu[i] <- killed[i] - post.mu[i]
}


#Treatment level
for(i in 1:Ntreats){
 #a
post.t.logit.a[i] ~ dnorm(t.logit.a[i], t.tau.a)
post.t.a[i] <- exp(post.t.logit.a[i])/(1+exp(post.t.logit.a[i]))

#h
post.t.log.h[i] ~ dnorm(t.log.h[i], t.tau.h)
post.t.h[i] <- exp(max(min(post.t.log.h[i],20),-20))
}

for(i in 1:n){
prob.t[i] <- max(0.0001,min(0.9999,T/(1/post.t.a[treat[i]] + post.t.h[treat[i]]*initial[i])))
post.t[i] ~ dbin(prob.t[i],initial[i])
e.t[i] <- killed[i] - post.t[i]
e.mu.t[i] <- post.t[i] - post.mu[i]
}



#Individual level
for(i in 1:num.ind){
#a
post.ind.logit.a[i] ~ dnorm(loga[i], tau_int.a)
post.ind.a[i] <- exp(post.ind.logit.a[i])/(1+exp(post.ind.logit.a[i]))

#h
post.ind.log.h[i] ~ dnorm(logh[i], tau_int.h)
post.ind.h[i] <- exp(max(min(post.ind.log.h[i],20),-20))
}

for(i in 1:n){
prob.ind[i] <- max(0.0001,min(0.9999,T/(1/post.ind.a[id[i]] + post.ind.h[id[i]]*initial[i])))
post.ind[i] ~ dbin(prob.ind[i],initial[i])
e.ind[i] <- killed[i] - post.ind[i]
e.ind.t[i] <- post.ind[i] - post.t[i]
e.ind.mu[i] <- post.ind[i] - post.mu[i]
}




}
                 ", file = here("code", "JAGs_models/heirarchical_jags.txt"))



model.loc=here("code","JAGs_models/heirarchical_jags.txt")
jags.params=c("a", "h", "t.a", "t.h", "mu.a", "mu.h", 
              "sigma_int.a", "sigma_int.h", "sigma_t.a", "sigma_t.h", 
              "loga", "logh", "t.logit.a", "t.log.h", "sigma_mu.a", "sigma_mu.h", "mu.tau.a", "mu.tau.h", "post.mu", "post.t", "post.ind", "e.mu", "e.t", "e.mu.t", "e.ind", "e.ind.t", "e.ind.mu"
)


tind <- distinct(df, treatment, id) %>%
  arrange(id)

tind <- as.factor(as.vector(tind$treatment))

id.vec <- distinct(df, id) %>%
  arrange(id)
id.vec <- as.factor(as.vector(id.vec$id))

jags.data = list("initial"= df$initial,
                 "killed" = df$killed,
                 "id" = df$id,
                 "num.ind" = length(unique(df$id)),
                 "n" = length(df$initial), 
                 "tind" = tind, 
                 "T" = 48, 
                 "Ntreats" = length(unique(df$treatment)), 
                 "treat" = df$treatment
) # named list

n.chains = 3
n.burnin = 250000/10
n.thin = 10
n.iter = 500000/10
model = R2jags::jags(jags.data,parameters.to.save=jags.params,inits=NULL,
                     model.file=model.loc, n.chains = n.chains, n.burnin=n.burnin,n.thin=n.thin, n.iter=n.iter, DIC=TRUE)


# pp <- as.mcmc(model) %>%
#   spread_draws(post.mu[junk] | junk) %>%
#   sample_draws(1000) %>%
#   gather_draws()
# 
# out <- data.frame(killed = y, initial = df$initial, t(apply(post.mu, 2, PI)))
# 
# ggplot(out, aes(x = initial, killed))+
#   geom_jitter()+
#   geom_ribbon(aes(ymin = X5., ymax = X94.))
# 
# out <- data.frame(killed = y, initial = df$initial, id = df$id, t(apply(post.ind, 2, PI, prob = 0.75)))
# ggplot(out, aes(x = initial, killed))+
#   geom_jitter()+
#   geom_ribbon(aes(ymin = X12., ymax = X88.))+
#   facet_wrap(~id)


y = df$killed
trials = df$initial
# The following code was based on the recommendations of Gelman and Pardoe 2006. It seems to work for the level of the data, see print(c(rsquared.y_mu, rsquared.y_t, rsquared.y_ind)), but not for comparisons between heiarcical levels. 


post.mu <- MCMCchains(model, params = "post.mu")
e.y_mu <- MCMCchains(model, params = "e.mu")
rsquared.y_mu <- 1 - mean (apply (e.y_mu, 1, var)) / var (y)
#rsquared.y_mu <- mean (apply (e.y_mu, 1, var)) / (var (y) + mean(apply(e.y_mu, 1, var)))

post.t <- MCMCchains(model, params = "post.t")
e.y_t <- MCMCchains(model, params = "e.t")
rsquared.y_t <- 1 - mean (apply (e.y_t, 2, var)) / var (y)
e.mu_t <- MCMCchains(model, params = "e.mu.t")
rsquared.mu_t <- 1 - mean (apply (e.mu_t, 2, var)) / mean(apply(post.t, 2, var))

post.ind <- MCMCchains(model, params = "post.ind")
e.y_ind <- MCMCchains(model, params = "e.ind")
rsquared.y_ind <- 1 - mean (apply (e.y_ind, 2, var)) / var (y)
#rsquared.y_ind <- mean (apply (e.y_ind, 2, var)) / (var(y) + mean (apply (e.y_ind, 2, var)))

e.t_ind <- MCMCchains(model, params = "e.ind.t")
rsquared.t_ind <- 1 - mean (apply (e.t_ind, 2, var)) / mean(apply(post.ind, 2, var))
e.mu_ind <- MCMCchains(model, params = "e.ind.mu")
rsquared.mu_ind <- 1 - mean (apply (e.mu_ind, 2, var)) / mean(apply(post.ind, 2, var))

print(c(rsquared.y_mu, rsquared.y_t, rsquared.y_ind))
print(c(rsquared.mu_ind, rsquared.t_ind))


# Technique based on Gelman et al. 2019
    # as of 8/21/20, I use these estimates in the results, but I'm worried that they are not correct! See bayes_R2.R in code folder for original code that I adapted.

bayes_R2_res <- function(y, trials, ypred) {
  #ypred <- ypred %*% diag(trials) # I did not include this because I don't understand why this step is necessary. The ides is to estimate the residual variance / residual variance + residual error variance. So why multiple the vectors of ypred by the diag of the trials? 
  if(is.matrix(y) == TRUE){
    e <- y - ypred
  }else{
  e <- -1 * sweep(ypred, 2, y)
  }
  var_ypred <- apply(ypred, 1, var)
  var_e <- apply(e, 1, var)
  var_ypred / (var_ypred + var_e)
}

#variance of data explained by population level predictions
pop.mean <- mean(bayes_R2_res(y = y, trials = trials, ypred = post.mu ))
hist(bayes_R2_res(y = y, trials = trials, ypred = post.mu ))

#variance of data explained by treatment level predictions
t.mean <- mean(bayes_R2_res(y = y, trials = trials, ypred = post.t ))
hist(bayes_R2_res(y = y, trials = trials, ypred = post.t ))

#variance of data explained by individual level prediction
ind.mean <- mean(bayes_R2_res(y = y, trials = trials, ypred = post.ind ))
hist(bayes_R2_res(y = y, trials = trials, ypred = post.ind ))

#variance of the population level estimate explained by the individual level prediction
mean(bayes_R2_res(y = post.mu, trials = trials, ypred = post.ind))
hist(bayes_R2_res(y = post.mu, trials = trials, ypred = post.ind))

#variance of the treatment level estimates explained by the individual level prediction
mean(bayes_R2_res(y = post.t, trials = trials, ypred = post.ind))
hist(bayes_R2_res(y = post.t, trials = trials, ypred = post.ind))

#variance of the population level estimates explained by the treatment level prediction
mean(bayes_R2_res(y = post.mu, trials = trials, ypred = post.t))
hist(bayes_R2_res(y = post.mu, trials = trials, ypred = post.t))


ypred <- post.mu %*% diag(trials)
e <- -1 * sweep(ypred, 2, y)
var_ypred <- apply(ypred, 1, var)
var_e <- apply(e, 1, var)
mean(var_ypred / (var_ypred + var_e))




test <- c(10,21,6,2)

y %*% diag(test)
y <- y[, 1]

y <- cbind(c(9, 16, 2, 0), c(1,5,4,2))
ypred <- matrix(nrow = 4, c(8, 17, 1, 1, 
                9, 16, 2, 1, 
                6, 15, 1, 1, 
                7, 20, 0, 0), byrow = T)

ypred2 <- matrix(nrow = 4, rbinom(16, size = 21, prob = 0.5), byrow = T)

ypred - ypred2

e <- -1 * sweep(ypred, 2, ypred2)

trials <- rowSums(y)
y <- y[, 1]
ypred <- ypred %*% diag(trials)

e <- -1 * sweep(ypred, 2, y)
var_ypred <- apply(ypred, 1, var)
var_e <- apply(e, 1, var)
var_ypred / (var_ypred + var_e)








e.y_mu <- -1*sweep(post.mu, 2, y)
e.mu_t <- -1*sweep(post.t, 2, cbind(post.mu, post.mu, post.mu))
e.y_ind <- -1*sweep(post.ind, 2, rep(y, times = jags.data$num.ind))
e.t_ind <- -1*sweep(post.ind, 2, post.t)
e.mu_ind <- -1*sweep(post.ind, 2, post.mu)

test <- matrix(nrow = 2, ncol = 2, data = rbinom(4, 10, 0.5))
test2 <- matrix(nrow = 2, ncol = 2*3, data = rbinom(12, 10, 0.5))

sweep(test2, 2, cbind(test, test, test))

rep(test, times = 3)



# test <- as.mcmc(model) %>%
#   recover_types(df) %>%
#   spread_draws(post.mu)
# 
# df.treat <- as.mcmc(model) %>%
#   recover_types(df) %>%
#   spread_draws(t.a[treatment], t.h[treatment])
# 
# df.ind <- as.mcmc(model) %>%
#   recover_types(df) %>%
#   spread_draws(a[id], h[id])
# 
# df.pop <- as.mcmc(model) %>%
#   recover_types(df) %>%
#   spread_draws(mu.a, mu.h, mu.tau.a, mu.tau.h, sigma_int.a, sigma_int.h, sigma_t.a, sigma_t.h)
# 
# write.csv(df.pop, here::here("data/cleaned/posteriors", "posteriors_population.csv"), row.names = F)
# write.csv(df.treat, here::here("data/cleaned/posteriors", "posteriors_treatments.csv"), row.names = F)
# write.csv(df.ind, here::here("data/cleaned/posteriors", "posteriors_individuals.csv"), row.names = F)


















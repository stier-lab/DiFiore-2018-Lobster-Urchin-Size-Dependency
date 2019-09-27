

###############################################
## Example
###############################################

x <- 1:1000

y <- 2 * rnorm(1000, x, 10)

plot(y ~ x)


lm <- lm(y ~ x)
summary(lm)


jagsscript = cat("
model{
  
  # likelihood function
for(i in n){
  y[i] ~ dnorm(mu[i], tau)
  mu[i] <- b*x[i]
}
  # priors
  b ~ dnorm(0, 0.001)
  
  #Variances
  
  sigma ~ dunif(0,10)
  tau <- 1/(sigma*sigma)
  
}
", file = here("code", "lmexample.txt"))




model.loc=here("code", "lmexample.txt")
jags.params=c("b", "sigma")


# Test against simulated dataset.
jags.data = list("x" = x, 
                 "y" = y, 
                 "n" = length(x)
) # named list

n.chains = 3
n.burnin = 5000
n.thin = 2
n.iter = 10000
model = R2jags::jags(jags.data, parameters.to.save=jags.params,inits=NULL,
                     model.file=model.loc, n.chains = n.chains, n.burnin=n.burnin,
                     n.thin=n.thin, n.iter=n.iter, DIC=TRUE)
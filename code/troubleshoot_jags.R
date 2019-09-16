############################################################
## Setup
############################################################

library(tidyverse)
library(ggplot2)
library(tidyr)
library(dplyr)
library(R2jags)
library(rjags)
library(MCMCvis)
library(ggpubr)
library(cowplot)
library(here)
library(ggpubr)
library(bbmle)
library(lme4)



#############################################################
## Simulate data
#############################################################


# write the holling II funcitonal response in terms of the probability a prey is killed
mortrisk <- function(N0,h,a){
  risk <- 1/(1/a + h*N0)                ## Holling risk (per capita)
  pmin(1.0,risk)                         ## bound risk <= 1
}

holling2=function(N,a,h) {
  a*N/(1+a*h*N)
}



set.seed(18) ## set random-number seed for reproducibility
## The data will simulate a case where initial attack rate and handling time vary between individual lobsters

simdata <- function(nind = 10){
  
  N0=c(2,3,5,10,16,26)
  test.vals <- expand.grid(N0 = N0,
                           ind=1:nind)
  ## attack rate varies randomly by individual with a median of 0.5
  ##  and proportional variation of approx 10%
  a <- rlnorm(nind,meanlog=log(0.5),sdlog=0.1)
  ## handling time varies randomly by individual with a median of 0.1
  ##  and proportional variation of approx 10%
  h <- rlnorm(nind,meanlog=log(0.1),sdlog=0.1)
  p <- with(test.vals,mortrisk(N0=N0,
                               a=a[ind],
                               h=h[ind]))
  z <- rbinom(nrow(test.vals),prob=p,size=test.vals$N0)
  data.frame(test.vals,killed=z, a = rep(a, each = length(N0)), h = rep(h, each = length(N0)))
}
x <- simdata(10)[,1:3]

# plot it up

names(x) <- c("initial", "ind", "killed")
library("lattice")

xyplot(killed~initial|as.factor(ind), data = x)

################################################################
## Fit with ML no random effects
################################################################

## Now fit the data using maximum likelihood without individual effects
eps <- 1e-4 ## used to bound probabilities between 0 and 1
##Using the classical unstructured Type II functional response
classic <- bbmle::mle2(killed~dbinom(prob=pmax(eps,
                                        pmin(1-eps,1/(1/a+h*initial))),
                              size=initial),start=list(a=.01,h=.001),data=x)

# Ok, so we are able to recover the attach rate, and handling time that generated the data, using a classic ML approach with no effect of individual


#################################################################
## Heirarchical bayesian model -- individual fits
#################################################################

# in this model I am using an informative prior based on the a recent study of the functional response between lobsters and urchins (Dunn et al. 2018, Ecology). Basically, I used the bootstrapped estimates of the functional repsonse fits from this paper to estimate the shape and scale parameters of a gamma distrubtion.

shape.a <- 10.25596
scale.a <- 0.02022754

shape.h <- 11.3787
scale.h <- 0.06735276


jagsscript = cat("model{
                 
                 # hyperpriors
                 #a
                 mu.a ~ dgamma(shape.a, 1/scale.a) #mean hyperparameter for random a
                 sigma_int.a ~dunif(0,10) #SD hyperparatmer for random a
                 tau_int.a <- 1/(sigma_int.a*sigma_int.a)
                 
                 #h
                 mu.h ~ dgamma(shape.h, 1/scale.h) #mean hyperparameter for random h
                 sigma_int.h ~dunif(0,10) #SD hyperparatmer for random h
                 tau_int.h <- 1/(sigma_int.h*sigma_int.h)
                 
                 # priors
                 #Individual attack rates and handling times vary according to a log normal distribution. 
                 
                 for(i in 1:num.ind){
                 a[i] ~ dlnorm(mu.a, tau_int.a)
                 h[i] ~ dlnorm(mu.h, tau_int.h)
                 }
                 
                 
                 #likelihood
                 
                 for(i in 1:n){
                 
                 prob[i] <- max(0.0001,min(0.9999,1/(1/a[id[i]] + h[id[i]]*initial[i])))
                 
                 killed[i] ~ dbin(prob[i],initial[i])
                 
                 }
                 
                 }", file = here("code", "heirarchical_jags3.txt"))

model.loc=here("code","heirarchical_jags3.txt")
jags.params=c("a", "h", "mu.a", "mu.h", "sigma_int.a", "sigma_int.h")


# Run the model against the simulated dataset
jags.data = list("initial"= x$initial,
                 "killed" = x$killed,
                 "id" = x$ind,
                 "num.ind" = length(unique(x$ind)),
                 "n" = length(x$initial), 
                 "scale.a" = scale.a,
                 "shape.a" = shape.a,
                 "scale.h" = scale.h,
                 "shape.h" = shape.h
) # named list

n.chains = 3
n.burnin = 10000
n.thin = 2
n.iter = 50000
model = R2jags::jags(jags.data,parameters.to.save=jags.params,inits=NULL,
                     model.file=model.loc, n.chains = n.chains, n.burnin=n.burnin,
                     n.thin=n.thin, n.iter=n.iter, DIC=TRUE)

print(model, dij = 3)

a = MCMCsummary(model,params='a')
h = MCMCsummary(model,params='h')

mu.a = MCMCsummary(model,params='mu.a')
mu.h = MCMCsummary(model,params='mu.h')


# Plot individual level fits
d <- par(mfrow = c(3,4), mar = c(4,4,2,1))
for(i in 1:10){
  plot(killed ~ jitter(initial),data = x[as.numeric(x$ind) == i, ], xlab="Number of prey",ylab="Number consumed", ylim = c(0,12))
  curve(holling2(x,a[i,1],h[i,1]),add=TRUE,col=1,lty=1) #true curve
}
par(d)

# Plot population level fits and compare with classic ML approach
plot(killed ~ jitter(initial),data = x, xlab="Number of prey",ylab="Number consumed", ylim = c(0,12), 
     main = "Polulation level estimates of the functional response \nfit to simulated data")
        #classic curve
        curve(holling2(x,classic@coef[1], classic@coef[2]),add=TRUE, col = "red") 
        
        #jags random effects curve
        curve(holling2(x,mu.a[,1],mu.h[,1]),add=TRUE,col= "darkgreen")
        
        # true curve
        curve(holling2(x,0.5,0.1),add=TRUE, lty = 4)
        
        legend("topleft", legend = c("ML (no random effects) Model", "JAGS random effects model", "True curve"), col = c("red", "darkgreen", "black"), lty = c(1,1,4))








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
mortrisk <- function(N0,h,a,T){
  risk <- T/(1/a + h*N0)                ## Holling risk (per capita)
  pmin(1.0,risk)                         ## bound risk <= 1
}

holling2=function(N,a,h,T) {
  a*T*N/(1+a*h*N)
}


set.seed(18) ## set random-number seed for reproducibility
## The data will simulate a case where initial attack rate and handling time vary between individual lobsters

simdata <- function(nind = 30, ngroups = 3, T = 1){
  
  N0=c(2,3,5,10,16,26)
  test.vals <- expand.grid(N0 = N0,
                           ind=1:nind)
  test.vals$group <- rep(1:ngroups, each = length(N0)*(nind/ngroups))
  a <- list()
  h <- list()
  
for(i in 1:ngroups){
  group.a <- c(0.8, 0.5, 0.1)
  group.h <- c(1/26, 1/10, 1/4)
  ## attack rate varies randomly by individual with a median of 0.5
  ##  and proportional variation of approx 10%
  a[[i]] <- rlnorm(nind/ngroups,meanlog=log(group.a[i]),sdlog=0.01)
  ## handling time varies randomly by individual with a median of 0.1
  ##  and proportional variation of approx 10%
  h[[i]] <- rlnorm(nind/ngroups,meanlog=log(group.h[i]),sdlog=0.01)
}

  ar <- unlist(a)
  ht <- unlist(h)
  
  p <- with(test.vals,mortrisk(N0=N0,
                               a=ar[ind],
                               h=ht[ind], 
                               T = T))
  z <- rbinom(nrow(test.vals),prob=p,size=test.vals$N0)
  data.frame(test.vals,killed=z, a = rep(ar, each = length(N0)), h = rep(ht,  each = length(N0)))
}
x <- simdata(nind = 12, ngroups = 3)

# plot it up

names(x) <- c("initial", "ind", "group", "killed", "a", "h")
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

summary(classic)



#################################################################
## Heirarchical bayesian model -- individual fits
#################################################################

jagsscript = cat("

model{


# PRIORS

# Treatment level effects
      for(i in 1:Ntreats){
        #a
        t.logit.a[i] ~ dnorm(0, t.tau.a)
        t.a[i] <- exp(t.logit.a[i])/(1+exp(t.logit.a[i]))
        
        #h
        t.h[i] ~ dlnorm(0, t.tau.h)
      
      }


# Individual level effects
      for(i in 1:num.ind){
        
        # a
        loga[i] ~ dnorm(t.a[tind[i]], tau_int.a)
        a[i] <- exp(loga[i])/(1+exp(loga[i]))
        
        # h
        logh[i] ~ dnorm(log(t.h[tind[i]]), tau_int.h)
        h[i] <- exp(max(min(logh[i],20),-20))
      }

# functional response likelihood

#For good estimates at the individual level
      for(i in 1:n){
        
        prob[i] <- max(0.0001,min(0.9999,T/(1/a[id[i]] + h[id[i]]*initial[i])))
        
        killed[i] ~ dbin(prob[i],initial[i])
      
      }

# # For good estimates at the treatment level
#       for(i in 1:n){
#         
#         prob[i] <- max(0.0001,min(0.9999,T/(1/t.a[t[i]] + t.h[t[i]]*initial[i])))
#         
#         killed[i] ~ dbin(prob[i],initial[i])
#         
#       }


# Variances for all levels
    sigma_int.a ~dunif(0,10)
    tau_int.a <- 1/(sigma_int.a*sigma_int.a)
    sigma_int.h ~dunif(0,10)
    tau_int.h <- 1/(sigma_int.h*sigma_int.h)

    sigma_t.a ~ dunif(0,10)
    t.tau.a <- 1/(sigma_t.a*sigma_t.a)
    sigma_t.h ~ dunif(0,10)
    t.tau.h <- 1/(sigma_t.h*sigma_t.h)

}

", file = "heirarchical_troubleshoot.txt")

model.loc= "heirarchical_troubleshoot.txt"
jags.params=c("a", "h", "t.a", "t.h")


# Run the model against the simulated dataset
jags.data = list("initial"= x$initial,
                 "killed" = x$killed,
                 "id" = as.factor(x$ind),
                 "num.ind" = length(unique(x$ind)),
                 "n" = length(x$initial), 
                 "Ntreats" = length(unique(x$group)),
                 "T" = 1, 
                 "t" = x$group,
                 "tind" = as.factor(as.vector(distinct(x, group, ind)[,2]))
) # named list

n.chains = 4
n.burnin = 25000
n.thin = 2
n.iter = 50000
model = R2jags::jags(jags.data,parameters.to.save=jags.params,inits=NULL,
                     model.file=model.loc, n.chains = n.chains, n.burnin=n.burnin,
                     n.thin=n.thin, n.iter=n.iter, DIC=TRUE)

MCMCtrace(model)

print(model, dij = 3)

a = MCMCsummary(model,params=c('a'), round = 4)
h = MCMCsummary(model,params=c('h'), round = 4)

t.a = MCMCsummary(model,params=c('t.a'), round = 4)
t.h = MCMCsummary(model,params=c('t.h'), round = 4)

# Plot individual level fits
d <- par(mfrow = c(3,4), mar = c(4,4,2,1))
for(i in 1:12){
  plot(killed ~ jitter(initial),data = x[as.numeric(x$ind) == i, ], xlab="Number of prey",ylab="Number consumed", ylim = c(0,26))
  curve(holling2(x,a[i,1],h[i,1], T=1),add=TRUE,col=1,lty=1) #true curve
}
par(d)

# Plot treatment level fits
d <- par(mfrow = c(1,3), mar = c(4,4,2,1))
for(i in 1:3){
  plot(killed ~ jitter(initial),data = x[as.numeric(x$group) == i, ], xlab="Number of prey",ylab="Number consumed", ylim = c(0,26))
  curve(holling2(x,t.a[i,1],t.h[i,1], T=1),add=TRUE,col=1,lty=1) #true curve
}
par(d)







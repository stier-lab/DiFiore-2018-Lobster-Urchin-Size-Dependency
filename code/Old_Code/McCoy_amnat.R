## Simulate data for prey mortality risk assuming a predator with a Holling type II Functional response
## and for which prey risk is power-Ricker function of prey size. Handling time is assumed to be independent of prey size.

## Load needed libraries
library(bbmle)
library(coda)

use_R2WinBUGS <- FALSE  		##Specify whether you will use WinBUGS or JAGs
use_R2jags <- TRUE

## Predation risk function: power-Ricker as a function of size
mortrisk <- function(N0,size,h,a,d,gamma){
  ar <- a*(size/d)^gamma*(exp(1-size/d)) ## power-Ricker attack rate
  risk <- ar/(1+ar*N0*h)                 ## Holling risk (per capita)
  pmin(1.0,risk)                         ## bound risk <= 1
}

## Generate data
set.seed(1001) ## set random-number seed for reproducibility

## The data will simulate a case where initial attack rate "a" varies
## across 6 replicate blocks
simdata <- function(nblock){
  test.vals <- expand.grid(N0=seq(10,100,by=10),
                           size=seq(5,50,by=5),
                           block=1:nblock)
  ## attack rate varies randomly by block with a median of 0.5
  ##  and proportional variation of approx 10%
  a <- rlnorm(nblock,meanlog=log(0.5),sdlog=.10)
  p <- with(test.vals,mortrisk(N0=N0,
                               size=size,
                               a=a[block],
                               d=6,h=0.01,
                               gamma=2))
  z <- rbinom(nrow(test.vals),prob=p,size=test.vals$N0)
  data.frame(test.vals,killed=z)
}
x <- simdata(6)

## Plot results ...
with(x,plot(jitter(N0),killed,col=size,pch=block))

## Now fit the data using maximum likelihood without block effects
eps <- 1e-4 ## used to bound probabilities between 0 and 1

##Using the classical unstructured Type II functional response
classic <- mle2(killed~dbinom(prob=pmax(eps,
                                        pmin(1-eps,1/(1/a+h*N0))),
                              size=N0),start=list(a=.01,h=.001),data=x)

##Using a size dependent functional response based on the power-Ricker
powerRicker <- mle2(killed~dbinom(prob=pmax(eps,
                                            pmin(1-eps,
                                                 1/(1/(a*(size/d)^gamma*(exp(1-size/d)))+h*N0))),
                                  size=N0),
                    start=list(a=.3,h=.001,d=5,gamma=2), data=x)

## profile(powerRicker)  
## Compare model fits using AIC
AICtab(classic,powerRicker, weights=TRUE, sort=TRUE)

predframe <- data.frame(N0=20,size=5:50)
plot(predframe$size,predict(powerRicker,newdata=predframe),type="l")
lines(predframe$size,predict(classic,newdata=predframe),col=2)
with(subset(x,N0==20),
     points(killed~size))

## Now estimate the parameters of the size dependent functional response using MCMC approach incorporating random effects
## Extract starting values based on likelihood fit above and add a starting value for the random effect

xstart <- c(as.list(coef(powerRicker)),list(sd.a=0.1))
#x$block=factor(x$block)
##Create a list of starting values that will be used for each chain based on xstart above
xstart2 = list(
  as.list(unlist(xstart)*runif(length(xstart),0.009,.011)),
  as.list(unlist(xstart)*runif(length(xstart),0.09,.11)),
  as.list(unlist(xstart)*runif(length(xstart),0.9,1.1)))

## Now create a data list for BUGS/JAGS
dat <- c(as.list(x), N=nrow(x), nblock=6)

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
    
b0 <- jags(data=dat, #this line calls jags and records computation time
               
               inits=c(list(xstart),xstart2), #this line initializes the starting values
               n.chains=4, #specify the number of chains to run
               n.iter=10000, #specify the number of iterations
               parameters=c(names(xstart),"avec"), #specify which output you want to extract
               ## for R2jags, bugs_model code above must be stored in an external file -- named "PR_jags.txt" in this example
               model.file=here("code", "PR_jags.txt"),
               working.directory=getwd())
b <- b0$BUGSoutput
  
  print(b,digits=3)
  ## View diagnostic plots etc. to check convergence
  plot(b) ##View plot of parameter estimates and errors
  b2 <- as.mcmc.list(b)
  devpos <- which(colnames(b2[[1]])=="deviance")
  b3 <- b2[,-devpos] ## drop deviance estimates
  class(b3) <- "mcmc.list"
  xyplot(b3,layout=c(3,4),asp="fill") ##Visual inspection of Markov chain behavior
  densityplot(b3,layout=c(3,4),asp="fill") ##Density plot for parameter estimates
  gelman.diag(b2) #provides the Gelman-Rubin statistics
  
  
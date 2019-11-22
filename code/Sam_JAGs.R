###############################################################
## test on Sam's data
###############################################################
library(here)
source(here("code","functions.R"))
source(here("code","setup.R"))

#uncode the following lines if you don't run through the repo!

# ##########################################
# ## Libraries
# ##########################################
# library(tidyverse)
# library(ggplot2)
# library(tidyr)
# library(dplyr)
# library(R2jags)
# library(rjags)
# library(MCMCvis)
# library(ggpubr)
# library(cowplot)
# library(here)
# library(ggpubr)
# library(bbmle)
# library(lme4)
# library(lmerTest)
# 
# 
# holling2=function(N,a,h,P,T) {
#   a*N*P*T/(1+a*h*N)
# }

sam <- read.csv(here("data/samdata", "fr_data.csv"))
sam$id <- as.numeric(sam$lobster_id)

sam$temp2 <- as.factor(paste("t", sam$temp, sep = ""))

s <- sam[,c( "temp2", "lobster_id","Initial", "Killed")]
names(s) <- c("temp", "id", "initial", "killed")

s <- arrange(s, temp, id, initial)

lattice:: xyplot(killed~initial|as.factor(id), data = s)


jagsscript = cat("

model{


# PRIORS

# Population level estimates
  
       mu.logit.a ~ dnorm(0, 0.01)
       mu.a <- exp(mu.logit.a)/(1+exp(mu.logit.a))
       mu.h ~ dnorm(0, 0.01)T(0,)
                 

# Treatment level effects
    for(i in 1:Ntreats){
        #a
        t.logit.a[i] ~ dnorm(mu.logit.a, t.tau.a)
        t.a[i] <- exp(t.logit.a[i])/(1+exp(t.logit.a[i]))
        
        #h
        t.log.h[i] ~ dnorm(mu.h, t.tau.h)
        t.h[i] <- exp(max(min(t.log.h[i],20),-20))
    
    }


# Individual level effects
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
                 
                 ", file = here("code", "heirarchical_jagsSAM.txt"))


model.loc=here("code", "heirarchical_jagsSAM.txt")
jags.params=c("a", "h", "t.a", "t.h", "mu.a", "mu.h")


tind <- distinct(s, temp, id) %>%
  arrange(id)

tind <- as.factor(as.vector(tind$temp))
# Test against simulated dataset.
jags.data = list("initial"= s$initial,
                 "killed" = s$killed,
                 "id" = s$id,
                 "num.ind" = length(unique(s$id)),
                 "n" = length(s$initial), 
                 "t" = s$temp, 
                 "Ntreats" = length(unique(s$temp)), 
                 "tind" = tind,
                 "T" = 24
) # named list


n.chains = 3
n.burnin = 10000
n.thin = 2
n.iter = 25000
model = R2jags::jags(jags.data, parameters.to.save=jags.params,inits=NULL,
                     model.file=model.loc, n.chains = n.chains, n.burnin=n.burnin,
                     n.thin=n.thin, n.iter=n.iter, DIC=TRUE)

a = MCMCsummary(model,params=c('a'), round = 4)
h = MCMCsummary(model,params=c('h'), round = 4)

t.a = MCMCsummary(model,params=c('t.a'), round = 4)
t.h = MCMCsummary(model,params=c('t.h'), round = 4)

mu.a = MCMCsummary(model,params='mu.a')
mu.h = MCMCsummary(model,params='mu.h')

ids <- sort(unique(sam$lobster_id))

#c(bottom, left, top, right)

png("figures/Sam_fits.png", width = 2000, height = 2000, res = 150)
d <- par(mfrow = c(5,5), mar = c(4,4,1,1))
for(i in 1:22){
  plot(I(killed/24) ~ jitter(initial),data = s[as.numeric(s$id) == i, ], xlab="Number of prey offered",ylab="Consumption rate", main = paste(ids[i]), ylim = c(0,1.6))
  curve(holling2(x,a[i,1],h[i,1],P=1,T=1),add=TRUE,lty=1) #true curve
  text(x = 10, y = 60, label = paste("a =", round(a[i,1], 3), "\n", "h = ", round(h[i,1], 3)))
  abline(a = 0,b = 1, lty = 4)
}
par(d)
dev.off()

d <- par(mfrow = c(2,2), mar = c(4,4,1,1))
for(i in 1:4){
  plot(I(killed/24) ~ jitter(initial),data = s[as.numeric(s$temp) == i, ], xlab="Number of prey offered",ylab="Consumption rate", ylim = c(0,1.6), main = paste(levels(s$temp)[i]))
  curve(holling2(x,t.a[i,1],t.h[i,1],P=1,T=1),add=TRUE,lty=1) #true curve
  # text(x = 10, y = 60, label = paste("a =", round(a[i,1], 3), "\n", "h = ", round(h[i,1], 3)))
  # abline(a = 0,b = 1, lty = 4)
}
par(d)

MCMCplot(model, params = c("a"), rank= T)
MCMCplot(model, params = c("h"), rank= T)


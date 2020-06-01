###############################################################
## test on Sam's data
###############################################################
library(here)
source(here("code","functions.R"))
source(here("code","setup.R"))

sam <- read.csv(here("data/samdata", "fr_data.csv"))
sam$id <- as.numeric(sam$lobster_id)

sam$temp2 <- as.factor(paste("t", sam$temp, sep = ""))

s <- sam[,c( "temp2", "lobster_id","Initial", "Killed")]
names(s) <- c("temp", "id", "initial", "killed")

s <- arrange(s, temp, id, initial)

lattice:: xyplot(killed~initial|as.factor(id), data = s)


jagsscript = cat("

model{

# This model assumed 4 temperature treatments an 11, 16, 21, 26. This model has three levels of heirarchy: a population level estiamte of the FR, a treatment level estimate of the FR, and an indivual level estiamte of the FR.

# PRIORS

# Treatment level effects
for(i in 1:Ntreats){
#a
t.logit.a[i] ~ dnorm(0, t.tau.a)
t.a[i] <- exp(t.logit.a[i])/(1+exp(t.logit.a[i]))

#h
t.h[i] ~ dlnorm(0, tau_int.h)

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

for(i in 1:n){

prob[i] <- max(0.0001,min(0.9999,T/(1/t.a[t[i]] + t.h[t[i]]*initial[i])))

killed[i] ~ dbin(prob[i],initial[i])

}

# Estimate the scaling relationship
  for(i in Ntreats){
    y[i] ~ dnorm(t.a[i], taub)
    log(t.a[i]) <- b*x[i]
  }
# priors
  b ~ dnorm(0, 0.001)


# Variances for all levels
sigma_int.a ~dunif(0,10)
tau_int.a <- 1/(sigma_int.a*sigma_int.a)
sigma_int.h ~dunif(0,10)
tau_int.h <- 1/(sigma_int.h*sigma_int.h)

sigma_t.a ~ dunif(0,10)
t.tau.a <- 1/(sigma_t.a*sigma_t.a)
sigma_t.h ~ dunif(0,10)
t.tau.h <- 1/(sigma_t.h*sigma_t.h)

sigmab ~ dunif(0,10)
taub <- 1/(sigmab*sigmab)

}
                 

                 ", file = here("code", "heirarchical_jagsSAMtreatment.txt"))


model.loc=here("code", "heirarchical_jagsSAMtreatment.txt")
jags.params=c("t.a", "t.h", "b")


# Test against simulated dataset.
jags.data = list("initial"= s$initial,
                 "killed" = s$killed,
                 "id" = s$id,
                 "num.ind" = length(unique(s$id)),
                 "n" = length(s$initial), 
                 "t" = s$temp, 
                 "Ntreats" = length(unique(s$temp)), 
                 "tind" = as.factor(as.vector(distinct(s, temp, id)[,1])), 
                 "T" = 24, 
                 "x" = c(-40.84099, -40.13476, -39.45255, -38.79314)
) # named list


n.chains = 3
n.burnin = 10000
n.thin = 2
n.iter = 25000
model = R2jags::jags(jags.data, parameters.to.save=jags.params,inits=NULL,
                     model.file=model.loc, n.chains = n.chains, n.burnin=n.burnin,
                     n.thin=n.thin, n.iter=n.iter, DIC=TRUE)

a = MCMCsummary(model,params=c('a', 't.a'), round = 4)
h = MCMCsummary(model,params=c('h', 't.h'), round = 4)

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

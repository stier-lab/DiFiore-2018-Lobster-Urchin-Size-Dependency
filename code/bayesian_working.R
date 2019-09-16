########################################################################
## Simulate data
########################################################################

mortrisk <- function(N0,h,a){
  risk <- 1/(1/a + h*N0)                ## Holling risk (per capita)
  pmin(1.0,risk)                         ## bound risk <= 1
}


set.seed(18) ## set random-number seed for reproducibility
## The data will simulate a case where initial attack rate "a" varies
## across 6 replicate blocks
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

ah <- distinct(simdata(10)[,4:5])
## Plot results ...
with(x,plot(jitter(N0),killed,col=ind))

ggplot(x, aes(x = N0, y = killed))+
  geom_jitter()+
  geom_line(aes(col = as.factor(ind)))




## Now fit the data using maximum likelihood without block effects
eps <- 1e-4 ## used to bound probabilities between 0 and 1
##Using the classical unstructured Type II functional response
classic <- mle2(killed~dbinom(prob=pmax(eps,
                                        pmin(1-eps,1/(1/a+h*N0))),
                              size=N0),start=list(a=.01,h=.001),data=x)


names(x) <- c("initial", "ind", "killed")
library("lattice")

xyplot(killed~initial|as.factor(ind), data = x)

#################################################################
## Model 3: Heirarchical bayesian model -- individual fits
#################################################################

jagsscript = cat("model{
                 
                 # hyperpriors
                 #a
                 #mu.a ~ dgamma(shape.a, 1/scale.a) #mean hyperparameter for random a
                 logit_mu.a ~ dnorm(0, 1)
                 mu.a <- exp(logit_mu.a)/(1+exp(logit_mu.a))
                 sigma_int.a ~dunif(0,10) #SD hyperparatmer for random a
                 tau_int.a <- 1/(sigma_int.a*sigma_int.a)
                 
                 #h
                 #mu.h ~ dgamma(shape.h, 1/scale.h) #mean hyperparameter for random h
                 logit_mu.h ~ dnorm(0, 1)
                 mu.h <- exp(logit_mu.h)/(1+exp(logit_mu.h))
                 sigma_int.h ~dunif(0,10) #SD hyperparatmer for random h
                 tau_int.h <- 1/(sigma_int.h*sigma_int.h)
                 
                 # priors
                 #Individual attack rates and handling times vary according to a log normal distribution. 
                 
                 for(i in 1:num.ind){
                 loga[i] ~ dnorm(logit_mu.a, tau_int.a)
                 logh[i] ~ dnorm(logit_mu.h, tau_int.h)
                 a[i] <- exp(loga[i])/(1+exp(loga[i]))
                 h[i] <- exp(logh[i])/(1+exp(logh[i]))# exp(max(min(logh[i],20),-20))
                 }
                 
                 #likelihood
                 
                 for(i in 1:n){
                 
                 prob[i] <- max(0.0001,min(0.9999,1/(1/a[id[i]] + h[id[i]]*initial[i])))
                 
                 killed[i] ~ dbin(prob[i],initial[i])
                 
                 }
                 
                 }", file = "heirarchical_jags3.txt")

model.loc=here("code","heirarchical_jags3.txt")
jags.params=c("a", "h", "mu.a", "mu.h", "sigma_int.a", "sigma_int.h")

# model.loc=here("code", "BH_jags.txt")
# jags.params=c("a", "h")



# Test against simulated dataset.
jags.data = list("initial"= x$initial,
                 "killed" = x$killed,
                 # "P" = 1,
                 # "T" = 1,
                 "id" = x$ind,
                 "num.ind" = length(unique(x$ind)),
                 "n" = length(x$initial), 
                 # "scale.a" = scale.a,
                 # "shape.a" = shape.a,
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

out <- model$BUGSoutput

out$mean
a = MCMCsummary(model,params='a') # the alpha estimate here is often bounding up against zero
h = MCMCsummary(model,params='h')

mu.a = MCMCsummary(model,params='mu.a')
mu.h = MCMCsummary(model,params='mu.h')

plot(killed ~ jitter(initial),data = x, xlab="Number of prey",ylab="Number consumed", ylim = c(0,26))
curve(holling2(x,mu.a[1], mu.h[1],P=1,T=1),add=TRUE,col=1,lty=1) #true curve

temp <- data.frame(distinct(df, id, size, as.numeric(id))) 
names(temp) <- c("id.car", "size", "id" )
df.a <- temp %>% left_join(data.frame(id = seq(1,46), a)) %>% mutate(parameter = "a")
df.h <- temp %>% left_join(data.frame(id = seq(1,46), h)) %>% mutate(parameter = "h")

model.output <- rbind(df.a, df.h)


for(i in 1:10){
  plot(killed ~ jitter(initial),data = x[as.numeric(x$ind) == i, ], xlab="Number of prey",ylab="Number consumed", ylim = c(0,26))
  curve(holling2(x,a[i],h[i],P=1,T=1),add=TRUE,col=1,lty=1) #true curve
}



##############################################################
## Data cloning
##############################################################

dc <- as.data.frame(sapply(x, rep.int, times=100))


jags.data = list("initial"= dc$initial,
                 "killed" = dc$killed,
                 # "P" = 1,
                 # "T" = 1,
                 "id" = dc$ind,
                 "num.ind" = length(unique(dc$ind)),
                 "n" = length(dc$initial), 
                 "scale.a" = scale.a, 
                 "shape.a" = shape.a, 
                 "scale.h" = scale.h, 
                 "shape.h" = shape.h) # named list

n.chains = 3
n.burnin = 10000
n.thin = 2
n.iter = 30000

model = jags(jags.data,parameters.to.save=jags.params,inits=NULL,
             model.file=model.loc, n.chains = n.chains, n.burnin=n.burnin,
             n.thin=n.thin, n.iter=n.iter, DIC=TRUE)



###################################
## ML Nonlinear w/ random effects
###################################

library(nlme)

x$ind <- as.factor(x$ind)

fm1 <- nlme(killed ~ 1/(1/a + h*initial),
            data = x,
            fixed = a + h ~ 1,
            # random = 1 ~ 1 | ind,
            groups = ~ ind,
            start = c(a = 0.51 , h = 0.05))
summary(fm1)
fm2 <- update(fm1, random = pdDiag(Asym + lrc ~ 1))
summary(fm2)





#####################################
## Get data
#####################################

df <- read.table(here("data/cleaned","loburc_cleaned.csv"), header = T, sep = ",")


jagsscript = cat("model{
                 
                 # hyperpriors
                 #a
                 logit_mu.a ~ dnorm(0, 1)
                 mu.a <- exp(logit_mu.a)/(1+exp(logit_mu.a))
                 
                 #h
                 logit_mu.h ~ dnorm(0, 1)
                 mu.h <- exp(logit_mu.h)/(1+exp(logit_mu.h))
              
                 # Variance parameters
                 sigma_int.a ~dunif(0,10) #SD hyperparatmer for random a
                 tau_int.a <- 1/(sigma_int.a*sigma_int.a)
                 sigma_int.h ~dunif(0,10) #SD hyperparatmer for random h
                 tau_int.h <- 1/(sigma_int.h*sigma_int.h)


                 # priors
                 #Individual attack rates and handling times vary according to a log normal distribution. 

                 for(i in 1:num.treatments){
                 loga.r[i] ~ dnorm(logit_mu.a, tau_int.a)
                 logh.r[i] ~ dnorm(logit_mu.h, tau_int.h)
                 a.r[i] <- exp(loga.r[i])/(1+exp(loga.r[i]))
                 h.r[i] <- exp(logh.r[i])/(1+exp(logh.r[i]))
                 }
                 
                 for(i in 1:num.ind){
                 loga[i] ~ dnorm(a.r[treatment[i]], tau_int.a)
                 logh[i] ~ dnorm(h.r[treatment[i]], tau_int.h)
                 a[i] <- exp(loga[i])/(1+exp(loga[i]))
                 h[i] <- exp(logh[i])/(1+exp(logh[i]))# exp(max(min(logh[i],20),-20))
                 }
                 
                 #likelihood
                 
                 for(i in 1:n){
                 
                 prob[i] <- max(0.0001,min(0.9999,1/(1/a[id[i]] + h[id[i]]*initial[i])))
                 
                 killed[i] ~ dbin(prob[i],initial[i])
                 
                 }
                 
                 }", file = here("code", "heirarchical_jags4.txt"))



model.loc=here("code","heirarchical_jags4.txt")
jags.params=c("a", "h", "a.r", "h.r", "mu.a", "mu.h", "sigma_int.a", "sigma_int.h")



killed <- df$num_consumed
initial <- df$num_offered
id <- df$id
treatment <- df$treatment

# Test against simulated dataset.
jags.data = list("initial"= initial,
                 "killed" = killed,
                 "treatment" = treatment,
                 "id" = id,
                 "num.ind" = length(unique(id)),
                 "n" = length(initial), 
                 "num.treatments" = length(unique(treatment))
) # named list

n.chains = 3
n.burnin = 5000
n.thin = 10
n.iter = 10000
model = R2jags::jags(jags.data,parameters.to.save=jags.params,inits=NULL,
                     model.file=model.loc, n.chains = n.chains, n.burnin=n.burnin,
                     n.thin=n.thin, n.iter=n.iter, DIC=TRUE)

a = MCMCsummary(model,params='a') # the alpha estimate here is often bounding up against zero
h = MCMCsummary(model,params='h')

mu.a = MCMCsummary(model,params='mu.a')
mu.h = MCMCsummary(model,params='mu.h')

ids <- unique(df$id)

par(mfrow = c(3,3), mar = c(4,4, 1,1))
for(i in 1:46){
  plot(num_consumed ~ jitter(num_offered),data = df[as.numeric(df$id) == i, ], xlab="Number of prey",ylab="Number consumed", ylim = c(0,26), main = paste(ids[i]))
  curve(holling2(x,a[i,1],h[i,1],P=1,T=1),add=TRUE,col=1,lty=1) #true curve
}

mod <- distinct(df, id, .keep_all = T)
mod$a <- a[,1]
mod$h <- h[,1]


ggplot(mod, aes(x = size, y = a))+ 
  geom_point()+
  facet_wrap(~treatment)

ggplot(mod, aes(x = size, y = h))+ 
  geom_point()+
  facet_wrap(~treatment)

plot(a ~ size, mod)
plot(I(1/h) ~ size, mod)


MCMCplot(model, params = "a", rank = T)
MCMCplot(model, params = "h", rank = T)


###############################################################
## test on Sam's data
###############################################################

sam <- read.csv("~/Downloads/fr_data.csv")
sam$id <- as.numeric(sam$lobster_id)

s <- sam[,c("Initial", "Killed", "id")]
names(s) <- c("initial", "killed", "ind")

lattice:: xyplot(killed~initial|as.factor(ind), data = s)

killed <- s$killed
initial <- s$initial
id <- s$ind

model.loc=here("code","heirarchical_jags4.txt")
jags.params=c("a", "h", "mu.a", "mu.h", "sigma_int.a", "sigma_int.h")

# model.loc=here("code", "BH_jags.txt")
# jags.params=c("a", "h")



# Test against simulated dataset.
jags.data = list("initial"= initial,
                 "killed" = killed,

                 "id" = id,
                 "num.ind" = length(unique(id)),
                 "n" = length(initial)#, 
                 # "scale.a" = scale.a,
                 # "shape.a" = shape.a,
                 # "scale.h" = scale.h,
                 # "shape.h" = shape.h
) # named list

n.chains = 3
n.burnin = 10000
n.thin = 2
n.iter = 50000
model = R2jags::jags(jags.data, parameters.to.save=jags.params,inits=NULL,
                     model.file=model.loc, n.chains = n.chains, n.burnin=n.burnin,
                     n.thin=n.thin, n.iter=n.iter, DIC=TRUE)

a = MCMCsummary(model,params='a') # the alpha estimate here is often bounding up against zero
h = MCMCsummary(model,params='h')

mu.a = MCMCsummary(model,params='mu.a')
mu.h = MCMCsummary(model,params='mu.h')

ids <- sort(unique(sam$lobster_id))

#c(bottom, left, top, right)

png("figures/Sam_fits.png", width = 2000, height = 2000, res = 150)
d <- par(mfrow = c(5,5), mar = c(4,4,1,1))
for(i in 1:22){
  plot(killed ~ jitter(initial),data = s[as.numeric(s$ind) == i, ], xlab="Number of prey",ylab="Number consumed", ylim = c(0,60), main = paste(ids[i]))
  curve(holling2(x,a[i,1],h[i,1],P=1,T=1),add=TRUE,col=1,lty=1) #true curve
  abline(a = 0,b = 1, lty = 4)
}
par(d)
dev.off()

ms <- read.csv("~/Downloads/nls_multstart_params.csv")
msa <- ms[ms$term == "a",c("lobster_id", "term", "estimate")]

jags <- data.frame(jags.a = a[,1], lobster_id = ids)

jags %>% left_join(msa)


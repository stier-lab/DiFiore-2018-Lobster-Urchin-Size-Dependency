
####################################
## Source prior code 
####################################

source(here("code","setup.R"))
source(here("code","functions.R"))
source(here("code","Dunn_priorgeneration.R"))


#####################################
## Get data
#####################################

df <- read.table(here("data/cleaned","loburc_cleaned.csv"), header = T, sep = ",")


#Write Jags model

jagsscript = cat("
                 
                 model{
                 
                 #likelihood
                 
                 for(i in 1:n){
                 
                 #prob[i] <- a*P*T/(1+a*h*initial[i])
                 prob[i] <- max(0.0001,min(0.9999,a*P*T/(1+a*h*initial[i])))
                 killed[i] ~ dbin(prob[i],initial[i])
                 
                 }
                 
                 #priors (lots of different priors...)
                 # a ~ dunif(0.001, 2)
                 # h ~ dunif(0.001, 2)
                
                 #priors based on Dunn et al. 2018 (Ecology)
                 a ~ dgamma(shape.a, 1/scale.a)
                 h ~ dgamma(shape.h, 1/scale.h)
                 
                 # a ~ dnorm(0, 0.01)
                 # h ~ dnorm(0, 0.01)
                 
                 # h ~ dlnorm(0,0.01)
                 # a ~ dlnorm(0,1)
                 
                 # a ~ dgamma(0.01, 0.01)
                 # h ~ dgamma(0.01, 0.01)
                 
                 # a ~ dunif(0.2, 0.3)
                 # h ~ dunif(0,0.2)
                 
                 }
                 
                 ",file="BH_jags.txt")

model.loc=here("code","BH_jags.txt") # name of the txt file
jags.params=c("a", "h")


n.chains = 3
n.burnin = 5000
n.thin = 10
n.iter = 10000

########################################
## Data fits
########################################

#global fit

killed <- df$num_consumed
initial <- df$num_offered

#data going into the model
jags.data = list("initial"= initial,
                 "killed" = killed,
                 "P" = 1, 
                 "T" = 1, 
                 n = length(initial), 
                 scale.a = scale.a, 
                 shape.a = shape.a, 
                 scale.h = scale.h, 
                 shape.h = shape.h) # named list


model.overall = jags(jags.data,parameters.to.save=jags.params,inits=NULL,
             model.file=model.loc, n.chains = n.chains, n.burnin=n.burnin,
             n.thin=n.thin, n.iter=n.iter, DIC=TRUE)


# Write function to fit jags code to subsets

fit.jags <- function(pred_size, urc_size, n.chains = 3, n.burnin = 10000, n.thin = 2, n.iter = 20000){
  temp <- df %>% filter(class == pred_size, treatment == urc_size)
  killed <- temp$num_consumed
  initial <- temp$num_offered
  
  jags.data = list("initial"= initial,
                   "killed" = killed,
                   "P" = 1, 
                   "T" = 1, 
                   n = length(initial), 
                   scale.a = scale.a, 
                   shape.a = shape.a, 
                   scale.h = scale.h, 
                   shape.h = shape.h) # named list
  
  jags(jags.data,parameters.to.save=jags.params,inits=NULL,
               model.file=model.loc, n.chains = n.chains, n.burnin=n.burnin,
               n.thin=n.thin, n.iter=n.iter, DIC=TRUE)
}



# Fit to each size/treatment combination and plot.

mod1 <- fit.jags("jumbo", "urc_large")
mod2 <- fit.jags("large", "urc_large")
mod3 <- fit.jags("medium", "urc_large")
mod4 <- fit.jags("small", "urc_large")
mod5 <- fit.jags("jumbo", "urc_medium")
mod6 <- fit.jags("large", "urc_medium")
mod7 <- fit.jags("medium", "urc_medium")
mod8 <- fit.jags("small", "urc_medium")
mod9 <- fit.jags("jumbo", "urc_small")
mod10 <- fit.jags("large", "urc_small")
mod11 <- fit.jags("medium", "urc_small")
mod12 <- fit.jags("small", "urc_small")


#################################################################
## Build heirarchical bayesian model
#################################################################

# The goal is to build a mixed effects bayesian model that estimates a population level attack rate and handling time, and individual attack rates and handling times for each lobster. 

#Build list with all necessary data contained as vectors in named slots in the list

id <- as.numeric(df$id)


dat <- list(initial = initial, killed = killed, id = id, n = length(initial), num.ind = length(unique(id)))

jagsscript = cat("
                 
                 model{
                 
                 #likelihood
                 
                 for(i in 1:n){
                 
                 prob[i] <- max(0.0001,min(0.9999,ar[id[i]]*P*T/(1+ar[id[i]]*ht[id[i]]*initial[i])))
                 killed[i] ~ dbin(prob[i],initial[i])
                 
                  
                 }
                 
                 #Individual attack rates and handling times vary according to a normal distribution. 
                
                  for(i in 1:num.ind){
                      ar[i] ~ dlnorm(a, tau1)
                      ht[i] ~ dlnorm(h, tau2)
                 }
                 
                
                
                 #Priors on the individual-level variation 
                 tau1 <- 1/ (sd * sd)
                 sd ~ dunif(0, 50)
                 
                 tau2 <- 1/ (sd * sd)
                 sd2 ~ dunif(0, 50)


                 # hyperpriors based on Dunn et al. 2018 (Ecology)
                 a ~ dgamma(shape.a, 1/scale.a)
                 h ~ dgamma(shape.h, 1/scale.h)
                 
                 
                 }
                 
                 ",file=here("code", "heirarchical_jags.txt"))

model.loc=here("code","heirarchical_jags.txt")
jags.params=c("a", "h", "ar", "ht")

jags.data = list("initial"= dat$initial,
                 "killed" = dat$killed,
                 "P" = 1, 
                 "T" = 1,
                 "id" = dat$id,
                 "num.ind" = dat$num.ind,
                 "n" = dat$n, 
                 "scale.a" = scale.a, 
                 "shape.a" = shape.a, 
                 "scale.h" = scale.h, 
                 "shape.h" = shape.h) # named list

n.chains = 3
n.burnin = 10000
n.thin = 2
n.iter = 25000
model = jags(jags.data,parameters.to.save=jags.params,inits=NULL,
             model.file=model.loc, n.chains = n.chains, n.burnin=n.burnin,
             n.thin=n.thin, n.iter=n.iter, DIC=TRUE)


#traceplot(model)






####################################
## Source prior code 
####################################

source(here("code","setup.R"))
source(here("code","functions.R"))
source(here("code","Dunn_priorgeneration.R"))


#####################################
## Get data
#####################################

df <- read.table("data/cleaned/loburc_cleaned.csv", header = T, sep = ",")


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
                 
                 ",file="code/BH_jags.txt")

model.loc="code/BH_jags.txt" # name of the txt file
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


model = jags(jags.data,parameters.to.save=jags.params,inits=NULL,
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

jagsscript = cat("
                 
                 model{
                 
                 #likelihood
                 
                 for(i in 1:n){
                 
                 #prob[i] <- a*P*T/(1+a*h*initial[i])
                 prob[i] <- max(0.0001,min(0.9999,a*P*T/(1+a*h*initial[i])))
                 killed[i] ~ dbin(prob[i],initial[i])
                 
                 }
                 
                 #hyperpriors


                 #priors based on Dunn et al. 2018 (Ecology)
                 a ~ dgamma(shape.a, 1/scale.a)
                 h ~ dgamma(shape.h, 1/scale.h)
                 
                 
                 }
                 
                 ",file="code/heirarchical_jags.txt")

model.loc="code/heirarchical_jags.txt"

jagsscript = cat("
model {  
                 U ~ dnorm(0, 0.01); #mean population growth rate
                 tauQ~dgamma(0.001,0.001); #error associated with growth rate
                 Q <- 1/tauQ;
                 
                 # Estimate the initial population abundance
                 X[1] ~ dnorm(3,0.01); # vague normal prior 
                 
                 # Autoregressive process for remaining years
                 for(i in 2:nYears) {
                 predX[i] <- X[i-1] + U; # our population growth model in log space
                 X[i] ~ dnorm(predX[i], tauQ); # adding the error to the growth model for each time step
                 }
                 
                 # Observation model
                 # The Rs are different in each survey (i.e. the error is different)
                 for(i in 1:nSurveys) {
                 tauR[i]~dgamma(0.001,0.001); #different error terms (R's) for each survey
                 R[i] <- 1/tauR[i];
                 }
                 for(i in 1:nYears) {
                 for(j in 1:nSurveys) {
                 Y[i,j] ~ dnorm(X[i],tauR[j]); # relate the latent state variable to the known observations allowing their errors to vary independently
                 }
                 }
                 }  
                 
                 ", 
                 file = "marss_lingcod_true.txt")












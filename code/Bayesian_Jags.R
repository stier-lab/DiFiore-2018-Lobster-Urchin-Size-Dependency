####################################
## Libraries
####################################

library(R2jags)
library(rjags)
library(MCMCvis)


#####################################
## Get data
#####################################

df <- read.table("data/cleaned/loburc_cleaned.csv", header = T, sep = ",")


#Write Jags model

jagsscript = cat("
                 
                 model{
                 
                 #likelihood
                 
                 for(i in 1:length(n)){
                 
                 #prob[i] <- a*P*T/(1+a*h*initial[i])
                 prob[i] <- max(0.0001,min(0.9999,a*P*T/(1+a*h*initial[i])))
                 killed[i] ~ dbin(prob[i],initial[i])
                 
                 }
                 
                 #priors (lots of different priors...)
                 a ~ dunif(0.001, 1)
                 h ~ dunif(0.001, 1)
                 
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
                 n = length(initial)) # named list


model = jags(jags.data,parameters.to.save=jags.params,inits=NULL,
             model.file=model.loc, n.chains = n.chains, n.burnin=n.burnin,
             n.thin=n.thin, n.iter=n.iter, DIC=TRUE)

# Write function to fit jags code to subsets

fit.jags <- function(pred_size, urc_size, n.chains = 3, n.burnin = 5000, n.thin = 10, n.iter = 10000){
  temp <- df %>% filter(class == pred_size, treatment == urc_size)
  killed <- temp$num_offered
  initial <- temp$num_consumed
  
  jags(jags.data,parameters.to.save=jags.params,inits=NULL,
               model.file=model.loc, n.chains = n.chains, n.burnin=n.burnin,
               n.thin=n.thin, n.iter=n.iter, DIC=TRUE)
}

plot.jags <- function(model){
  temp <- data.frame(Killed = killed, Initial = initial)
  plot(Killed ~ Initial,data = temp, xlab="initial density",ylab="final density")
  a = MCMCsummary(model,params='a')[1] # the alpha estimate here is often bounding up against zero
  h = MCMCsummary(model,params='h')[1]
  curve(holling2(x,a,h,P=1,T=1),add=TRUE,col=1,lty=1) #true curve
  
}

# Write for loop to estimate parameters for each size combination

mod1 <- fit.jags("jumbo", "urc_large")
plot.jags(mod1)
mod2 <- fit.jags("large", "urc_large")
plot.jags(mod2)
mod3 <- fit.jags("medium", "urc_large")
plot.jags(mod3)
mod4 <- fit.jags("small", "urc_large")
plot.jags(mod4)
mod5 <- fit.jags("jumbo", "urc_medium")
plot.jags(mod5)
mod6 <- fit.jags("large", "urc_medium")
plot.jags(mod6)
mod7 <- fit.jags("medium", "urc_medium")
plot.jags(mod7)
mod8 <- fit.jags("small", "urc_medium")
plot.jags(mod8)
mod9 <- fit.jags("jumbo", "urc_small")
plot.jags(mod9)
mod10 <- fit.jags("large", "urc_small")
plot.jags(mod10)
mod11 <- fit.jags("medium", "urc_small")
plot.jags(mod11)
mod12 <- fit.jags("small", "urc_small")
plot.jags(mod12)
























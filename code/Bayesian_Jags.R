
####################################
## Source prior code 
####################################

library(here)
source(here("code","setup.R"))
source(here("code","functions.R"))
source(here("code","Dunn_priorgeneration.R"))
source(here("code", "Simulated_datasetforJags.R"))


#####################################
## Get data
#####################################

df <- read.table(here("data/cleaned","loburc_cleaned.csv"), header = T, sep = ",")

#############################################
# Model 1: Basic
#############################################

# Model 1 is a basic bayesian model that estimates the number of individuals consumed. It does not have any hierarchical structure. 

jagsscript = cat("
                 
                 model{
                 
                 #likelihood
                 
                 for(i in 1:n){
                 
                 max(0.0001,min(0.9999,1/(1/a+h*initial[i])))
                 #prob[i] <- max(0.0001,min(0.9999,a*P*T/(1+a*h*initial[i])))
                 killed[i] ~ dbin(prob[i],initial[i])
                 
                 }
                 
                
                 #priors based on Dunn et al. 2018 (Ecology)
                 a ~ dgamma(shape.a, 1/scale.a)
                 h ~ dgamma(shape.h, 1/scale.h)
                 
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
id <- as.numeric(df$id)

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
model.overall

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

###################################################################################
## Model 2: Heirarchical bayesian model -- addative random effect of individual
###################################################################################

jagsscript = cat("
                 
                 model{
                 
                 #likelihood
                 
                 for(i in 1:n){
                 
                 prob[i] <- max(0.0001,min(0.9999,(1/(1/a+h*initial[i]) + ind[id[i]]) ))
                 killed[i] ~ dbin(prob[i],initial[i]) 
                 
                 }

                 for(j in 1:num.ind){
                    ind[j] ~ dnorm(0, tau)
                 }
                 
                 tau <- 1/ (sd * sd)
                 sd ~ dgamma(0.01, 0.01)
                 
                 #priors based on Dunn et al. 2018 (Ecology)
                 a ~ dgamma(shape.a, 1/scale.a)
                 h ~ dgamma(shape.h, 1/scale.h)
                 
                 }
                 
                 ",file=here("code", "Heirarchical_additive.txt"))

model.loc=here("code","Heirarchical_additive.txt") # name of the txt file
jags.params=c("a", "h")

#data going into the model
jags.data = list("initial"= initial,
                 "killed" = killed,
                 "id" = id,
                 n = length(initial),
                 num.ind = length(unique(id)),
                 scale.a = scale.a, 
                 shape.a = shape.a, 
                 scale.h = scale.h, 
                 shape.h = shape.h) # named list

n.chains = 3
n.burnin = 10000
n.thin = 2
n.iter = 25000
model.additive = jags(jags.data,parameters.to.save=jags.params,inits=NULL,
                     model.file=model.loc, n.chains = n.chains, n.burnin=n.burnin,
                     n.thin=n.thin, n.iter=n.iter, DIC=TRUE)


#################################################################
## Model 3: Heirarchical bayesian model -- individual fits
#################################################################

# The goal is to build a mixed effects bayesian model that estimates a population level attack rate and handling time, and individual attack rates and handling times for each lobster. 

#Build list with all necessary data contained as vectors in named slots in the list

jagsscript = cat("
                 
                 model{
                 
                 #likelihood
                 
                 for(i in 1:n){
                 
                 prob[i] <- max(0.0001,min(0.9999,1/(1/a[id[i]] + h[id[i]]*initial[i])))
                 #prob[i] <- max(0.0001,min(0.9999,a[id[i]]*P*T/(1+a[id[i]]*h[id[i]]*initial[i])))
                 killed[i] ~ dbin(prob[i],initial[i])
                 
                  
                 }
                 
                 #Individual attack rates and handling times vary according to a normal distribution. 
                
                  for(i in 1:num.ind){
                      a[i] ~ dnorm(mu.a, tau1)
                      h[i] ~ dnorm(mu.h, tau2)
                 }
                 
                
                
                 #Priors on the individual-level variation 
                 tau1 <- 1/ (sd * sd)
                 sd ~ dgamma(0.01, 0.01)
                 
                 tau2 <- 1/ (sd2 * sd2)
                 sd2 ~ dgamma(0.01, 0.01)


                 # hyperpriors based on Dunn et al. 2018 (Ecology)
                 mu.a ~ dgamma(shape.a, 1/scale.a)
                 mu.h ~ dgamma(shape.h, 1/scale.h)
                 
                 
                 }
                 
                 ",file=here("code", "heirarchical_jags.txt"))

model.loc=here("code","heirarchical_jags.txt")
jags.params=c("mu.a", "mu.h", "a", "h")

jags.data = list("initial"= initial,
                 "killed" = killed,
                 # "P" = 1, 
                 # "T" = 1,
                 "id" = id,
                 "num.ind" = length(unique(id)),
                 "n" = length(initial), 
                 "scale.a" = scale.a, 
                 "shape.a" = shape.a, 
                 "scale.h" = scale.h, 
                 "shape.h" = shape.h) # named list

n.chains = 3
n.burnin = 10000
n.thin = 2
n.iter = 10^6
model = jags(jags.data,parameters.to.save=jags.params,inits=NULL,
             model.file=model.loc, n.chains = n.chains, n.burnin=n.burnin,
             n.thin=n.thin, n.iter=n.iter, DIC=TRUE)






# Test against simulated dataset.
jags.data = list("initial"= x$initial,
                 "killed" = x$killed,
                 #"P" = 1, 
                 #"T" = 1,
                 "id" = x$ind,
                 "num.ind" = length(unique(x$ind)),
                 "n" = length(x$initial), 
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





a = MCMCsummary(model,params='a') # the alpha estimate here is often bounding up against zero
h = MCMCsummary(model,params='h')

plot(killed ~ jitter(initial),data = df, xlab="Number of prey",ylab="Number consumed", ylim = c(0,26))
curve(holling2(x,a[1],h[1],P=1,T=1),add=TRUE,col=1,lty=1) #true curve

temp <- data.frame(distinct(df, id, size, as.numeric(id))) 
names(temp) <- c("id.car", "size", "id" )
df.a <- temp %>% left_join(data.frame(id = seq(1,46), a)) %>% mutate(parameter = "a")
df.h <- temp %>% left_join(data.frame(id = seq(1,46), h)) %>% mutate(parameter = "h")

model.output <- rbind(df.a, df.h)

model.output %>% filter(Rhat < 1.3) %>% 
ggplot(aes(x = size, y = mean))+
  geom_point()+
  facet_wrap(~parameter)

for(i in 1:46){
plot(num_offered ~ jitter(num_consumed),data = df[as.numeric(df$id) == i, ], xlab="Number of prey",ylab="Number consumed", ylim = c(0,26))
curve(holling2(x,a[i],h[i],P=1,T=1),add=TRUE,col=1,lty=1) #true curve
}





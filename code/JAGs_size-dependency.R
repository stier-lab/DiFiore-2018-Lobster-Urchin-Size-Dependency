library(here)
source(here("code", "setup.R"))
source(here("code", "functions.R"))


#####################################
## Get data
#####################################


df <- read.table(here("data/cleaned","loburc_cleaned.csv"), header = T, sep = ",")

##################################################################
## the model
##################################################################

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
sigma_t.h ~ dunif(0,2)
t.tau.h <- 1/(sigma_t.h*sigma_t.h)

}", file = here("code", "heirarchical_jagsLU.txt"))



model.loc=here("code","heirarchical_jagsLU.txt")
jags.params=c("a", "h", "t.a", "t.h", "mu.a", "mu.logit.a", "mu.h")




jags.data = list("initial"= df$initial,
                 "killed" = df$killed,
                 "id" = df$id,
                 "num.ind" = length(unique(df$id)),
                 "n" = length(df$initial), 
                 "tind" = as.factor(as.vector(distinct(df, treatment, id)[,2])), 
                 "T" = 48, 
                 "Ntreats" = length(unique(df$treatment))
) # named list

n.chains = 3
n.burnin = 10000
n.thin = 2
n.iter = 25000
model = R2jags::jags(jags.data,parameters.to.save=jags.params,inits=NULL,
                     model.file=model.loc, n.chains = n.chains, n.burnin=n.burnin,n.thin=n.thin, n.iter=n.iter, DIC=TRUE)

a = MCMCsummary(model,params='a', round = 4)
h = MCMCsummary(model,params='h', round = 4)

t.a = MCMCsummary(model,params='t.a', round = 4)
t.h = MCMCsummary(model,params='t.h', round = 4)

ids <- distinct(df, id, size, treatment, udiam) %>% arrange(id)

d <- par(mfrow = c(5,5), mar = c(4,4, 1,1))
for(i in 1:46){
  plot(I(killed/48) ~ jitter(initial),data = df[as.numeric(df$id) == i, ], xlab="Number of prey",ylab="Consumption rate (prey consumed per predator per hour)", ylim = c(0,1), main = paste(as.numeric(unique(df$id))[i]))
  curve(holling2(x,a[i,1],h[i,1],P=1,T=1),add=TRUE,col=1,lty=1) #true curve
  text(x = 10, y = 0.6, label = paste("a = ", round(a[i,1], 3), "\nh =", round(h[i,1], 3), "\n lsize = ", ids[i, 2], "\n usize =", levels(df$id)[i], sep = ""))
}
par(d)


treats <- levels(df$treatment)
d <- par(mfrow = c(1,3), mar = c(4,4, 1,1))
for(i in 1:3){
  plot(I(killed/48) ~ jitter(initial),data = df[as.numeric(df$treatment) == i, ], xlab="Number of prey",ylab="Consumption rate (prey consumed per predator per hour)", ylim = c(0,1), main = paste(treats[i]))
  curve(holling2(x,t.a[i,1],t.h[i,1],P=1,T=1),add=TRUE,col=1,lty=1) #true curve
}
par(d)


mod <- distinct(df, id, size, treatment, mr, mc) %>% arrange(id)
mod$a <- a[,1]
mod$h <- h[,1]
mod$max.intake <- 1/mod$h
levels(mod$treatment) <- c("urc_small", "urc_medium", "urc_large")
#mod <- mod[mod$id != "N07", ]


lm1 <- lm(log(a) ~ log(mr) + log(mc), data = mod)

lm1 <- lm(log(a) ~ log(mr) + log(I(mc/mr)) + I(mc/mr), data = mod)

mod$out <- predict(lm1)
plot(out ~ log(I(mod$mc/mod$mr)))

lm2 <- lm(log(h) ~ log(mr) + log(mc), data = mod)

ggplot(mod, aes(x = mc/mr, y = a))+ 
  geom_point(aes(color = treatment, size = treatment), pch = 21)+
  geom_smooth(method = "lm", aes(color = treatment))

ggplot(mod, aes(x = mc, y = h))+ 
  geom_point(aes(color = treatment, size = treatment), pch = 21)




ggplot(mod, aes(x = mc, y = a))+ 
  geom_point()+
  facet_wrap(~treatment)

ggplot(mod, aes(x = log(mc), y = log(a)))+ 
  geom_point()+
  facet_wrap(~treatment)


ggplot(mod, aes(x = mc, y = h))+ 
  geom_point()+
  facet_wrap(~treatment)

ggplot(mod, aes(x = mc, y = max.intake))+ 
  geom_point()+
  facet_wrap(~treatment, scales = "free")

mod$size.ratio <- mod$mc / mod$mr

ggplot(mod, aes(x = log(size.ratio), y = log(a)))+
  geom_point()










MCMCplot(model, params = "a", rank = T)
MCMCplot(model, params = "h", rank = T)

################################### Attempt 2

killed <- df$num_consumed
initial <- df$num_offered
dat <- data.frame(initial, killed, id = as.factor(id))

dat <- dat[!duplicated(t(apply(dat[c("initial", "id")], 1, sort))), ]
table(dat$initial, dat$id)

dat <- dat %>%
  filter(initial != 27) %>%
  spread(initial, killed)
dat[is.na(dat)] <- 0

# initial <- dat$initial
initial <- sort(unique(df$num_offered))
killed <- as.matrix(dat[, 2:length(names(dat))])

jagsscript = cat("
                 
                 model{
                 
                 #likelihood
                 
                 for (i in 1:num.ind){
                 for(j in 1:num.den){
                 
                 prob[i, j] <- max(0.0001,min(0.9999,a[i]*P*T/(1+a[i]*h[i]*initial[j])))
                 killed[i,j] ~ dbin(prob[i,j],initial[j])
                 
                 }
                 
                 #Individual attack rates and handling times vary according to a normal distribution. 
                 a[i] ~ dnorm(mu.a, tau.a)
                 h[i] ~ dnorm(mu.h, tau.h)
                 }
                 
                 #Priors on the individual-level variation 
                 tau.a <- 1/ (sd * sd)
                 sd ~ dunif(0, 10)
                 
                 tau.h <- 1/ (sd.h * sd.h)
                 sd.h ~ dunif(0, 10)
                 
                 
                 # hyperpriors based on Dunn et al. 2018 (Ecology)
                 mu.a ~ dgamma(shape.a, 1/scale.a)
                 mu.h ~ dgamma(shape.h, 1/scale.h)
                 
                 }
                 
                 ",file=here("code", "heirarchical_jags2.txt"))

model.loc=here("code","heirarchical_jags2.txt")
jags.params=c("mu.a", "mu.h", "a", "h")

jags.data = list("initial"= initial,
                 "killed" = killed,
                 "P" = 1, 
                 "T" = 1,
                 "num.ind" = length(dat$id),
                 "num.den" = ncol(killed),
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









x <- x %>%
  spread(initial, killed)

# initial <- dat$initial
initial <- colnames(x[-1])
killed <- as.matrix(x[, 2:length(names(x))])

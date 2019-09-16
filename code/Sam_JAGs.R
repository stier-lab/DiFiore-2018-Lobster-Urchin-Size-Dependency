###############################################################
## test on Sam's data
###############################################################

sam <- read.csv(here("data/samdata", "fr_data.csv"))
sam$id <- as.numeric(sam$lobster_id)

sam$temp2 <- as.factor(paste("t", sam$temp, sep = ""))

s <- sam[,c( "temp2", "lobster_id","Initial", "Killed")]
names(s) <- c("temp", "id", "initial", "killed")

s <- arrange(s, temp, id, initial)

lattice:: xyplot(killed~initial|as.factor(id), data = s)

killed <- s$killed
initial <- s$initial
id <- s$id
t <- s$temp
tind <- as.vector(distinct(s, temp, id)[,1])

jagsscript = cat("

model{
  
  # This model assumed 4 temperature treatments an 11, 16, 21, 26. This model has three levels of heirarchy: a population level estiamte of the FR, a treatment level estimate of the FR, and an indivual level estiamte of the FR.
  
  # PRIORS
  
  # Global population estiamte
      #a
      logit_mu.a ~ dnorm(0, 1)
      mu.a <- exp(logit_mu.a)/(1+exp(logit_mu.a))

      #h
      logit_mu.h ~ dnorm(0, 1)
      mu.h <- exp(logit_mu.h)/(1+exp(logit_mu.h))

  # Treatment level effects
      for(i in 1:Ntreats){
        #a
        t.logit.a[i] ~ dnorm(mu.a, t.tau.a)
        t.a[i] <- exp(t.logit.a[i])/(1+exp(t.logit.a[i]))

        #h
        t.logit.h[i] ~ dnorm(mu.h, t.tau.h)
        t.h[i] <- exp(t.logit.h[i])/(1+exp(t.logit.h[i]))

      }

  
  # Individual level effects
      for(i in 1:num.ind){

        # a
        loga[i] ~ dnorm(t.a[t[i]], tau_int.a)
        a[i] <- exp(loga[i])/(1+exp(loga[i]))

        # h
        logh[i] ~ dnorm(t.h[t[i]], tau_int.h)
        h[i] <- exp(logh[i])/(1+exp(logh[i]))
      }
  
  # functional response likelihood
  
      for(i in 1:n){
        
        prob[i] <- max(0.0001,min(0.9999,1/(1/a[id[i]] + h[id[i]]*initial[i])))
        
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


# Test against simulated dataset.
jags.data = list("initial"= initial,
                 "killed" = killed,
                 "id" = id,
                 "num.ind" = length(unique(id)),
                 "n" = length(initial), 
                 "t" = t, 
                 "Ntreats" = length(unique(t))
) # named list

n.chains = 3
n.burnin = 10000
n.thin = 2
n.iter = 25000
model = R2jags::jags(jags.data, parameters.to.save=jags.params,inits=NULL,
                     model.file=model.loc, n.chains = n.chains, n.burnin=n.burnin,
                     n.thin=n.thin, n.iter=n.iter, DIC=TRUE)

a = MCMCsummary(model,params='a') # the alpha estimate here is often bounding up against zero
h = MCMCsummary(model,params='h')

mu.a = MCMCsummary(model,params='mu.a')
mu.h = MCMCsummary(model,params='mu.h')

ids <- sort(unique(sam$lobster_id))

#c(bottom, left, top, right)

#png("figures/Sam_fits.png", width = 2000, height = 2000, res = 150)
d <- par(mfrow = c(5,5), mar = c(4,4,1,1))
for(i in 1:22){
  plot(killed ~ jitter(initial),data = s[as.numeric(s$ind) == i, ], xlab="Number of prey",ylab="Number consumed", ylim = c(0,60), main = paste(ids[i]))
  curve(holling2(x,a[i,1],h[i,1],P=1,T=1),add=TRUE,lty=1) #true curve
  text(x = 10, y = 60, label = paste("a =", round(a[i,1], 3), "\n", "h = ", round(h[i,1], 3)))
  abline(a = 0,b = 1, lty = 4)
}
par(d)
#dev.off()



mod <- distinct(sam, lobster_id, temp, id) %>% arrange(id)
mod$a <- a[,1]
mod$h <- h[,1]
mod$max.intake <- 1/mod$h


p1 <- ggplot(mod, aes(x = temp, y = log(a)))+
  geom_point()
p2 <- ggplot(mod, aes(x = temp, y = max.intake))+
  geom_point()

cowplot::plot_grid(p1, p2, ncol = 2)

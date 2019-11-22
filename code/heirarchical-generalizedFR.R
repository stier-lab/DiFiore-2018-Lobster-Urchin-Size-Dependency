flex <- function(initial, b,h,q,T){
  b*initial^(q+1) / (1 + b * h * initial ^(q+1))
}


jagsscript = cat("

model{
# PRIORS

# Population level estimates

mu.logit.b ~ dnorm(0, 0.01)
mu.b <- exp(mu.logit.b)/(1+exp(mu.logit.b))
mu.log.h ~ dnorm(0, 0.01)
mu.h <- exp(max(min(mu.log.h, 20), -20))
# mu.log.q ~ dnorm(0, 0.01)T(0,)
# mu.q <- exp(max(min(mu.log.q, 20), -20))
mu.q ~ dnorm(0, 0.01)T(0,)


# Treatment level effects
for(i in 1:Ntreats){
#a
t.logit.b[i] ~ dnorm(mu.logit.b, t.tau.b)
t.b[i] <- exp(t.logit.b[i])/(1+exp(t.logit.b[i]))

#h
t.log.h[i] ~ dnorm(mu.log.h, t.tau.h)
t.h[i] <- exp(max(min(t.log.h[i],20),-20))

#q 
# t.log.q[i] ~ dnorm(mu.log.q, t.tau.q)
# t.q[i] <- exp(max(min(t.log.q[i],20),-20))
t.q[i] ~ dnorm(mu.q, t.tau.q)T(0,)
}


# Individual level effects
for(i in 1:num.ind){

# a
logb[i] ~ dnorm(t.logit.b[tind[i]], tau_int.b)
b[i] <- exp(logb[i])/(1+exp(logb[i]))

# h
logh[i] ~ dnorm(t.log.h[tind[i]], tau_int.h)
h[i] <- exp(max(min(logh[i],10),-20))

# # q
# logq[i] ~ dnorm(t.log.q[tind[i]], tau_int.q)
# q[i] <- exp(max(min(logq[i],10),-20))
}
# 
# functional response likelihood

for(i in 1:n){

prob[i] <- max(0.0001,min(0.9999,
              ((b[id[i]]*initial[i]^(1 + t.q[treats[i]])*T) /
              (1 + b[id[i]]*h[id[i]]*initial[i]^(1 + t.q[treats[i]])))/initial[i]))

killed[i] ~ dbin(prob[i],initial[i])

}



# Variances for all levels
sigma_int.b ~dunif(0,10)
tau_int.b <- 1/(sigma_int.b*sigma_int.b)
sigma_int.h ~dunif(0,10)
tau_int.h <- 1/(sigma_int.h*sigma_int.h)

sigma_t.b ~ dunif(0,10)
t.tau.b <- 1/(sigma_t.b*sigma_t.b)
sigma_t.h ~ dunif(0,10)
t.tau.h <- 1/(sigma_t.h*sigma_t.h)

sigma_t.q ~ dunif(0,10)
t.tau.q <- 1/(sigma_t.q*sigma_t.q)
# sigma_int.q ~dunif(0,10)
# tau_int.q <- 1/(sigma_int.q*sigma_int.q)
}
", file = here("code", "flexible.txt"))



model.loc=here("code","flexible.txt")
jags.params=c("b", "h", "q", "t.b", "t.h", "t.q", "mu.q", "mu.b", "mu.logit.b", "mu.h")

tind <- distinct(df, newtreat, id) %>%
  arrange(id)

tind <- as.factor(as.vector(tind$newtreat))


jags.data = list("initial"= df$initial,
                 "killed" = df$killed,
                 "id" = df$id,
                 "num.ind" = length(unique(df$id)),
                 "n" = length(df$initial), 
                 "tind" = tind, 
                 "T" = 48, 
                 "Ntreats" = length(unique(df$newtreat)), 
                 "treats" = df$newtreat
) # named list
n.chains = 3
n.burnin = 25000
n.thin = 10
n.iter = 50000
model = R2jags::jags(jags.data,parameters.to.save=jags.params,inits=NULL,
                     model.file=model.loc, n.chains = n.chains, n.burnin=n.burnin,n.thin=n.thin, n.iter=n.iter, DIC=TRUE)

b = MCMCsummary(model,params='b', round = 4)
h = MCMCsummary(model,params='h', round = 4)
#q = MCMCsummary(model, params = "q", round = 4)

t.b = MCMCsummary(model,params='t.b', round = 4)
t.h = MCMCsummary(model,params='t.h', round = 4)
t.q = MCMCsummary(model,params='t.q', round = 4)

mu.b = MCMCsummary(model,params='mu.b', round = 4)
mu.h = MCMCsummary(model,params='mu.h', round = 4)
mu.q = MCMCsummary(model,params='mu.q', round = 4)

ids <- distinct(df, id, size, treatment, udiam) %>% arrange(id)

d <- par(mfrow = c(5,5), mar = c(4,4, 1,1))
for(i in 1:46){
  plot(I(killed/48) ~ jitter(initial),data = df[as.numeric(df$id) == i, ], xlab="Number of prey",ylab="Consumption rate (prey consumed per predator per hour)", ylim = c(0,1), main = paste(as.numeric(unique(df$id))[i]))
  curve(flex(x,b[i,1],h[i,1],t.q[tind[i,1]],T=1),add=TRUE,col=1,lty=1) 
  # text(x = 10, y = 0.6, label = paste("a = ", round(a[i,1], 3), "\nh =", round(h[i,1], 3), "\n lsize = ", ids[i, 2], "\n usize =", levels(df$id)[i], sep = ""))
}
par(d)


# treats <- levels(df$treatment)
# d <- par(mfrow = c(1,3), mar = c(4,4, 1,1))
# for(i in 1:3){
#   plot(I(killed/48) ~ jitter(initial),data = df[as.numeric(df$treatment) == i, ], xlab="Number of prey",ylab="Consumption rate (prey consumed per predator per hour)", ylim = c(0,1), main = paste(treats[i]))
#   curve(flex(x,t.b[i,1],t.h[i,1],t.q[i,1],T=1),add=TRUE,col=1,lty=1)
# }
# par(d)
# 
# plot(I(killed/48) ~ jitter(initial),data = df, xlab="Number of prey",ylab="Consumption rate (prey consumed per predator per hour)", ylim = c(0,1))
# curve(flex(x,mu.b[1],mu.h[1],mu.q[1],T=1),add=TRUE,col=1,lty=1)


cols <- rep(viridis::viridis(3),3)
linew <- rep(c(3,2,1), 3)
linetp <- rep(c(4,3,1), 3)

#png(here("figures", "lob-urc-treatmentlevel-FR.png"), width = 1000*3, height = 333*3, res = 300)
#svg(here("figures", "lob-urc-treatmentlevel-FR.svg"), width = 10, height = 3.33)
d <- par(mfrow = c(1,3), mar = c(4.5,5,3,1), las = 1, bty = "n"
         #, cex.axis= 1.2, cex.lab = 2, cex.main = 2
)

plot(I(killed/48) ~ jitter(initial), data = df[df$treatment == "urc_small", ], ylim = c(0, 0.8), xlim = c(0,27), ylab = "Number consumed", xlab = "", main = "Small urchins (1.0-3.0 cm)")
for(i in 7:9){
  curve(flex(x,t.b[i,1],t.h[i,1],t.q[i,1],T=1),add=TRUE, col = cols[i], lwd = linew[i], lty = linetp[i], to = 26)
  # curve(holling2(x,t.a[i,3],t.h[i,5],P=1,T=48),add=TRUE, col = cols[i], lty = 2, to = 26)
  # curve(holling2(x,t.a[i,5],t.h[i,3],P=1,T=48),add=TRUE, col = cols[i], lty = 2, to = 26)
}

e <- par(mar = c(4.5,2, 3,1))
plot(I(killed/48) ~ jitter(initial), data = df[df$treatment == "urc_medium", ], ylim = c(0, 0.8), xlim = c(0,27), ylab = "", xlab = "Number of prey offered", main = "Medium urchins (3.0-5.0 cm)")
for(i in 4:6){
  curve(flex(x,t.b[i,1],t.h[i,1],t.q[i,1],T=1),add=TRUE, col = cols[i], lwd = linew[i], lty = linetp[i], to = 26)
}

plot(I(killed/48) ~ jitter(initial), data = df[df$treatment == "urc_large", ], ylim = c(0, 0.8), xlim = c(0,27), ylab = "", xlab = "", main = "Large urchins (5.0-7.0 cm)")
for(i in 1:3){
  curve(flex(x,t.b[i,1],t.h[i,1],t.q[i,1],T=1),add=TRUE, col = cols[i], lwd = linew[i], lty = linetp[i], to = 26)
}
legend(x = 5, y = 20, legend = c("Large lobsters (>9.0 cm)", "Medium lobsters (7.0-9.0 cm)", "Small lobsters (5.0-7.0 cm)"), col = cols[1:3], lwd = linew[1:3], lty = linetp[1:3], bty = "n", xjust = 0, yjust = 1)

par(e)
par(d)
#dev.off()





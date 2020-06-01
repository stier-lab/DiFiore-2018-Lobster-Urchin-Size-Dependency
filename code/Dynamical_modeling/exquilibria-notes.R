
library(deSolve)

# creating a function takes the same for as all the functions you use

# within the function called sir.fun, supply arguments t, x, and params
sir.fun <- function (t, x, params) { 
  # starting conditions

  V = x[1]
  #P = x[2]

  r = params[1] 
  k = params[2]
  a = params[3]
  h0 = params[4]
  mv = params[5]
  mp = params[6]
  alphaV = params[7]
  alphaP = params[8]
  # f = params[9]
  # q = params[10]
  P = params[9]
  
  # equations to solve
  # state variables get put into these boxes
  #dVdt = r*V*(1-(V/k)) - a*V*P
  dVdt = r*V*(1-(V/k)) - ((a*V*P) / (1 + (a*exp(h0 + alphaV*log(mv) + alphaP*log(mp))*V)))
  #dPdt = r*P*(1-(V/kp)) - q*P
  #dPdt = f*a*V*P - q*P
  
  # return an output
  return(list(dVdt))
}

# program initial conditions
V0 <- 100
P0 <- 10
initial_values = c(V = V0)

r = 1.01
k = 10000
a = mean(df.ind$median.a)
h0 = 12.20539
mv = max(df.ind$mr, na.rm = T)
mp = min(df.ind$mc, na.rm = T)
alphaV = -2.04202
alphaP = 1.34103
P = 100

parameters = c(r, k, a, h0, mv, mp, alphaV, alphaP, P)

times = seq(0, 60, length.out = 20000)


results1 = lsoda(initial_values, times, sir.fun, parameters)
head(results)

parameters = c(r, k, a, h0, mv = min(df.ind$mr, na.rm = T), mp = max(df.ind$mc, na.rm = T), alphaV, alphaP, P)
results2 = lsoda(initial_values, times, sir.fun, parameters)

# ODE PLOT
# indicating specific columns to plot
plot(x = results[, "time"], y = results[, "V"], type = 'l', col=Ncol, las = 1, lwd=2, xlab = 'Time', ylab = 'Number of Individuals', main = "ODE Solver")
abline(h = 0)
abline(h = 1000, col = "red")
lines(x = results[, "time"], y = results[, "V"], lwd = 2)
lines(x = results2[, "time"], y = results2[, "V"], lwd = 2, col = "green")
legend(x = max(tset)*0.1, y = N*1.2, legend=c('N','S','I','R'),lwd= 2, col=c(Ncol,Scol,Icol,Rcol), horiz=TRUE)

plot(x = tset, y = N.simu1, type = 'l', col=Ncol, las = 1, lwd=2, xlab = 'Time', ylab = 'Number of Individuals', ylim = c(0,N*1.2), main = "for() loop")
abline(h = 0)
lines(x = tset, y = S.simu1, col = Scol, lwd = 2)
lines(x = tset, y = I.simu1, col = Icol, lwd = 2)
lines(x = tset, y = R.simu1, col = Rcol, lwd = 2)
legend(x = max(tset)*0.1, y = N*1.2, legend=c('N','S','I','R'),lwd= 2, col=c(Ncol,Scol,Icol,Rcol), horiz=TRUE)





r = 1.01
k = 10000
a = mean(df.ind$median.a)
h0 = 12.20539
mv = min(df.ind$mr, na.rm = T)
mp = max(df.ind$mc, na.rm = T)
alphaV = -2.04202
alphaP = 1.34103
P = 10
V = 0:26

dvdt <- function(r, V, k, a, P, h0, mv, mp, alphaV, alphaP){
  dVdt = r*V*(1-(V/k)) - ((a*V*P) / (1 + (a*exp(h0 + alphaV*log(mv) + alphaP*log(mp))*V)))
}

# dvdt <- function(a, P, h0, mv, mp, alphaV, alphaP){
#   ((a*V*P) / (1 + (a*h0*(mv^alphaV)*(mp^alphaP)*V)))
# }

out <- dvdt(r = r, V = V, k = k, a= a, P = P, h0 = h0, mv = mv, mp = mp, alphaV = alphaV, alphaP = alphaP)
out2 <- dvdt(r = r, V = V, k = k, a= a, P = P, h0 = h0, mv = max(df.ind$mr, na.rm =T), mp = min(df.ind$mc, na.rm = T), alphaV = alphaV, alphaP = alphaP)
out3 <- r*V*(1-(V/k))
out4 <- ((a*V*P) / (1 + (a*exp(h0 + alphaV*log(mv) + alphaP*log(mp))*V)))

plot(out ~ V, type = "n")
lines(out ~ V)
lines(out2 ~ V, col = "green")
lines(out3 ~ V, col = "blue")
lines(out4 ~ V, col = "brown")

plot(out4 ~ V)

r*V*(1-(V/k)) - ((a*V*P) / (1 + (a*exp(h0 + alphaV*log(mv) + alphaP*log(mp))*V)))



out3 <- dvdt(a = a, P = 1, h0 = h0, mv = mv, mp = mp, alphaV = alphaV, alphaP = alphaP)
out4 <- 

  
plot(out3 ~ V, type = "n")
lines(out3 ~ V)
lines(out2 ~ V, col = "green")

bd.FR <- function(n, a, p, t, mc, mr, r, k){
  
  h.log <- coef(mte1.h)[1] + coef(mte1.h)[2]*log(mc) + coef(mte1.h)[3]*log(mr)
  h <- exp(h.log)
  
  #predicted consumption rate = 
  cr <- a*n*p*t/(1+a*h*n)
  
  # urchin growth rate = 
  ur <- r*n*(1-(n/k))
  
  dudt <- ur - cr
  dudt
}

n = seq(-10^6, 10^6, length.out = 1000)

out1 <- bd.FR(n = n, a = mean(df.ind$median.a, na.rm = T), p = 0.16, mc = max(df.ind$mc, na.rm = T), mr = min(df.ind$mr, na.rm = T), t = 1, r = 0.1, k = 30)
out2 <- bd.FR(n = n, a = mean(df.ind$median.a, na.rm = T), p = 100, mc = min(df.ind$mc, na.rm = T), mr = max(df.ind$mr, na.rm = T), t = 1, r = 1.01, k = 10^4)
out3 <- bd.FR(n = n, a = mean(df.ind$median.a, na.rm = T), p = 100, mc = mean(df.ind$mc, na.rm = T), mr = mean(df.ind$mr, na.rm = T), t = 1, r = 1.01, k = 10^4)

plot(out1 ~ n, type = "n", ylim = c(-1000000, 1000000))
lines(out1 ~ n, col = "blue")
lines(out2 ~ n, col = "red")
lines(out3 ~ n, col = "black")
abline(a = 0, b = 0)

dudt <- vector()
input <- seq(-1000, 1000, by = 0.1)

for(i in 1:length(input)){
  dudt[i] <- bd.FR(n = input[i],a = mean(df.ind$median.a, na.rm = T), p = mean(lob.a$density, na.rm = T), mc = max(df.ind$mc, na.rm = T), mr = min(df.ind$mr, na.rm = T), t = 1, r = 0.1, k = 10)
  
}

temp <- data.frame(input, dudt)

temp[temp$dudt > -1*0.1 & temp$dudt < 0.1, ]


dudt[dudt > -1*0.1 & dudt < 0.1]



logistic <- function(n, k, r){
  
  # urchin growth rate = 
  dudt <- r*n*(1-(n/k))
  dudt
}

out <- vector()
input <- seq(-100, 100)
for(i in 1:length(input)){
  out[i] <- logistic(n = input[i], k = 30, r = 0.1 )
}

temp <- data.frame(input, out)

temp[temp$out > -1*0.1 & temp$out < 0.1, ]








# dvdt will be at equilibrium when p = r/a * (1 - V/k) * (1 + a * exp(h0 + alphaV*log(mv) + alphaP*log(mp))*V). Thus dvdt will go negative when p is < that value and be positive when p > that value. This equilibria depends only on p, mc, r, 

pred <- summary(lob.a$density)
pred <- c(0, 0.25, 1)
n = seq(0, 37, length.out = 1000)

mat <- matrix(ncol = length(pred), nrow = length(n))

for(i in 1:length(pred)){
temp <- bd.FR(n = n, a = mean(df.ind$median.a, na.rm = T), p = pred[i], mc = max(df.ind$mc, na.rm = T), mr = min(df.ind$mr, na.rm = T), t = 1, r = 0.01, k = 30)

mat[,i] <- temp
}


matplot(x = n, y = mat, type = "l")
abline(a = 0, b = 0)
abline(v = 30)


mc <- c(min(df.ind$mc, na.rm = T), mean(df.ind$mc,na.rm = T), max(df.ind$mc, na.rm = T))
mr <- c(max(df.ind$mr, na.rm = T), mean(df.ind$mr,na.rm = T), min(df.ind$mr, na.rm = T))


n = seq(0, 37, length.out = 1000)

mat <- matrix(ncol = length(mc), nrow = length(n))

for(i in 1:length(mc)){
  temp <- bd.FR(n = n, a = mean(df.ind$median.a, na.rm = T), p = 0.16, mc = mc[i], mr = mr[i], t = 1, r = 0.1, k = 30)
  
  mat[,i] <- temp
}


matplot(x = n, y = mat, type = "l")
abline(a = 0, b = 0)
abline(v = 30)


















coef(mte1.h)[2]*log(mc) + coef(mte1.h)[3]*log(mr)
ratio <- df.ind$mc^coef(mte1.h)[2] * df.ind$mr^coef(mte1.h)[3]


predict <- expand.grid(n = seq(0, max(urc.a$density), length.out = 1000), p = summary(lob.a$density))
predict$mc <- rep(seq(min(df.ind$mc, na.rm = T), max(df.ind$mc, na.rm = T), length.out = 1000), times = 6)
predict$mr <- rep(seq(max(df.ind$mr, na.rm = T), min(df.ind$mr, na.rm = T), length.out = 1000), times = 6)


predict <- expand.grid(n = seq(0, max(urc.a$density), length.out = 100), mc = seq(min(df.ind$mc^coef(mte1.h)[2], na.rm = T), max(df.ind$mc^coef(mte1.h)[2], na.rm = T), length.out = 10), mr = seq(max(df.ind$mr^coef(mte1.h)[3], na.rm = T), min(df.ind$mr^coef(mte1.h)[3], na.rm = T), length.out = 10))

# predict <- data.frame(n = rep(seq(0, max(urc.a$density), length.out = 10), times = 100), mc = rep(predict$mc, times = 10), mr = rep(predict$mr, times = 10))

predict$dudt <- bd.FR(n = predict$n, a = mean(df.ind$median.a, na.rm = T), p = 0.16, mc = predict$mc, mr = predict$mr, t = 1, r = 1.01, k = 30)

predict$ratio <- predict$mc * predict$mr

ggplot(predict, aes(x = n, y = ratio))+
  geom_raster(aes(fill = dudt))

plot(mc ~ jitter(mr), predict)
plot(ratio ~ jitter(n), predict)
plot(dudt ~ I(mc^coef(mte1.h)[2]*mr^coef(mte1.h)[3]), predict)


ratio <- df.ind$mc^coef(mte1.h)[2] * df.ind$mr^coef(mte1.h)[3]


bd.FR <- function(n, a, p, t, ratio, r, k){
  
  h <- coef(mte1.h)[1]*ratio
  
  #predicted consumption rate = 
  cr <- a*n*p*t/(1+a*h*n)
  
  # urchin growth rate = 
  ur <- r*n*(1-(n/k))
  
  dudt <- ur - cr
  dudt
}

predict <- expand.grid(n = seq(0, max(urc.a$density), length.out = 100), ratio=seq(min(ratio, na.rm = T), max(ratio, na.rm = T), length.out = 100))

predict$dudt <- bd.FR(n = predict$n, a = mean(df.ind$median.a, na.rm = T), p = 1, ratio = predict$ratio, t = 1, r = 0.1, k = 30)

ggplot(predict, aes(x = n, y = ratio))+
  geom_raster(aes(fill = dudt))+
  geom_contour(aes(z = dudt))










################


quad <- function(r = 0.1, a = mean(df.ind$median.a, na.rm = T), h0 = , mv, mp, alphaV, alphaP, V, P, k){
  h = h0*mv^alphaV*mp^alphaP
  a = 
  b = 
  c = 
    
  quad1 <- -1*b + sqrt(b^2 - )
  
  
  
  out <- (-1*r*a*h/k)*V^2 + (r*a*h - r/k)*V + (r - a*P)
  out
}












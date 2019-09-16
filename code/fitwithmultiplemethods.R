
###########################################################
## fit jags function
############################################################

fit.jags <- function(data, n.chains = 3, n.burnin = 10000, n.thin = 2, n.iter = 30000){
  killed <- data$killed
  initial <- data$initial
  
  jags.data = list("initial"= initial,
                   "killed" = killed,
                   "P" = 1, 
                   "T" = 1, 
                   n = length(initial), 
                   scale.a = scale.a, 
                   shape.a = shape.a, 
                   scale.h = scale.h, 
                   shape.h = shape.h) # named list
  
  model.loc=here("code","BH_jags.txt") # name of the txt file
  jags.params=c("a", "h")
  
  jags(jags.data,parameters.to.save=jags.params,inits=NULL,
       model.file=model.loc, n.chains = n.chains, n.burnin=n.burnin,
       n.thin=n.thin, n.iter=n.iter, DIC=TRUE)
  
}


##############################################################
## Individual level ML estimates using different packages
##############################################################

# This function will fit the holling type II using 4 different packages/approaches. DF must be in the form of names(df) <- c("killed", "initial", "ind"). When plot = F it returns a data.frame. When plot = T it returns plots (n = length(x$initial)) for each individual with curves for each fit.


fit.individual <- function(x, plot = F){
  
if(plot == F){
inds <- sort(unique(x$ind))
out <- list()
df.a <- data.frame(parameter = "a", ind = inds, mle2 = NA, frair = NA, nls = NA, JAGs = NA)
df.h <- data.frame(parameter = "h", ind = inds, mle2 = NA, frair = NA, nls = NA, JAGs = NA)

for(i in inds){
  tryCatch({
  # MLE2 package
  temp <- mle2(killed~dbinom(prob=pmax(eps,
                             pmin(1-eps,1/(1/a+h*initial))),
                   size=initial),start=list(a=.01,h=.001),data=x[x$ind == i,])
  
  df.a[i,3] <- coef(temp)[1]
  df.h[i, 3] <- coef(temp)[2]
  
  # Friar package
  
  temp2 <- frair_fit(killed~initial,data=x[x$ind == i, ],
                  response='hollingsII',
                  start=list(a = 0.01, h = 0.001), fixed=list(T = 1))
  
  df.a[i,4] <- coef(temp2)[1]
  df.h[i, 4] <- coef(temp2)[2]
  
  # nls
  
  temp3 <- nls(
    killed ~ (a*initial)/(1+a*initial*h),
    start = c(a = 0.01, h = 0.001), data = x[x$ind == i, ])
  
  df.a[i,5] <- coef(temp3)[1]
  df.h[i, 5] <- coef(temp3)[2]
  
  # jags
  
  temp4 <- fit.jags(data = x[x$ind == i, ])
  
  df.a[i,6] <- temp4$BUGSoutput$mean[1]
  df.h[i, 6] <- temp4$BUGSoutput$mean[3]

  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

df <- rbind(df.a, df.h) %>% gather(package,estimate, -c(parameter, ind)) %>% spread(parameter, estimate)
df$packag.id <- as.numeric(as.factor(df$package))
df
} else{

# make the plot
col <- c("red", "green", "blue", "gray", "black")

for(i in 1:length(inds)){
  plot(killed ~ jitter(initial),data = x[as.numeric(x$ind) == i, ], xlab="Number of prey",ylab="Number consumed", ylim = c(0,15))
  curve(holling2(x,df$a[df$ind == i & df$package == "frair"],df$h[df$ind == i & df$package == "frair"],P=1,T=1),add=TRUE,col = col[1],lty=1) #true curve
  curve(holling2(x,df$a[df$ind == i & df$package == "JAGs"],df$h[df$ind == i & df$package == "JAGs"],P=1,T=1),add=TRUE,col=col[2],lty=2) #true curve
  curve(holling2(x,df$a[df$ind == i & df$package== "mle2"],df$h[df$ind == i & df$package== "mle2"],P=1,T=1),add=TRUE,col=col[3],lty=3) #true curve
  curve(holling2(x,df$a[df$ind == i & df$package== "nls"],df$h[df$ind == i & df$package== "nls"],P=1,T=1),add=TRUE,col=col[4],lty=4) #true curve
  
  legend("topleft", legend = c("frair", "jags", "mle2", "nls"), col = c("red", "green", "blue", "gray"), lty = c(1,2,3,4))
}

}}



########################################################
## read in the lobster data 
#########################################################

df <- read.table(here("data/cleaned","loburc_cleaned.csv"), header = T, sep = ",")

x <- data.frame(killed = df$num_consumed, initial = df$num_offered, ind = as.numeric(df$id))

out <- fit.individual(x = x, plot = F)



df <- fit.individual(x, plot = F) 
fit.individual(x , plot = T)




temp <- df %>% dplyr::select(-packag.id) %>%
  gather(parameter, estimate, -c(ind, package)) %>%
  
  mutate(estimate = round(estimate, 3)) %>% 
           spread(package, estimate)





















































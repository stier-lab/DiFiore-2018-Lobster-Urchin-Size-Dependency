#--------------------------------------------------------------------------------------------
## Fit flexible functional resonse w/ prey depletion
#--------------------------------------------------------------------------------------------


# Setup
library(here)
source(here("code", "setup.R"))
source(here("code", "functions.R"))
library(frair)

df <- read.table(here("data/cleaned","loburc_cleaned.csv"), header = T, sep = ",") %>%
  mutate(lobcat = ifelse(size <=70, "small", 
                         ifelse(size >70 & size <=90, "medium", 
                                ifelse(size > 90, "large", NA))), 
         newtreat = as.factor(paste(treatment, lobcat, sep = "-"))) %>%
  arrange(newtreat, id)


# Using frair

# this model assumes that the FR can range from a type II to a type III, and accounts for the fact that prey were not replaced. I fit for both experimental time and in units of hours. 
m1 <- frair_fit(killed~initial, data=df[df$newtreat == "urc_large-small", ], response='flexp', 
                   start=list(b = 0.4, h = 0.06, q = 0.), fixed=list(T=1))
plot(m1)
lines(m1)
summary(m1$fit)

m1.hours <- frair_fit(killed~initial, data=df[df$newtreat == "urc_small-large", ], response='flexp', 
                start=list(b = 0.002, h = 5, q = 1.2), fixed=list(T=48))
plot(m1.hours)
lines(m1.hours)
summary(m1.hours$fit)

m2 <- frair_fit(killed~initial, data=df[df$treatment == "urc_medium", ], response='flexp', 
                start=list(b = 0.4, h = 0.12, q = 0.89), fixed=list(T=1))
plot(m2)
lines(m2)

m2.hours <- frair_fit(killed~initial, data=df[df$treatment == "urc_medium", ], response='flexp', 
                start=list(b = 0.3/48, h = 0.3*48, q = 0.6), fixed=list(T=48))
plot(m2.hours)
lines(m2.hours)

m3 <- frair_fit(killed~initial, data=df[df$treatment == "urc_large", ], response='flexp', 
                start=list(b = 0.4, h = 0.12, q = 0.89), fixed=list(T=1))
plot(m3)
lines(m3)

m3.hours <- frair_fit(killed~initial, data=df[df$treatment == "urc_large" & df$initial > 1, ], response='flexp', 
                start=list(b = 0.01/48, h = 1.1*48, q = 3), fixed=list(T=48))
plot(m3.hours)
lines(m3.hours)

# Ok so fit it across all individuals... see if we can get it to go using mle2

fit.individual <- function(x, plot = F, T = 48){
  
  if(plot == F){
    inds <- sort(unique(x$ind))
    manual_fit <- list()
    df.b <- data.frame(parameter = "b", ind = inds, mle2 = NA)
    df.q <- data.frame(parameter = "q", ind = inds, mle2 = NA)
    df.h <- data.frame(parameter = "h", ind = inds, mle2 = NA)
    
    for(i in inds){
      tryCatch({
        # MLE2 package
        cntrl <- list(trace = 3, maxit = 1000)
        manual_fit[[i]] <- coef(mle2(flexp_nll, start=list(b = 0.01, h = 1, q = 0.8), fixed=list(T = T),
                           method='SANN', data=data.frame(X = x$initial[x$ind == inds[i]], Y = x$killed[x$ind == inds[i]])))
        
        temp <- mle2(flexp_nll, start=list(b = manual_fit[[i]][1], q = manual_fit[[i]][2], h = manual_fit[[i]][3]), fixed=list(T = T), 
                     method='BFGS', data=data.frame(X = x$initial[x$ind == inds[i]], Y = x$killed[x$ind == inds[i]]), control = cntrl)
        
        df.b[i,3] <- coef(temp)[1]
        df.q[i, 3] <- coef(temp)[2]
        df.h[i, 3] <- coef(temp)[3]
        
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    }
    
    df <- rbind(df.b, df.q, df.h) %>% gather(package,estimate, -c(parameter, ind)) %>% spread(parameter, estimate)
    df$packag.id <- as.numeric(as.factor(df$package))
    df
  } else{
    
    # make the plot
    col <- c("red", "green", "blue", "gray", "black")
    
    for(i in 1:length(inds)){
      plot(killed ~ jitter(initial),data = x[as.numeric(x$ind) == i, ], xlab="Number of prey",ylab="Number consumed", ylim = c(0,15))
      curve(holling2(x,df$a[df$ind == i & df$package == "mle2"],df$h[df$ind == i & df$package == "frair"],P=1,T=1),add=TRUE,col = col[1],lty=1) #true curve

    }
    
  }}

x <- data.frame(killed = df$killed, initial = df$initial, ind = as.numeric(df$id))

out <- fit.individual(x = x, plot = F)



# The role of T
## Users need to be aware that changing T will change 
## the units of fitted coefficients.  
## For example, with the Gammarus dataset:
g_T1 <- frair_fit(formula = eaten~density, data = gammarus, 
                  response = "rogersII", 
                  start = list(a = 2, h = 0.1), fixed = list(T = 1))
g_Td <- frair_fit(formula = eaten~density, data = gammarus, 
                  response = "rogersII", 
                  start = list(a = 1, h = 0.1), fixed = list(T = 40/24))
g_Th <- frair_fit(formula = eaten~density, data = gammarus, 
                  response = "rogersII", 
                  start = list(a = 0.05, h = 4), fixed = list(T = 40))
diff_t <- round(rbind(coef(g_T1), coef(g_Td), coef(g_Th)), 2)
row.names(diff_t) <- c("g_T1 (Experimental Time)", "g_Td (Days)", "g_Th (Hours)")
print(diff_t)


x <- data.frame(killed = df$killed, initial = df$initial, treat = as.numeric(df$newtreat))

fit.treats <- function(x, plot = F, T = 48){
  
  if(plot == F){
    treats <- sort(unique(x$treat))
    manual_fit <- list()
    df.b <- data.frame(parameter = "b", treat = treats, mle2 = NA)
    df.q <- data.frame(parameter = "q", treat = treats, mle2 = NA)
    df.h <- data.frame(parameter = "h", treat = treats, mle2 = NA)
    
    for(i in treats){
      tryCatch({
        # MLE2 package
        cntrl <- list(trace = 3, maxit = 1000)
        manual_fit[[i]] <- coef(mle2(flexp_nll, start=list(b = 0.01, h = 1, q = 0.8), fixed=list(T = T),
                                     method='SANN', data=data.frame(X = x$initial[x$treat == treats[i]], Y = x$killed[x$treat == treats[i]])))
        
        temp <- mle2(flexp_nll, start=list(b = manual_fit[[i]][1], q = manual_fit[[i]][2], h = manual_fit[[i]][3]), fixed=list(T = T), 
                     method='BFGS', data=data.frame(X = x$initial[x$treat == treats[i]], Y = x$killed[x$treat== treats[i]]), control = cntrl)
        
        df.b[i,3] <- coef(temp)[1]
        df.q[i, 3] <- coef(temp)[2]
        df.h[i, 3] <- coef(temp)[3]
        
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    }
    
    df <- rbind(df.b, df.q, df.h) %>% gather(package,estimate, -c(parameter, treat)) %>% spread(parameter, estimate)
    df$packag.id <- as.numeric(as.factor(df$package))
    df
  } else{
    
    # make the plot
    col <- c("red", "green", "blue", "gray", "black")
    
    for(i in 1:length(inds)){
      plot(killed ~ jitter(initial),data = x[as.numeric(x$ind) == i, ], xlab="Number of prey",ylab="Number consumed", ylim = c(0,15))
      curve(holling2(x,df$a[df$ind == i & df$package == "mle2"],df$h[df$ind == i & df$package == "frair"],P=1,T=1),add=TRUE,col = col[1],lty=1) #true curve
      
    }
    
  }}

out <- fit.treats(x = x, plot = F)

























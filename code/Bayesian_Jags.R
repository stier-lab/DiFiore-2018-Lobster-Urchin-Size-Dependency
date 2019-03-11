
####################################
## Source prior code 
####################################

source("code/setup.R")
source("code/functions.R")


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
                 a ~ dunif(0.001, 2)
                 h ~ dunif(0.001, 2)
                 
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

fit.jags <- function(pred_size, urc_size, n.chains = 3, n.burnin = 10000, n.thin = 2, n.iter = 20000){
  temp <- df %>% filter(class == pred_size, treatment == urc_size)
  killed <- temp$num_consumed
  initial <- temp$num_offered
  
  jags.data = list("initial"= initial,
                   "killed" = killed,
                   "P" = 1, 
                   "T" = 1, 
                   n = length(initial)) # named list
  
  jags(jags.data,parameters.to.save=jags.params,inits=NULL,
               model.file=model.loc, n.chains = n.chains, n.burnin=n.burnin,
               n.thin=n.thin, n.iter=n.iter, DIC=TRUE)
}



# Fit to each size/treatment combination and plot.

d <- par(mfrow = c(3,4))
mod1 <- fit.jags("jumbo", "urc_large")
plot.jags(mod1, "jumbo", "urc_large")
mod2 <- fit.jags("large", "urc_large")
plot.jags(mod2, "large", "urc_large")
mod3 <- fit.jags("medium", "urc_large")
plot.jags(mod3, "medium", "urc_large")
mod4 <- fit.jags("small", "urc_large")
plot.jags(mod4, "small", "urc_large")
mod5 <- fit.jags("jumbo", "urc_medium")
plot.jags(mod5,"jumbo", "urc_medium")
mod6 <- fit.jags("large", "urc_medium")
plot.jags(mod6,"large", "urc_medium")
mod7 <- fit.jags("medium", "urc_medium")
plot.jags(mod7,"medium", "urc_medium")
mod8 <- fit.jags("small", "urc_medium")
plot.jags(mod8,"small", "urc_medium")
mod9 <- fit.jags("jumbo", "urc_small")
plot.jags(mod9,"jumbo", "urc_small")
mod10 <- fit.jags("large", "urc_small")
plot.jags(mod10,"large", "urc_small")
mod11 <- fit.jags("medium", "urc_small")
plot.jags(mod11,"medium", "urc_small")
mod12 <- fit.jags("small", "urc_small")
plot.jags(mod12,"small", "urc_small")
par(d)

# Make plots of parameter fits

mods <- paste("mod", seq(1,12), sep = "") #vector of model names

mod.coefs <- lapply(mods, function(x) {
  MCMCsummary(get(x), params = c("a", "h"))
}) #apply the MCMCsummary() function to extract summaries of a and h for each model

names(mod.coefs) <- mods #name the list 

mod.coef <- list_df2df(mod.coefs) #convert the list to a data.frame using the nifty qdap package function

names(mod.coef)[c(1,4,5,6)] <- c("model", "c2.5", "c50", "c97.5") #fix the name
mod.coef$parameter <- rep(c("a", "h"), times = 12) # add variable for parameter


# make correct labels for each model
treat <- data.frame(model = mods, lob.size = rep(c("jumbo", "large", "medium", "small"), times = 3), urc.size = rep(c("large", "medium", "small"), each = 4))
treat$combo <- paste(treat$lob.size, treat$urc.size, sep = "_")

mod.coef <- mod.coef %>% left_join(treat)

pred.labels <- rev(rep(c("small", "medium", "large", "jumbo"), each = 3))
prey.labels <- rev(rep(c("small", "medium", "large"), times = 4))

# plot coefficients and 95% confidence intervals
ggplot(mod.coef)+
  geom_linerange(aes(x = as.numeric(as.factor(combo)), ymin = c2.5, ymax = c97.5), lwd = 1) +
  geom_point(aes(x = as.numeric(as.factor(combo)), y = c50), size = 2, lwd = 1, shape = 21)+
  scale_x_continuous(breaks = 1:length(pred.labels), labels = pred.labels, 
                     sec.axis = sec_axis(~., breaks = 1:length(pred.labels), 
                                         labels = prey.labels)) +
  coord_flip() +
  facet_wrap(~parameter, ncol = 1)+
  theme_bw()+
  theme(axis.title.y = element_blank(), axis.text = element_text(size = 12), axis.title = element_text(size = 14))



####################################
## Build Violin Plot
####################################

#build df to use in plot
out <- list()

for(i in mods){
  temp <- MCMCchains(get(i))[, c("a", "h")]
  out[[i]] <- temp
} #extract information for chains of each model

chains <- list_df2df(out)
names(chains)[1] <- "model"

chains <- chains %>% gather(parameter, estimate, -model) %>%
  left_join(treat)

ggplot(chains)+
  geom_violin(aes(x = combo, y = estimate))+
  scale_x_discrete(labels = pred.labels)+
  sec_axis(~ . * 1, breaks = 1:length(prey.labels), labels = prey.labels)+
  #geom_boxplot(aes(x = combo, y = estimate), width=0.1)+
  coord_flip() +
  facet_wrap(~parameter, ncol = 1)+
  theme_bw()+
  theme(axis.title.y = element_blank(), axis.text = element_text(size = 12), axis.title = element_text(size = 14))











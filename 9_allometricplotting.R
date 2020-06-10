#---------------------------------------------------------------------
## Get data
#--------------------------------------------------------------------

# Posteriors
df.ind <- read.csv(here::here("data/cleaned/posteriors", "allometric_individual.csv")) %>% as_tibble()
df.pop <- read.csv(here::here("data/cleaned/posteriors", "allometric_population.csv")) %>% as_tibble()

# Experimental data
df <- read.table(here("data/cleaned","loburc_cleaned.csv"), header = T, sep = ",") 




# Function to estimate mean and CI's of posteriors

# assumes no variance between individuals...
p_link <- function(N, mc, mr, data = df.pop){
  loga <- with(data, mu.alpha.a + beta1.a*log(mc) + beta2.a*log(mr))
  logh <- with(data, mu.alpha.h + beta1.h*log(mc) + beta2.h*log(mr))
  a <- exp(loga)
  h <- exp(logh)
  
  a*N*48 / (1 + a*h*N)
}


# Posterior predictive simulation
a_sim <- with(df.pop, rnorm( length(df.pop$mu.alpha.a) , mu.alpha.a , 1/sqrt(var.a.ind)))
h_sim <- with(df.pop, rnorm(length(df.pop$mu.alpha.h), mu.alpha.h, 1/sqrt(var.h.ind)))

p_sim <- function(N, mc, mr, data = df.pop){
  loga <- with(data, a_sim + beta1.a*log(mc) + beta2.a*log(mr))
  logh <- with(data, h_sim + beta1.h*log(mc) + beta2.h*log(mr))
  a <- exp(loga)
  h <- exp(logh)
  
  a*N*48 / (1 + a*h*N)
  # prob <- 48/(1/a + h*N)
  # killed <- rbinom(length(prob), size = N, prob = prob)
  # return(killed)
  
}

allometric_CI <- function(mc, mr, prob = 0.95, ...){
  temp.mr <- as.numeric(mr)
  temp.mc <- as.numeric(mc)

N.vec <- seq(0, 26, length.out = 100)
p_sim_output <- sapply(N.vec, function(i) p_sim(i, mc = temp.mc, mr = temp.mr))
p_mu <- apply(p_sim_output, 2 ,mean, na.rm = T)
p_ci <- t(apply( p_sim_output , 2 , PI, prob = prob))

return(data.frame(N = N.vec, mu = p_mu, mu.lower = p_ci[,1], mu.upper = p_ci[,2]))
}

predict <- expand.grid(mc = c(median(df$mc, na.rm = T), quantile(df$mc, probs = c(0.1, 0.9))), mr = unique(df$mr))
#predict$treatment <- rep(c("Medium urchins", "Large urchins", "Small urchins"), each = 3)

forplot <- predict %>%
  purrr::pmap(allometric_CI) %>%
  set_names(paste(predict$mc, predict$mr, sep = "-")) %>%
  do.call(rbind, .) %>%
  mutate(mc = rep(predict$mc, each = 100),
         mr = rep(predict$mr, each = 100)) %>%
  rename(initial = N)

forplot$treatment <- rep(c("urc_medium", "urc_large", "urc_small"), each = 300)
forplot <- arrange(forplot, mc)
forplot$lob.sizeclass <- rep(rev(c("Large", "Medium", "Small")), each = 300)



plot4 <- ggplot(df, aes(x = initial, y = killed))+
  geom_jitter(aes(size = treatment), pch = 21, show.legend = F, stroke = 1)+
  scale_size_manual(values = c(4,2.5,1.5))+
  #scale_fill_manual(values = c('#d53e4f','#fc8d59','#fee08b'))+
  geom_ribbon(data = forplot, aes(ymin = mu.lower, ymax = mu.upper, y = mu, group = lob.sizeclass), color = "gray", alpha = 0.25)+
  geom_line(data = forplot, aes(x = initial, y = mu, color = lob.sizeclass), size = 1.5)+
  scale_color_manual(values = rev(c('#d53e4f','#fc8d59','#fee08b')))+
  facet_wrap(~treatment)+
  geom_hline(yintercept = 0, lty = 4, color = "gray25")+
  labs(x = "Number of prey offered", y = "Number of prey consumed", color = "")

ggsave(here::here("figures/", "allometric_fr.png"), plot4, width = 10, height = 4)



#-------------------------------------------------------------------------
## Posterior vs. prior plots
#-------------------------------------------------------------------------

tau.jags <- function(tau){
  var = 1/tau
  sd = sqrt(var)
  return(c(var, sd))
}



prior.beta1.h <- rnorm(length(df.pop$beta1.h), mean = -0.75, sd = tau.jags(rgamma(length(df.pop$beta1.h), shape = 5, rate = 1)))
plot(density(prior.beta1.h))


prior.beta2.h <- rnorm(length(df.pop$beta1.h), mean = 0.75, sd = tau.jags(rgamma(length(df.pop$beta1.h), shape = 5, rate = 1)))
plot(density(prior.beta2.h))




df.plot <- data.frame(prior.beta1.h, posterior.beta1.h = df.pop$beta1.h, 
                   prior.beta2.h, posterior.beta2.h = df.pop$beta2.h) %>%
  gather(dist, value) %>%
  separate(dist, into = c("model", "scaling_exponent", "parameter"), sep = "[.]") %>%
  mutate(parameter = paste(scaling_exponent, parameter, sep = "."))

meta <- df.plot %>% group_by(model, parameter) %>% 
  summarize(mean = mean(value))
meta$text <- c(NA, NA, "MTE\n-0.75", "MTE\n0.75")

p5 <- ggplot(df.plot, aes(x = value, fill = model))+
  geom_density(alpha = 0.5)+
  facet_wrap(~parameter)+
  coord_cartesian(xlim = c(-3, 2))+
  geom_vline(data = meta, aes(xintercept = mean), lty = 4, color = "gray30")+
  geom_text(data = meta, aes(x = c(0, 0, -0.3, 0.3), y = 4, label = text))+
  labs(x = "Allometric scaling exponent", y = "Density", fill = "")+
  theme(legend.position = c(0.5,0.5))

ggsave(here::here("figures/", "posterior-prior.png"), p5, width = 8, height = 3)





#------------------------------------------------------------------------
## Scrap
#------------------------------------------------------------------------


temp <- p_sim(N = 2, mc = mean(df$mc, na.rm = T), mr = mean(df$mr, na.rm = T))











temp <- matrix(c(1, 2, 3, NA, 5, 6), nrow = 3, ncol = 2)
apply(temp, 2 ,mean, na.rm = TRUE)



temp <- p_sim(N = 26, mc = mean(df$mc, na.rm = T), mr = mean(df$mr, na.rm = T))

killed <- rbinom(length(temp), size = 10, prob = temp)
p_mu <- apply( killed , 2 , mean, na.rm = T )















allometric_CI <- function(mc, mr, data, prob = 0.95){
  temp.mr <- as.numeric(mr)
  temp.mc <- as.numeric(mc)
  mu.link <- function(mc, mr = temp.mr, alpha = data[, "alpha"], beta1 = data[, "beta1"], beta2 = data[, "beta2"]){
    exp(alpha)*mc^beta1*mr^beta2
  } # defines a function to predict the prey killed at combination of a and h in the posteriors
  
  mc.seq <- seq( from=min(df$mc) , to=max(df$mc) , length.out = 100) #define a sequence of consumer masses
  mu <- sapply( mc.seq, mu.link) # apply the mu.link funciton to each N in the sequence
  
  mu.median <- apply( mu , 2 , median ) # calculate the median predicted value for each N
  mu.PI <- t(apply( mu , 2 , PI , prob=prob )) # calculate the credible interval for each value of N
  
  return(data.frame(mc.seq = mc.seq, mu = mu.median, mu.lower = mu.PI[,1], mu.upper = mu.PI[,2]))
}


calculate.CI <- function(a.parameter, h.parameter, prob = 0.95){
  require(rethinking)
  
  mu.link <- function(N, a, h){
    df.model[,a.parameter]*N/(1+df.model[,a.parameter]*df.model[,h.parameter]*N)
  } # defines a function to predict the prey killed at combination of a and h in the posteriors
  
  N.seq <- seq( from=0 , to=60 , length.out = 100 ) # define a sequence of initial densities
  mu <- sapply( N.seq , mu.link) # apply the mu.link funciton to each N in the sequence
  
  mu.median <- apply( mu , 2 , median ) # calculate the median predicted value for each N
  mu.PI <- t(apply( mu , 2 , PI , prob=prob )) # calculate the credible interval for each value of N
  
  return(data.frame(N.seq = N.seq, mu = mu.median, mu.lower = mu.PI[,1], mu.upper = mu.PI[,2])) # return a data.frame to organize output!
  
}





ggplot(df, aes(x = initial, y = killed/48))+
  geom_jitter(aes(size = treatment, fill = treatment), pch = 21, show.legend = F, stroke = 1)+
  scale_size_manual(values = c(4,2.5,1.5))+
  scale_fill_manual(values = c('#d53e4f','#fc8d59','#fee08b'))+
  #geom_line(data = newdat, aes(x = initial, y = killed, color = lobcat), size = 1.5, show.legend = F)+
  #scale_color_manual(values = c('#d53e4f','#fc8d59','#fee08b'))+
  facet_wrap(~treatment)

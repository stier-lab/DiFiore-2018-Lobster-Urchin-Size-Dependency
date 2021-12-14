library(tidyverse)

set.seed(10001)
#set.seed(10)

# Build example code to ask for help 

nsites <- 10
ndraws <- 100 # number of draws from the body size distribution of predators and prey

# estimate size distributions for each site, where the mean size differs between sites

pred.size <- matrix(ncol = nsites, nrow = ndraws)
prey.size <- matrix(ncol = nsites, nrow = ndraws)

# To emulate the behavior in my data, I allow the mean body size to differ between sites and crank down the within site variation in body size to almost zero. 

mean.pred.size <- runif(n = nsites, min = 4, max = 7)
mean.prey.size <- runif(n = nsites, min = 0.5, max = 1.5)

for(i in 1:nsites){
  pred.size[,i] <- rlnorm(n = ndraws, meanlog = mean.pred.size[i], sd = 1)
  prey.size[,i] <- rlnorm(n = ndraws, meanlog = mean.prey.size[i], sd = 1)
}
# 
# pred.density <- rnorm(n = nsites, mean = 1, sd = 0.5)^2
# prey.density <- rnorm(n = nsites, mean = 5, sd = 0.5)^2

# If the technique works that I should be able to eliminate variation in body size in the "real" data, crank up between site variation in density, and then show that most of the variation is due to density not body size.
# for(i in 1:nsites){
#   pred.size[,i] <- rlnorm(n = ndraws, mean = runif(n = nsites, min = 4, max = 4), sd = 1)
#   prey.size[,i] <- rlnorm(n = ndraws, mean = runif(n = nsites, min = 0.5, max = 0.5), sd = 1)
# }

# estimate the density of predators and prey for each site, where

pred.density <- rnorm(n = nsites, mean = 0.1, sd = 1)^2
prey.density <- rnorm(n = nsites, mean = 0.1, sd = 1)^2


# Build some empty vectors 
out <- matrix(nrow = ndraws, ncol = nsites)
mat.a <- matrix(nrow = ndraws, ncol = nsites)
mat.h <- matrix(nrow = ndraws, ncol = nsites)

# Parameter values
a0. = 1
h0. = 1
beta1a. = 0.75
beta2a. = 0.75
beta1h. = -0.75
beta2h. = 0.5


# Simulate interaction strengths for each site
  # This equation is the same that I use in the paper. Here I have set the parameter values to parameter values expected from theory for simplicity.
for(i in 1:nsites){
  for(j in 1:ndraws){
    mat.a[j,i] <- exp(a0. + beta1a.*log(pred.size[j,i]) + beta2a.*log(prey.size[j,i]))
    mat.h[j,i] <- exp(h0. + beta1h.*log(pred.size[j,i]) + beta2h.*log(prey.size[j,i]))

    out[j,i] <- mat.a[j,i]*pred.density[i]*prey.density[i] / (1 + mat.a[j,i]*mat.h[j,i]*prey.density[i])
  }
}

# Clean up the data a bit
full <- data.frame(out)
names(full) <- paste("site", 1:10, sep = "_")
full <- full %>% pivot_longer(cols = site_1:site_10, names_to = "site", values_to = "interaction_strenth")

hist(log(full$interaction_strenth))

# So my main question is: 
  # HOW WOULD YOU PARTITION VARIATION DUE TO BODY SIZE (i.e. the combined variation due to both variation in predator size and prey size) AND DENSITY? 


#-------------------------------------------------------------------
## Strategy two
#-------------------------------------------------------------------
  
# The second strategy that I came up with is based on regression. Linear regression measures the variation in the data explained by the predictors. So at least in theory, if I set the simulation that includes both sources of variation as my response variable and then predict the full model using the results of a simulation that isolates variation due to body size then I should be able to use the R2 as a measure of variance explained by body size/density. 

    # Estimate the interactions if density was the same across all sites
    
    out <- matrix(nrow = ndraws, ncol = nsites)
    mat.a <- matrix(nrow = ndraws, ncol = nsites)
    mat.h <- matrix(nrow = ndraws, ncol = nsites)
    for(i in 1:nsites){
      for(j in 1:ndraws){
        mat.a[j,i] <- exp(a0. + beta1a.*log(pred.size[j,i]) + beta2a.*log(prey.size[j,i]))
        mat.h[j,i] <- exp(h0. + beta1h.*log(pred.size[j,i]) + beta2h.*log(prey.size[j,i]))
        
        out[j,i] <- mat.a[j,i]*max(pred.density)*min(prey.density) / (1 + mat.a[j,i]*mat.h[j,i]*min(prey.density))
      }
    }
    
    duetobodysize <- data.frame(out)
    names(duetobodysize) <- paste("site", 1:10, sep = "_")
    duetobodysize <- duetobodysize %>% pivot_longer(cols = site_1:site_10, names_to = "site", values_to = "interaction_strenth")
    
    # Estimate the interactions if body size was the same across all sites, as if every individual was the same size, and that size was the regional average.
    
      # Build a matrix of individuals by site that are all the same average size
    
    pred.size.mean <- matrix(ncol = nsites, nrow = ndraws)
    prey.size.mean <- matrix(ncol = nsites, nrow = ndraws)

    for(i in 1:nsites){
      pred.size.mean[,i] <- rep(mean(pred.size), ndraws)
      prey.size.mean[,i] <- rep(mean(prey.size), ndraws)
    }
    
    out <- matrix(nrow = ndraws, ncol = nsites)
    mat.a <- matrix(nrow = ndraws, ncol = nsites)
    mat.h <- matrix(nrow = ndraws, ncol = nsites)
    for(i in 1:nsites){
      for(j in 1:ndraws){
        mat.a[j,i] <- exp(a0. + beta1a.*log(pred.size.mean[j,i]) + beta2a.*log(prey.size.mean[j,i]))
        mat.h[j,i] <- exp(h0. + beta1h.*log(pred.size.mean[j,i]) + beta2h.*log(prey.size.mean[j,i]))
        
        out[j,i] <- mat.a[j,i]*pred.density[i]*prey.density[i] / (1 + mat.a[j,i]*mat.h[j,i]*prey.density[i])
      }
    }
    
    duetodensity <- data.frame(out)
    names(duetodensity) <- paste("site", 1:10, sep = "_")
    duetodensity <- duetodensity %>% pivot_longer(cols = site_1:site_10, names_to = "site", values_to = "interaction_strenth")
    
    # Partition the variance using linear regression techniques based on the different simulations.

df <- data.frame(site = full$site, full = full$interaction_strenth, duetobodysize = duetobodysize$interaction_strenth, duetodensity = duetodensity$interaction_strenth)

ggplot(df, aes(y = full, x = duetobodysize))+
  geom_point(aes(color = site))

ggplot(df, aes(y = full, x = duetodensity))+
  geom_point(aes(color = site))

summary(lm(full ~ duetobodysize, df))
summary(lm(full~ duetodensity, df))

summary(lm(full~duetobodysize * duetodensity, df))



lm.full <- lm(full~duetobodysize * duetodensity, df)
lm.1 <- lm(full~duetobodysize, df)
lm.2 <- lm(full~duetodensity, df)
anova(lm.full, lm.1, lm.2)
anova(lm.full)

summary(lm.full)$r.squared
summary(lm.1)$r.squared
summary(lm.2)$r.squared

temp <- vegan::varpart(df$full, df$duetobodysize, df$duetodensity, df$duetobodysize * df$duetodensity)


plot(temp)











#block




























#-------------------------------------------------------------------
## Strategy one
#-------------------------------------------------------------------

# The first strategy that I have used is essential to estimate the variance in the full simulation and compare that variance to the variance from reduced simulations that fix body size and density respectively at the mean across sites. My idea here is that there can only be variation due to body size or denisty in these simulations therefore variance due to body size + variance due to density = total variance. For example: 

# First estimate the interactions if body size was the same across all sites
mat.a <- vector()
mat.h <- vector()
out <- vector()

# Fix the body size to the mean across sites
for(i in 1:nsites){
  mat.a[i] <- exp(a0. + beta1a.*log(mean(pred.size[i])) + beta2a.*log(mean(prey.size[i])))
  mat.h[i] <- exp(h0. + beta1h.*log(mean(pred.size[i])) + beta2h.*log(mean(prey.size[i])))
  
  out[i] <- mat.a[i]*pred.density[i]*prey.density[i] / (1 + mat.a[i]*mat.h[i]*prey.density[i])
}

duetodensity <- out
# Then estimate the interactions if density was the same across all sites

out <- matrix(nrow = ndraws, ncol = nsites)
mat.a <- matrix(nrow = ndraws, ncol = nsites)
mat.h <- matrix(nrow = ndraws, ncol = nsites)
for(i in 1:nsites){
  for(j in 1:ndraws){
    mat.a[j,i] <- exp(a0. + beta1a.*log(pred.size[j,i]) + beta2a.*log(prey.size[j,i]))
    mat.h[j,i] <- exp(h0. + beta1h.*log(pred.size[j,i]) + beta2h.*log(prey.size[j,i]))
    
    out[j,i] <- mat.a[j,i]*mean(pred.density)*mean(prey.density) / (1 + mat.a[j,i]*mat.h[j,i]*mean(prey.density))
  }
}

duetobodysize <- data.frame(out)
names(duetobodysize) <- paste("site", 1:10, sep = "_")
duetobodysize <- duetobodysize %>% pivot_longer(cols = site_1:site_10, names_to = "site", values_to = "interaction_strenth")

# Then partition the variance.

var_duetodensity = var(duetodensity)
var_duetobodysize = var(duetobodysize$interaction_strenth)
var_total = var(full$interaction_strenth)

var_duetodensity / var_total
var_duetobodysize / var_total # Obviously, something isn't correct with this method. How can variance increase, if I have reduced the sources of variation going into the simulation? Similarly, if you remove the set.seed at line 3 in this code, you will get a totally different answer with each run of the code.

# Check it
var_duetodensity / var_total + var_duetobodysize / var_total


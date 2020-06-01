# ----------------------------------------------------
## Equilibrium analyis of experimental resutls
#-----------------------------------------------------

# If we assume that urchin poplulations grow a constant rate to some carrying capacity, and predators reduce urchin populations according to a size-dependent type II funcitonal repsonse, we can predict how variation in lobster and urchin size changes how strongly these species interact at equilibrium. Specifically, one theoretical way of assessing the interaction strength of interaction species is to determine the chance in equilibrium abundance with and without the predator present (Wootton and Emmerson 2004, DeLong et al. 2015). Therefore, we estimated the equilibirum biomass of urchins with predators present relative to when predators were absent.


# write a function to estimate the rate of change in the urchin population
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



# Allow predator density to vary at the largest asymetry in body sizes.

        pred <- summary(lob.a$density)[c(1,4,5)]
        #pred <- c(0, 0.25, 1)
        n = seq(0, 37, length.out = 1000)
        
        mat <- matrix(ncol = length(pred), nrow = length(n))
        
        for(i in 1:length(pred)){
          temp <- bd.FR(n = n, a = mean(df.ind$median.a, na.rm = T), p = pred[i], mc = max(df.ind$mc, na.rm = T), mr = min(df.ind$mr, na.rm = T), t = 1, r = 0.1, k = 30)
          
          mat[,i] <- temp
        }
        
        
        matplot(x = n, y = mat, type = "l", xlab = "Urchin density", ylab = "dU/dt", main = "Variation in predator density")
        abline(a = 0, b = 0, lty = 3)
        abline(v = 30, lty = 3)



# Allow relative body sizes to vary
        mc <- c(min(df.ind$mc, na.rm = T), mean(df.ind$mc,na.rm = T), max(df.ind$mc, na.rm = T))
        mr <- c(max(df.ind$mr, na.rm = T), mean(df.ind$mr,na.rm = T), min(df.ind$mr, na.rm = T))
        
        n = seq(0, 37, length.out = 1000)
        
        mat2 <- matrix(ncol = length(mc), nrow = length(n))
        
        for(i in 1:length(mc)){
          temp <- bd.FR(n = n, a = mean(df.ind$median.a, na.rm = T), p = 0.16, mc = mc[i], mr = mr[i], t = 1, r = 0.1, k = 30)
          
          mat2[,i] <- temp
        }
        
        
        matplot(x = n, y = mat2, type = "l")
        abline(a = 0, b = 0, lty = 3)
        abline(v = 30, lty = 3)
        
        
        d <- par(mfrow = c(2,2))
        
        matplot(x = n, y = mat, type = "l", xlab = "Urchin density", ylab = "dU/dt", main = "Variation in predator density")
        abline(a = 0, b = 0, lty = 3)
        abline(v = 30, lty = 3)
        
        matplot(x = n, y = mat2, type = "l", xlab = "Urchin density", ylab = "dU/dt", main = "Variation in relative body size at max predator density")
        abline(a = 0, b = 0, lty = 3)
        abline(v = 30, lty = 3)

        par(d)


# Now determine how interaction strength (urchin biomass with predators / urchin biomass without) varies with the ratio of predators to prey.

        input <- seq(0, 40, by = 0.01)
        mr <- seq(min(df.ind$mr), max(df.ind$mr), length.out = 100)
        mc <- seq(max(df.ind$mc, na.rm = T), min(df.ind$mc, na.rm = T), length.out = 100)
        vec <- vector()
        dudt <- matrix(nrow = length(input), ncol = length(mr))
        
        # This takes a while to run!
        for(i in 1:length(input)){
          for(j in 1:length(mr)){
          dudt[i,j] <- bd.FR(n = input[i],a = mean(df.ind$median.a, na.rm = T), p = max(lob.a$density, na.rm = T), mc = mc[j], mr = mr[j], t = 1, r = 0.05, k = 30)
          }
        }
        
        temp <- apply(dudt, MARGIN = 2, FUN = function(x){input[x * lead(x) < 0]})
        
        interaction <- data.frame(mr = mr, mc = mc, Uwl = as.vector(temp[1, ]), Unl = 30)
        interaction$int.strength = interaction$Uwl/interaction$Unl
        
        ggplot(interaction, aes(x = mc/mr, y = int.strength))+
          geom_point()


# few things to note... from this analysis we can conclude that dudt has 1 stable biologically relavant equilibrium point. In the absence of predators this equilibria is at K, otherwise predators force the equilibria < k. How strongly predators influence the equilibia is sensitive to r. 
        
# is it worth bringing in size dependency to r and k?
        
        







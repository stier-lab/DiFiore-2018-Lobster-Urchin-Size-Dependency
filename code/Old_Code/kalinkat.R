############################################################
## Size dependence in the functional response
###########################################################

# This code code fits predator and prey size as continuous covariates in the functional response. Derivation of the FR with size dependence based on Kalinkat.

###############################################################

# Get data

df <- read.table(here("data/cleaned","loburc_cleaned.csv"), header = T, sep = ",")

df <- df[,c("id", "size", "treatment", "num_offered", "num_left", "num_dead")]

df[is.na(df)] <- 0

df$killed <- ifelse((df$num_offered - df$num_dead) - df$num_left < 0, 0, (df$num_offered - df$num_dead) - df$num_left)
df <- df[,c("id", "size", "treatment", "num_offered", "killed")]
names(df) <- c("id", "size", "treatment", "initial", "killed")

df$udiam <- ifelse(df$treatment == "urc_medium", 40,
                   ifelse(df$treatment == "urc_large", 60, 20))



#a*(SIZE_C^b)*smearing_estimate*DENSITY
################################################################

# Get weights (g) for predators and prey

df$mr <- 0.000592598*df$udiam^2.872636198*1.01

df$mc <- 0.001352821*df$size^2.913963113*1

#################################################################

#  Hump-shaped allometric functional resopnse fit with ML


hump <- function(initial, mr, mc, h0, a0, mr.exp, mc.exp, beta, eps){
  a <- a0*(mr^beta)*(mc/mr)*(exp(eps*(mc/mr)))
  h <- h0*mr^mr.exp*mc^mc.exp
  a*initial/(1+a*h*initial)
}

hump.prop <- function(initial, mr, mc, h0, a0, mr.exp, mc.exp, beta, eps){
  a <- a0*(mr^beta)*(mc/mr)*(exp(eps*(mc/mr)))
  h <- h0*mr^mr.exp*mc^mc.exp
  -sum(dbinom(killed,prob=pmax(bound,
                               pmin(1-bound, 1/(1/a + h*initial))),size=initial))
  
}

bound <- 1e-9 ## used to bound probabilities between 0 and 1


fit.hump = mle2(hump.prop,
                start=list(h0 = 1, a0 = 1, mr.exp = 1, mc.exp = .75, beta=0.75 , eps= 0.01),
                data=df, 
                optimizer="optim",
                method="L-BFGS-B",
                control=list(maxit=10000)
                )
summary(fit.hump)


fitnls <- nls(hump, )















# Plot it

dat <- data.frame(mr = seq(min(df$mr), max(df$mr), length.out = 1000), mc = seq(max(df$mc), min(df$mc), length.out = 1000))
dat$ratio <- dat$mc/dat$mr

pred.hump.a <- function(mr, mc, a0, beta, eps){
  a <- a0*(mr^beta)*(mc/mr)*exp(eps*(mc/mr))
  a
}

pred.hump.h <- function(mr, mc, h0, mr.exp, mc.exp){
  h <- h0*mr^mr.exp*mc^mc.exp
  h
}


dat$a <- pred.hump.a(mr = dat$mr, mc = dat$mc, 
                     a0 = coef(fit.hump)[2], 
                     beta = coef(fit.hump)[5], 
                     eps = coef(fit.hump)[6])

dat$h <- pred.hump.h(mr = dat$mr, mc = dat$mc, 
                     h0 = coef(fit.hump)[1],
                     mr.exp = coef(fit.hump)[3], 
                     mc.exp = coef(fit.hump)[4])
par(mfrow = c(1,2))
plot(a~ratio, dat, xlab = "Mc/Mr")
plot(h~ratio, dat, xlab = "Mc/Mr")
df$ratio <- df$mc/df$mr


redtoblue <- colorRampPalette(rev(c("#b2182b", "#d6604d", "#f4a582", "#fddbc7", "#f7f7f7","#d1e5f0", "#92c5de","#4393c3","#2166ac")))

redtoblue2 <- colorRampPalette(rev(c('#d73027','#f46d43','#fdae61','#fee090','#ffffbf','#e0f3f8','#abd9e9','#74add1','#4575b4')))



ggplot(dat, aes(x = mr, y = mc))+
  geom_raster(aes(fill = a))+
  scale_fill_gradientn(colors = redtoblue2(20), guide = guide_colorbar(barheight = 15))+
  #scale_x_continuous(breaks = seq(0,20, by = 4)) +
  geom_contour(aes(z = "a"), color = "white", binwidth = 0.1)


+
  geom_point(data = data[data$ratio.stancor > 0, ], aes_string(x = "jitter.dist", y = focal, size = "ratio.stancor"), pch = 21, col = "gray80", bg = adjustcolor("black",alpha.f=0.6))+
  scale_size(range = c(2,6), breaks = seq(0, 1, by = 0.2), guide = guide_legend(title = expression("Observed \nproportion \ngrazed"), order = 1))+
  geom_point(data = data[data$ratio.stancor == 0, ], aes_string(x = "jitter.dist", y = focal), pch = 21, col = "gray80")+
  #scale_y_continuous(trans = scales::log_trans(), breaks = seq(0, 0.03, by = 0.005), labels = seq(0, 0.03, by = 0.005))+
  labs(x = xlab, y = ylab, fill = expression(paste("Predicted \ngrazing"))) +
  theme(legend.spacing = unit(1, "cm"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

}












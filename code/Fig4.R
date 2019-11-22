#-----------------------------------------------------------------------------#
## Figure 4                                                                   #
#-----------------------------------------------------------------------------#
md <- read.csv(here("data", "lob-metadata.csv"))
md$id.n <- as.numeric(md$id)

md <- md %>% arrange(id.n) %>% 
  mutate(name = id, 
         id  = NULL)

df.ind <- df.ind %>% left_join(md, by = c("ind" = "id.n"))

df.ind$ratio <- df.ind$mc / df.ind$mr
df.ind$name.ord <- factor(df.ind$name, levels = df.ind$name[order(df.ind$mc)])

#-----------------------------------------------------------------------------#
# attack rates                                                                #
#-----------------------------------------------------------------------------#
    
    # Tests of MTE expectations and comparison with Rall
    
    mte1.a <- lm(log(median.a) ~ log(mc) + log(mr), data = df.ind)
    summary(mte1.a)
    
    mte2.ratio <- lm(log(median.a) ~ log(I(mc/mr)), data = df.ind)
    summary(mte2.ratio)
    AIC(mte1.a, mte2.ratio)
    
    car::crPlots(mte1.a)

    # Tests of hump shaped relationships based on Kalinkat, McCoy, or Barrios-Oneill 
    
    # the Kalinkat exponential-ricker function 
    
    kal.a <- lm(log(median.a) ~ log(mr) + log(I(mc/mr)) + I(mc/mr), data = df.ind) # this should probably be fit with NLS because of colinearlity between the predictors
    AIC(mte1.a, kal.a)
    
    # is this really necessary? Doesn't seem to be any evidence, wait on this because getting the NLS fits to work is gonna be a pain in the ass particularly as visualization of the data doesn't suppport the more complex models.
    

# Ok so what does this mean. Well the slopes on each predictor (consumer mass and resource mass) represent the allometric scaling exponents. We find no relationship between attack rates and consumer mass and resource mass, or between the attack rates and the body size ratio. This means that there is no evidence for allometric scaling in attack rates with consumer mass, resource mass, or the body size ratio, despite considerable empirical and theoretical evidence based on MTE. Other work has suggested that the relationship would likely be hump shaped. Therefore, we fit two versions of power/exponential-ricker functions to the data based on Kalinkat et al. 2013 and Barrio-Oneill et al. 2016. We see no support for a hump-shaped relationship. 
    
#-----------------------------------------------------------------------------#
# handling times                                                              #
#-----------------------------------------------------------------------------#

    # Tests of MTE expectations and comparison with Rall
    
    mte1.h <- lm(log(median.h) ~ log(mc) + log(mr), data = df.ind)
    summary(mte1.h)
    
    # Does including consumer mass or resource mass add anything to the model? 
    
    mte2.h <- lm(log(median.h) ~ 1, df.ind)
    summary(mte2.h)
    
    mte3.h <- lm(log(median.h) ~ log(mc), df.ind)
    summary(mte2.h)
    
    mte4.h <- lm(log(median.h) ~ log(mr), df.ind)
    summary(mte2.h)
    
    AICtab(mte1.h, mte2.h, mte3.h, mte4.h)
    
    
    mte1.h.test <- lm(log(median.h) ~ scale(log(mc)) + scale(log(mr)), data = df.ind)
    summary(mte1.h.test)
    
    # Does the ratio describe the allometric scaling relationship better than independent scaling on the consumer mass or resource mass?
    mte.h.ratio <- lm(log(median.h) ~ log(ratio), df.ind)
    summary(mte.h.ratio)
    car::crPlots(mte1.h)
    
    AIC(mte1.h, mte.h.ratio) # Doesn't look like it. The MTE model describes the data better.
    
    # Is there any validity to the Barrios-Oneill formulation of the model? 
    
    bo.h <- nls(log(median.h) ~ c1*exp(c2*log(ratio)),
               start = c(c1 = 0.5, c2 = 1), data = df.ind)
    summary(bo.h)
    AIC(mte1.h, bo.h) # still seem like the MTE model desribes the data better. 

#-----------------------------------------------------------------------------#
# Ok, so plot it up                                                           #
#-----------------------------------------------------------------------------#

mod2 <- lm(log10(median.h) ~ log10(ratio), df.ind)
summary(mod2)

newdat <- data.frame(ratio = seq(min(df.ind$ratio, na.rm = T), max(df.ind$ratio, na.rm = T), length.out = 1000))

newdat$predicted <- predict(mod, newdata = newdat)
newdat$predicted2 <- predict(mod2, newdat = newdat, type = "response")

plot(log(median.h) ~ ratio, df.ind)
lines(I(10^predicted2) ~ ratio, newdat)


p1 <- ggplot(df.ind, aes(x = ratio, y = median.a))+
  geom_point(aes(fill = ratio), shape = 21, show.legend = F, size = 2)+
  scale_fill_viridis(alpha = 0.8)+
  scale_y_log10(breaks = c(0.01, 0.02, 0.03), labels = c(0.01, 0.02, 0.03))+
  coord_cartesian(ylim = c(0.01, 0.03))+
  labs(y = "Attack rate", x = "Predator:prey body size ratio")

p2 <- ggplot(df.ind, aes(x = ratio, y = median.h))+
  geom_point(aes(fill = ratio), shape = 21, show.legend = F, size = 2)+
  scale_fill_viridis(alpha = 0.8)+
  scale_y_log10(breaks = c(1,10,100,1000), labels =  c(1,10,100,1000))+
  #geom_line(data = newdat, aes(x = ratio, y = 10^(predicted)))+
  geom_line(data = newdat, aes(x = ratio, y = 10^(predicted2)))+
  labs(y = "Handling time", x = "Predator:prey body size ratio")

p12 <- cowplot::plot_grid(p1,p2, nrow = 1, align = "h")

ggsave(here("figures", "ha-bodysize-fit.png"), p12, width = 8.5, height = 8.5/2)



png(here("figures", "ha-bodysize-fit.png"), width = 1200, height = 600, res = 200)
d <- par(mfrow = c(1,2), mar = c(4,4,1,1))

plot(log10(a+1) ~ ratio, data = df[df$package == "JAGS", ], ylab = "log(a + 1)", xlab = "Predator:prey body size ratio")

plot(log10(h+1) ~ ratio, data = df[df$package == "JAGS", ], ylab = "log(h +1)", xlab = "Predator:prey body size ratio")
lines(predicted ~ ratio, data = newdat)

par(d)
dev.off()



x <- 1:100
y <- 10^(x*-1)
plot(y~x)
plot(log10(y) ~ x)
plot(log10(y)~I(10^x))
plot(y ~ log10(x))
log(-10)
log(0.0001)



mod2 <- lm(log10(y) ~ log10(x))
plot(log10(y) ~ log10(x))
plot(log10(y) ~ I(10))
plot(y ~ I(10^x))

plot(y ~ exp(log10(x)))

summary(mod2)

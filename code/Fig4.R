library(here)
source(here("code", "setup.R"))

#-----------------------------------------------------------------------------#
## Figure 4                                                                   #
#-----------------------------------------------------------------------------#

df.ind <- read.csv(here("data", "JAGSparam-estimates.csv"))
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
    
    kal.a <- nls(median.a ~ b0*mr^b.r*(mc/mr)*exp(eps*(mc/mr)),
                 start = c(b0 = 1.5, b.r = 1, eps = -0.01), 
                 data = df.ind, 
                 na.action = na.exclude)
    summary(kal.a)
    
    mte2.a <- nls(median.a ~ a0*(mr^r.exp)*(mc^c.exp), 
                  start = list(a0 = 0.01, r.exp = -0.002, c.exp = 0.07), 
                  data = df.ind)
    summary(mte2.a)
    
    
    AIC(mte2.a, kal.a) # hmmmmmmmmm....
    
    # Barrios-Oneill et al. 2016... generalized ricker functions from Persson et al. 1998 
    
    bo <- function(mc, mr, beta, gamma, prop){
      log10(a) ~ beta*((log10(mc/mr)/gamma) * exp(1-(log10(mc/mr)/gamma)))^prop
    }
    
    bo.a <- nls(log10(median.a) ~ beta*((log10(mc/mr)/gamma) * exp(1-(log10(mc/mr)/gamma)))^prop, 
                start = c(beta = 1 , gamma = 0.5, prop = 1), data = df.ind)
    summary(bo.a)
    
    
    mte3.a <- lm(log10(median.a) ~ log10(mc) + log10(mr), data = df.ind)
    summary(mte3.a)
    AIC(mte3.a, bo.a)
    
    # is this really necessary? Doesn't seem to be any evidence, wait on this because getting the NLS fits to work is gonna be a pain in the ass particularly as visualization of the data doesn't suppport the more complex models.
    

# Ok so what does this mean. Well the slopes on each predictor (consumer mass and resource mass) represent the allometric scaling exponents. We find no relationship between attack rates and consumer mass and resource mass, or between the attack rates and the body size ratio. This means that there is no evidence for allometric scaling in attack rates with consumer mass, resource mass, or the body size ratio, despite considerable empirical and theoretical evidence based on MTE. Other work has suggested that the relationship would likely be hump shaped. Therefore, we fit a generalized ricker function to the data based on Kalinkat et al. 2013 and Barrio-Oneill et al. 2016. We see no evidence of a hump-shaped relationship. Furthermore, AIC comparison suggests that the simpler model is more parsimonious. 
    
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

    # The goal here is to make a two panel plot of the relationship between attack rate (a) and handling time (b) as a funciton of consumer size at each resource size. Thus a two panel plot with three curves per panel. We are also intersted in plotting what we would have expected to get assuming cross-taxa or within taxa relationships. Thus each plot will include curves that represent the the expectations from Rall et al. 2012. 
    
    # First refit the model of handling time on a log10 sclae
    mte.h <- lm(log10(median.h) ~ log10(mc) + log10(mr), df.ind)
    summary(mte.h)
    
    # Predict the model for plotting purposes
    newdat <- expand.grid(mc = seq(min(df.ind$mc, na.rm = T), max(df.ind$mc, na.rm = T), length.out = 1000), mr = unique(df.ind$mr))
    newdat$median.h <- predict(mte.h, newdata = newdat)
    
    # Ok, so lets see what the relationship between handling time and body size is from the rall paper. These are the exponents for marine invertebrates from Table 1 Rall et al. 2012. Rall reports the intercepts as ln(h0) / ln(a0), and uses units of mg, m2/m3, s. Therefore, I exponentiated and converted to experimental units (hours, m2, g). 
    newdat$rall <- (exp(9.92) * (newdat$mc*1000)^(-0.76) * (newdat$mr*1000)^(0.76)) / (60*60)
    
    
    # add a treatment variable to newdat
    newdat$treatment <- as.factor(ifelse(newdat$mr < 5, "urc_small", ifelse(
      newdat$mr > 5 & newdat$mr <50, "urc_medium", "urc_large"
    )))
    
    
    # We don't see a relationship between attack rate and body size or consumer, resource, or ratio. What would it be though if lobsters following expectations? 
    newdat$rall.a <- (exp(-21.23) * (newdat$mc*1000)^(0.85) * (newdat$mr*1000)^(0.09)) * (60*60)
    
    # Lets add the prediction for attack rate based on the rall equation, despite the lack of significance... 
    mte.a <- lm(log10(median.a) ~ log10(mc) + log10(mr), df.ind)
    summary(mte.a)
    newdat$predicted.a <- predict(mte.a, newdata = newdat)


    p1 <- ggplot(df.ind, aes(x = mc, y = median.a))+
      geom_point(aes(fill = treatment), shape = 21, size = 2, show.legend = F)+
      geom_line(data = newdat, aes(x = mc, y = 10^(predicted.a), color = treatment, linetype = treatment))+
      scale_fill_manual(values = bigsur)+
      geom_line(data= newdat, aes(x = mc, y = rall.a, linetype = treatment), color = "black")+
      scale_color_manual(values = bigsur)+
      scale_y_log10()+
      labs(y = "Attack rate", x = "Consumer mass")+
      facet_wrap(~treatment)+
      theme(
        strip.background = element_blank(),
        strip.text.x = element_blank()
      )
    
    
    ggplot(df.ind, aes(x = mc, y = median.h))+
      geom_point(aes(fill = treatment), shape = 21, size = 2, show.legend = F)+
      facet_wrap(~treatment)
    
    
    p2 <- ggplot(df.ind, aes(x = mc, y = median.h))+
      geom_point(aes(fill = treatment), shape = 21, size = 2)+
      scale_fill_manual(values = bigsur)+
      scale_y_log10(breaks = c(1,10,100,1000), labels =  c(1,10,100,1000))+
      geom_line(data = newdat, aes(x = mc, y = 10^(median.h), color = treatment, linetype = treatment))+
      scale_color_manual(values = bigsur)+
      geom_line(data = newdat, aes(x = mc, y =rall, linetype = treatment), color = "black")+
      labs(y = "Handling time", x = "Consumer mass")+
      facet_wrap(~treatment)+
      theme(
        strip.background = element_blank(),
        strip.text.x = element_blank()
      )
    
    
    
    p12 <- cowplot::plot_grid(p1+theme(legend.position = "none"),p2+theme(legend.position = "none"), nrow = 1, align = "h")
    
    ggsave(here("figures", "ha-bodysize-fit.png"), p12, width = 8.5, height = 8.5/2)


    # plot of attack rate as a function of the body size ratio, as a function of Barrios-Oneil generalized Ricker function and the Kalinkat power exponential...
    newdat2 <- expand.grid(mr = min(df.ind$mr), ratio = seq(min(df.ind$ratio, na.rm = T), max(df.ind$ratio, na.rm = T))) 
    
    newdat2 <- data.frame(mr = seq(max(df.ind$mr), min(df.ind$mr), length.out = 1000), mc = seq(min(df.ind$mc, na.rm = T), max(df.ind$mc, na.rm = T), length.out  = 1000))
    
    newdat2$predict.bo <- predict(bo.a, newdata = newdat2)
    newdat2$predict.kal <- predict(kal.a, newdata = newdat2)
    
    #Generalized Ricker function
    p1.b <- ggplot(df.ind, aes(x = mc/mr, y = median.a))+
      geom_point(aes(fill = treatment), shape = 21, size = 2, show.legend = F)+
      scale_fill_manual(values = bigsur)+
      geom_line(data = newdat2, aes(x = mc/mr, y = 10^predict.bo), color = "gray50", linetype = "dashed")+
      scale_y_log10(breaks = c(0.01, 0.02, 0.03), labels = c(0.01, 0.02, 0.03))+
      coord_cartesian(ylim = c(0.0075, 0.04))+
      scale_x_log10()+
      labs(x = "Predator:prey body size ratio", y = "Attack rate")
    
    # ggplot(df.ind, aes(x = log10(mc/mr), y = log10(median.a)))+
    #   geom_point(aes(fill = treatment), shape = 21, size = 2, show.legend = F)+
    #   geom_line(data = newdat2, aes(x = log10(mc/mr), y = predict.bo))
    
    #Exponential Ricker function
    ggplot(df.ind, aes(x = ratio, y = median.a))+
      geom_point(fill = "gray50", shape = 21, size = 2, show.legend = F)+
      geom_line(data = newdat2, aes(x = mc/mr, y = predict.kal), show.legend = F)
  

p12.b <- cowplot::plot_grid(p1.b, p2+labs(fill = "Prey size", linetype = "Prey size", color = "Prey size"), rel_widths = c(0.33, 0.66))
ggsave(here("figures", "ha-bodysize-fit-v2.png"), p12.b, width = 8.5*1.25, height = 8.5*1.25/2)


#---------------------------------------------------------
## for poster
#---------------------------------------------------------


p1 <- ggplot(df.ind, aes(x = mc, y = median.a))+
  geom_point(aes(fill = treatment), shape = 21, size = 2, show.legend = F)+
  geom_line(data = newdat, aes(x = mc, y = 10^(predicted.a), color = treatment), lwd = 1.2)+
  scale_fill_manual(values = c('#AF8DC3','#C3AF8D','#8DC3AF'))+
  geom_line(data= newdat, aes(x = mc, y = rall.a), color = "white", lwd = 1.2)+
  scale_color_manual(values = c('#AF8DC3','#C3AF8D','#8DC3AF'))+
  scale_y_log10()+
  labs(y = "Attack rate", x = "Consumer mass")+
  facet_wrap(~treatment)+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(), 
    axis.text.x=element_text(angle = 45, hjust = 1))

#AF8DC3

#C3AF8D

#8DC3AF


p2 <- ggplot(df.ind, aes(x = mc, y = median.h))+
  geom_point(aes(fill = treatment), shape = 21, size = 2)+
  scale_fill_manual(values = c('#AF8DC3','#C3AF8D','#8DC3AF'))+
  scale_y_log10(breaks = c(1,10,100,1000), labels =  c(1,10,100,1000))+
  geom_line(data = newdat, aes(x = mc, y = 10^(median.h), color = treatment), lwd = 1.2)+
  scale_color_manual(values = c('#AF8DC3','#C3AF8D','#8DC3AF'))+
  geom_line(data = newdat, aes(x = mc, y =rall), color = "white", lwd = 1.2)+
  labs(y = "Handling time", x = "Consumer mass")+
  facet_wrap(~treatment)+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank(), 
    axis.text.x=element_text(angle = 45, hjust = 1)
  )



p12 <- cowplot::plot_grid(p1+theme(legend.position = "none"),p2+theme(legend.position = "none"), nrow = 1, align = "h")

ggsave(here("figures", "ha-bodysize-fit-poster.svg"), p12, width = 3.3*2, height = 10/3)


#--------------------------------------------------------------------------------
# Jensen's inequality
#--------------------------------------------------------------------------------

# Much of the research to date has evaluated the effect of body size taken at the mean of the predator and prey size distributions. But the relationship between body size and parameters of the FR is nonlinear. Therefore, we would expect that estimates of a and h taken at the mean body size should NOT equal mean estimates of a and h across variation in body size of individuals. This section of code explores these relationships. From the analysis to date, we have shown that there is no relationship between attack rate and body size. BUT there is a strong nonlinear relationship between handling time and body size. Building from these experiemntally derived relationships, we predict the handling time 1) assuming no size-dependence, 2) at the mean body size of predators and prey, and 3) with size-dependence. 

  # size-independent

  predict(mte1.h, newdat = list(mc = 1, mr = 1), se.fit = T) # this isn't correct. To do this accurately, I'd have to fit a null model to the data without any size dependency....

  # evaluate h at the mean of predator and prey mass

  v1 <- exp(predict(mte1.h, newdat = list(mc = mean(df.ind$mc, na.rm = T), mr = mean(df.ind$mr, na.rm = T))))
  
  v1.inv <- 1/v1

  
  # mean of h across variation in mc and mr
  
  v2 <- exp(mean(predict(mte1.h, newdat = list(mc = df.ind$mc, mr = df.ind$mr)), na.rm = T))
  v2.inv <- 1/v2


  v1.inv*24
  v2.inv*24

  
  #NOTE: need to figure out how to boostrap all of this!


#
  df.ind %>% arrange(mr, mc) %>%
    mutate(id = seq(1, length(mr))) %>%
  ggplot(aes(x = mc, y = 1/median.h))+
    geom_linerange(aes(ymin = 1/h.high, ymax = 1/h.low), alpha = 0.75, lwd = 1.5)+
    geom_point(aes(size = mr, color = treatment), show.legend = F)+
    coord_flip()+
    scale_y_log10()+
    facet_wrap(~treatment, scales = "free_y", nrow = 3)+
    geom_hline(yintercept = v1.inv, linetype = "dashed")+
    geom_hline(yintercept = v2.inv, linetype = "solid")+
    labs(x = "Lobster size (g)", y = expression(paste("Maximum consumption rate (ind. h"^"-1"*")")))



exp(9.92) * (mean(newdat$mc, na.rm = T)*1000)^(-0.76) * (mean(newdat$mr*1000)^(0.76)) / (60*60)














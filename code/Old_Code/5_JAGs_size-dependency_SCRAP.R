#------------------------------------------------------------------------
## Outdated code
#------------------------------------------------------------------------




a = MCMCsummary(model,params='a', round = 4)
h = MCMCsummary(model,params='h', round = 4)

df.ind <- data.frame(ind = 1:46, package = "JAGS", median.a = a[,4],
                     a.low = a[,3],
                     a.high = a[,5],
                     median.h = h[,4],
                     h.low = h[,3],
                     h.high = h[,5])
write.csv(df.ind, here("data", "JAGSparam-estimates.csv"), row.names = F)

t.a = MCMCsummary(model,params='t.a', round = 4)
t.h = MCMCsummary(model,params='t.h', round = 4)

mu.a = MCMCsummary(model,params='mu.a', round = 4)
mu.h = MCMCsummary(model,params='mu.h', round = 4)

ids <- distinct(df, id, size, treatment, udiam) %>% arrange(id)

d <- par(mfrow = c(5,5), mar = c(4,4, 1,1))
for(i in 1:46){
  plot(I(killed/48) ~ initial,data = df[as.numeric(df$id) == i, ], xlab="Number of prey",ylab="Consumption rate (prey consumed per predator per hour)", ylim = c(0,1), xlim = c(0, 30), main = ids[i,1])
  curve(holling2(x,a[i,4],h[i,4],P=1,T=1),add=TRUE,col=1,lty=1) #true curve
  text(x = 10, y = 0.6, label = paste("a = ", round(a[i,1], 3), "\nh =", round(h[i,1], 3), "\n lsize = ", ids[i, 2], "\n usize =", ids[i,3], sep = ""))
}
par(d)


treats <- levels(df$treatment)
d <- par(mfrow = c(3,3), mar = c(4,4, 1,1))
for(i in 1:3){
  plot(I(killed/48) ~ jitter(initial),data = df[as.numeric(df$treatment) == i, ], xlab="Number of prey",ylab="Consumption rate (prey consumed per predator per hour)", ylim = c(0,0.6), main = paste(treats[i]))
  curve(holling2(x,t.a[i,4],t.h[i,4],P=1,T=1),add=TRUE,col=1,lty=1) #true curve
}
par(d)


plot(I(killed/48) ~ jitter(initial),data = df, xlab="Number of prey",ylab="Consumption rate (prey consumed per predator per hour)", ylim = c(0,1))
curve(holling2(x,mu.a[1, 4],mu.h[1,4],P=1,T=48),add=TRUE,col=1,lty=1) #true curve

MCMCplot(model, params = "t.a")
MCMCplot(model, params = "t.h", rank = T, xlim = c(-10, 700))


#--------------------------------------------------------------------------
## Coefficient plot arranged by size
#--------------------------------------------------------------------------

# md <- read.csv(here("data", "lob-metadata.csv"))
# md$id.n <- as.numeric(md$id)
# 
# md <- md %>% arrange(id.n) %>% 
#   mutate(name = id, 
#          id  = NULL)
# 
# out <- df %>% left_join(md, by = c("id" = "id.n"))
# 
# out$ratio <- out$mc / out$mr
# out$name.ord <- factor(out$name, levels = out$name[order(out$mc)])
# 
# 
# ggplot(out[out$name != "N07", ])+
#   coord_flip() + 
#   geom_hline(yintercept = 0, colour = gray(1/2), lty = 2)+
#   geom_linerange(aes(x = name.ord, ymin = a.low, ymax = a.high), lwd = 1, position = position_dodge(width = 1/2))+
#   geom_point(aes(x = name.ord, y = mean.a), pch = 20, size = 3, lwd = 1, position = position_dodge(width = 0.5), shape = 21)+
#   theme_bw()+
#   theme(axis.title.y = element_blank(), axis.text = element_text(size = 12), axis.title = element_text(size = 14))
# 
# ggplot(out[out$name != "N07", ])+
#   coord_flip() + 
#   geom_hline(yintercept = 0, colour = gray(1/2), lty = 2)+
#   geom_linerange(aes(x = name.ord, ymin = log(h.low+1), ymax = log(h.high+1)), lwd = 1, position = position_dodge(width = 1/2))+
#   geom_point(aes(x = name.ord, y = log(mean.h+1)), pch = 20, size = 3, lwd = 1, position = position_dodge(width = 0.5), shape = 21)+
#   theme_bw()+
#   theme(axis.title.y = element_blank(), axis.text = element_text(size = 12), axis.title = element_text(size = 14))


# #-----------------------------------------------------------------------------
# ## Figure for lab meeting of FR fits
# #-----------------------------------------------------------------------------
# 
# s <- as.numeric(unique(df$id[df$treatment == "urc_small"]))
# m <- as.numeric(unique(df$id[df$treatment == "urc_medium"]))
# l <- as.numeric(unique(df$id[df$treatment == "urc_large"]))
# 
# 
# png(here("figures", "FR-individual-3panel.png"), width = 1000*2, height = 563*2, res = 300)
# d <- par(mfrow = c(1,3), mar = c(4.5,5,3,1), las = 1
#          #, cex.axis= 1.2, cex.lab = 2, cex.main = 2
#          )
# 
# plot(I(killed/48) ~ jitter(initial), data = df[df$treatment == "urc_small", ], ylim = c(0, 0.7), ylab = "Consumption rate (ind. per predator per hour)", xlab = "", main = "Small urchins (1-3 cm)")
# for(i in s){
#   curve(holling2(x,a[i,4],h[i,4],P=1,T=1),add=TRUE,lty=1, col = "#034e7b")
# }
# curve(holling2(x,t.a[3,4],t.h[3,4],P=1,T=1),add=TRUE,lty=1, col = "red")
# 
# e <- par(mar = c(4.5,2, 3,1))
# plot(I(killed/48) ~ jitter(initial), data = df[df$treatment == "urc_medium", ], ylim = c(0, 0.7), ylab = "", xlab = "Number of prey offered", main = "Medium urchins (3-5 cm)")
# for(i in m){
#   curve(holling2(x,a[i,4],h[i,4],P=1,T=1),add=TRUE,lty=1, col = "#3690c0")
# }
# curve(holling2(x,t.a[2,4],t.h[2,4],P=1,T=1),add=TRUE,lty=1, col = "red")
# 
# plot(I(killed/48) ~ jitter(initial), data = df[df$treatment == "urc_large", ], ylim = c(0, 0.7), ylab = "", xlab = "", main = "Large urchins (5-7 cm)")
# for(i in l){
#   curve(holling2(x,a[i,4],h[i,4],P=1,T=1),add=TRUE,lty=1, col = "#a6bddb")
# }
# curve(holling2(x,t.a[1,4],t.h[1,4],P=1,T=1),add=TRUE,lty=1, col = "red")
# par(e)
# par(d)
# dev.off()

#-----------------------------------------------------------------------------
## MCMCvis plots
#-----------------------------------------------------------------------------

prior <- rnorm(15000*3, 0, 10000)
logit.prior <- exp(prior)/(1+exp(prior))

prior.h <- rnorm(15000*3, 0 , 10000)
exp.prior <- sapply(prior.h, FUN = function(x){exp(max(min(x, 20), -20))})


test <- rnorm(1000, 0, (1/0.1^2))
hist(test)

test2 <- exp(test)/(1+exp(test))
hist(test2)

MCMCtrace(model, params = c('t.logit.a'), 
          ind = TRUE, priors = prior,
          Rhat = TRUE, n.eff = TRUE, 
          post_zm = F)


MCMCtrace(model, params = c('t.log.h'), 
          ind = TRUE, priors = prior.h,
          Rhat = TRUE, n.eff = TRUE, 
          post_zm = F)

MCMCtrace(model, params = c('t.h'), 
          ind = TRUE, priors = prior.h,
          Rhat = TRUE, n.eff = TRUE, 
          post_zm = F)

#-----------------------------------------------------------------------------
## Figure 3
#-----------------------------------------------------------------------------
cols <- rep(colors, 3)

cols <- rep(viridis::viridis(3),3)
linew <- rep(c(3,2,1), 3)
linetp <- rep(c(4,3,1), 3)


png(here("figures", "lob-urc-treatmentlevel-FR.png"), width = 1000*3, height = 333*3, res = 300)
#svg(here("figures", "lob-urc-treatmentlevel-FR.svg"), width = 10, height = 3.33)
d <- par(mfrow = c(1,3), mar = c(4.5,5,3,1), las = 1, bty = "n", bg = NA
         #, cex.axis= 1.2, cex.lab = 2, cex.main = 2
)

plot(killed ~ jitter(initial), data = df[df$treatment == "urc_small", ], ylim = c(0, 27), xlim = c(0,27), ylab = "Number consumed", xlab = "", main = "Small urchins (1.0-3.0 cm)")
for(i in 7:9){
  curve(holling2(x,t.a[i,4],t.h[i,4],P=1,T=48),add=TRUE, col = cols[i], lwd = linew[i], lty = linetp[i], to = 26)
  # curve(holling2(x,t.a[i,3],t.h[i,5],P=1,T=48),add=TRUE, col = cols[i], lty = 2, to = 26)
  # curve(holling2(x,t.a[i,5],t.h[i,3],P=1,T=48),add=TRUE, col = cols[i], lty = 2, to = 26)
}

e <- par(mar = c(4.5,2, 3,1))
plot(killed ~ jitter(initial), data = df[df$treatment == "urc_medium", ], ylim = c(0, 27), xlim = c(0,27), ylab = "", xlab = "Number of prey offered", main = "Medium urchins (3.0-5.0 cm)")
for(i in 4:6){
  curve(holling2(x,t.a[i,4],t.h[i,4],P=1,T=48),add=TRUE, col = cols[i], lwd = linew[i], lty = linetp[i], to = 26)
}

plot(killed ~ jitter(initial), data = df[df$treatment == "urc_large", ], ylim = c(0, 27), xlim = c(0,27), ylab = "", xlab = "", main = "Large urchins (5.0-7.0 cm)")
for(i in 1:3){
  curve(holling2(x,t.a[i,4],t.h[i,4],P=1,T=48),add=TRUE, col = cols[i], lwd = linew[i], lty = linetp[i], to = 26)
}
legend(x = 5, y = 20, legend = c("Large lobsters (>9.0 cm)", "Medium lobsters (7.0-9.0 cm)", "Small lobsters (5.0-7.0 cm)"), col = cols[1:3], lwd = linew[1:3], lty = linetp[1:3], bty = "n", xjust = 0, yjust = 1)

par(e)
par(d)
dev.off()



as.mcmc(model) %>%
  spread_draws(a[t.a])

formerge <- distinct(df, treatment,treatment) %>% arrange(treatment) %>% 
  mutate(i = as.numeric(treatment))

df.plot <- as.mcmc(model) %>%
  spread_draws(t.a[i], t.h[i]) %>%
  group_by(i) %>%
  gather(var, value, -c(.chain, .iteration, .draw, i)) %>%
  left_join(formerge)

a.plot <- df.plot %>% filter(var == "t.a")%>%
  ggplot(aes(y = forcats::fct_rev(lobcat), x = value)) +
  stat_pointintervalh(aes(color = forcats::fct_rev(lobcat)), show.legend = F)+
  scale_color_manual(values = viridis(3))+
  facet_wrap(~treatment, scales = "free")+
  labs(y = "", x = "Attack rate")+
  #theme_pubclean()+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
  )

ggsave(here("figures", "attackrate-coef.png"), device = "png", width = 8.5, height = 8.5/3)
ggsave(here("figures", "attackrate-coef.svg"), device = "svg", width = 8.5, height = 8.5/3)

h.plot <- df.plot %>% filter(var == "t.h")%>%
  ggplot(aes(y = forcats::fct_rev(lobcat), x = value)) +
  stat_pointintervalh(aes(color = forcats::fct_rev(lobcat)), show.legend = F)+
  scale_color_manual(values = viridis(3))+
  facet_wrap(~treatment, scales = "free")+
  scale_x_log10()+
  labs(y = "", x = "Handling time")+
  #theme_pubclean()+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
  )
ggsave(here("figures", "handlingtime-coef.png"), device = "png", width = 8.5, height = 8.5/3)
ggsave(here("figures", "handlingtime-coef.svg"), device = "svg", width = 8.5, height = 8.5/3)

treats <- unique(df.plot$treatment)
plotlist <- list()
for(i in 1:3){
  myfilepath <- paste("figures/coef", i, ".png", sep = "")
  
  p <- df.plot %>% filter(var == "t.a", treatment == !! treats[i]) %>%
    ggplot(aes(y = forcats::fct_rev(lobcat), x = value)) +
    stat_pointintervalh(aes(color = forcats::fct_rev(lobcat)), show.legend = F)+
    scale_color_manual(values = rev(viridis(3)))+
    labs(y = "", x = "")+
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank()
    )
  ggsave(myfilepath, p, width = 8.5/6*2, height = 8.5/3, bg = "transparent")
  plotlist[[i]] = p
}

plotlist <- list()
for(i in 1:3){
  myfilepath <- paste("figures/coef", i+3, ".png", sep = "")
  
  f <- df.plot %>% filter(var == "t.h", treatment == !! treats[i]) %>%
    ggplot(aes(y = forcats::fct_rev(lobcat), x = value)) +
    stat_pointintervalh(aes(color = forcats::fct_rev(lobcat)), show.legend = F)+
    scale_color_manual(values = rev(viridis(3)))+
    scale_x_log10()+
    labs(y = "", x = "")+
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank()
    )
  ggsave(myfilepath, f, width = 8.5/6*2, height = 8.5/3, bg = "transparent")
  plotlist[[i]] = f
}


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






bo.h <- function(ratio, c1,c2){
  c1*10^(c2*log10((ratio)))
}

mod <- nls(log10(median.h) ~ c1*exp(c2*log10(ratio)),
           start = c(c1 = 0.5, c2 = 1), data = df.ind)
summary(mod)

mod2 <- lm(log10(median.h) ~ log10(ratio), df.ind)
summary(mod2)

newdat <- data.frame(ratio = seq(min(df.ind$ratio, na.rm = T), max(df.ind$ratio, na.rm = T), length.out = 1000))

newdat$predicted <- predict(mod, newdata = newdat)
newdat$predicted2 <- predict(mod2, newdat = newdat, type = "response")

plot(log(median.h) ~ ratio, df.ind)
lines(I(10^predicted2) ~ ratio, newdat)

df.ind %>% dplyr::select(name.ord, median.a, median.h, ratio) %>%
  gather(param, value, -c(name.ord, ratio)) %>%
  ggplot(aes(x = ratio, y = value))+
  geom_point()+
  scale_y_log10()+
  facet_wrap(~param, scales = "free")

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

#----------------------------------------------------------------------------------
# for poster
#----------------------------------------------------------------------------------


newdat <- expand.grid(initial = seq(min(df$initial, na.rm = T), max(df$initial, na.rm = T), length.out = 100), newtreat = unique(df$newtreat)) %>% arrange(newtreat) %>%
  separate(newtreat, into = c("treatment", "lobcat"), sep = "[-]")
newdat$a <- rep(t.a$`50%`, each = 100)
newdat$h <- rep(t.h$`50%`, each = 100)
newdat$killed <- holling2(N = newdat$initial, a = newdat$a, h = newdat$h, P = 1, T = 48)

poster1 <- ggplot(df, aes(x = initial, y = killed))+
  geom_jitter(aes(size = treatment, fill = lobcat), pch = 21, show.legend = F, stroke = 1)+
  scale_size_manual(values = c(4,2.5,1.5))+
  scale_fill_manual(values = c('#d53e4f','#fc8d59','#fee08b'))+
  geom_line(data = newdat, aes(x = initial, y = killed, color = lobcat), size = 1.5, show.legend = F)+
  scale_color_manual(values = c('#d53e4f','#fc8d59','#fee08b'))+
  facet_wrap(~treatment)

ggsave(here("figures", "forposter1.svg"), poster1, device = "svg", width = 10, height = 10/3)












mu.logit.a <- rnorm(10000, -5.865, tau.jags(rgamma(10000, shape = 10, rate = 20)))
mu.a <- exp(mu.logit.a)/(1+exp(mu.logit.a))
d <- par(mfrow = c(2,2))
hist(mu.logit.a)
hist(mu.a, breaks = 1000)
hist(logit(df.pop$mu.a))
hist(df.pop$mu.a, breaks = 1000)
par(d)


mu.log.h <- rnorm(10000, 0, tau.jags(0.5))
mu.h <- exp(mu.log.h)
d <- par(mfrow = c(1,2))
hist(mu.log.h)
hist(mu.h, breaks = 1000)
par(d)

mu.log.h <- rnorm(10000, log(dunn), tau.jags(rgamma(10000, shape = 200, rate = 400)))
mu.h <- exp(mu.log.h)
d <- par(mfrow = c(2,2))
hist(mu.log.h)
hist(mu.h, breaks = 1000)
hist(log(df.pop$mu.h))
hist(df.pop$mu.h, breaks = 1000)
par(d)

mu.log.h <- rnorm(10000, 0, tau.jags(runif(10000, 0, 10)))
mu.h <- exp(mu.log.h)
d <- par(mfrow = c(1,2))
hist(mu.log.h)
hist(mu.h, breaks = 1000)
par(d)

hist(runif(10000, 0, 10))
hist(rgamma(10000, shape = 15, rate = 30))
hist(rlnorm(10000, meanlog = log(dunn), sdlog = 1))
mean(rlnorm(10000, meanlog = log(dunn), sdlog = 0.01))

dunn <- 0.741 * 24 #days


df.null <- read.csv(here::here("data/cleaned/posteriors/", "posteriors_null.csv"))

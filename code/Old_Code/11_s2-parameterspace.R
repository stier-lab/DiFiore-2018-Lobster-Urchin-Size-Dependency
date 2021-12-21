#-------------------------------------------
## Setup
#-------------------------------------------

library(here)
source(here("code", "10_purrr_extrapolation.R"))


# Simulation visualization of the interaction between beta.a and beta.h. 

df.exp <- read.table(here("data/cleaned","loburc_cleaned.csv"), header = T, sep = ",") 
df.pop <- read.csv(here::here("data/cleaned/posteriors", "allometric_population.csv")) %>% as_tibble()
df.mte <- read.csv(here::here("data/cleaned/posteriors", "mte_population.csv")) %>% as_tibble()


predict.fun <- function(mc. = mc, mr. = mr, urc.a. = urc.a, lob.a. = lob.a, beta1a., beta2a., beta1h., beta2h., h0. = 11, a0. = -4, T = 1, ...){
  
  loga <- a0. + beta1a.*log(mc.) + beta2a.*log(mr.)
  logh <- h0. + beta1h.*log(mc.) + beta2h.*log(mr.)
  a <- exp(loga)
  h <- exp(logh)
  a*urc.a.*lob.a.*T / (1 + a*h*urc.a.)
  
}

n = 100
lob.a = mean(s$lob.a)
urc.a = mean(s$urc.a)
mc = mean(df.exp$mc, na.rm = T)
mr = mean(df.exp$mr, na.rm = T)
beta1h = seq(-2, 0, length.out = n)
beta1a = seq(-0.1, 0.75, length.out = n)
h0 = mean(c(mean(df.pop$mu.alpha.h), mean(df.mte$mu.alpha.h)))
a0 = mean(c(mean(df.pop$mu.alpha.a), mean(df.mte$mu.alpha.a)))


df <- expand.grid(beta1h = beta1h, beta1a = beta1a)
df$pred_fixR <- predict.fun(beta1a. = df$beta1a, beta1h. = df$beta1h, beta2a. = 0, beta2h. = 0, h0. = h0, a0 = a0)

redtoblue2 <- colorRampPalette(rev(c('#d73027','#f46d43','#fdae61','#fee090','#ffffbf','#e0f3f8','#abd9e9','#74add1','#4575b4')))

expectations <- data.frame(beta1h = -0.75, beta1a = 0.58)
posteriors <- df.pop %>% mean_qi()

p1 <- ggplot(df, aes(x = beta1h, y = beta1a))+
  geom_raster(aes(fill = pred_fixR))+
  scale_fill_gradientn(colors = redtoblue2(20), guide = guide_colorbar(barheight = 15))+
  geom_contour(aes(z = pred_fixR), color = "white", binwidth = 0.001, alpha = 0.25)+
  annotate("point", x = -0.75, y = 0.58, color = "#ebcde0")+
  geom_linerange(data = posteriors, aes(y = beta1.a, xmin = beta1.h.lower, xmax = beta1.h.upper), inherit.aes = F)+
  geom_linerange(data = posteriors, aes(x = beta1.h, ymin = beta1.a.lower, ymax = beta1.a.upper), inherit.aes = F)+
  labs(x = TeX("$\\beta_{c,h}$"), 
       y = TeX("$\\beta_{c,a}$"), 
       fill = "")


beta2h = seq(0, 1.5, length.out = n)
beta2a = seq(-0.5, 0.5, length.out = n)
df <- expand.grid(beta2h = beta2h, beta2a = beta2a)
df$pred_fixC <- predict.fun(beta1a. = 0, beta1h. = 0, beta2a. = df$beta2a, beta2h. = df$beta2h, h0. = h0, a0 = a0)

redtoblue2 <- colorRampPalette(rev(c('#d73027','#f46d43','#fdae61','#fee090','#ffffbf','#e0f3f8','#abd9e9','#74add1','#4575b4')))

expectations <- data.frame(beta2h = -0.75, beta2a = 0.58)
posteriors <- df.pop %>% mean_qi()

p2 <- ggplot(df, aes(x = beta2h, y = beta2a))+
  geom_raster(aes(fill = pred_fixC))+
  scale_fill_gradientn(colors = redtoblue2(20), guide = guide_colorbar(barheight = 15))+
  geom_contour(aes(z = pred_fixC), color = "white", alpha = 0.25)+
  annotate("pointrange", xmin = 0, xmax = 1, x = 0.5, y = 0.33, color = "#ebcde0")+
  geom_linerange(data = posteriors, aes(y = beta2.a, xmin = beta2.h.lower, xmax = beta2.h.upper), inherit.aes = F)+
  geom_linerange(data = posteriors, aes(x = beta2.h, ymin = beta2.a.lower, ymax = beta2.a.upper), inherit.aes = F)+
  labs(x = TeX("$\\beta_{r,h}$"), 
       y = TeX("$\\beta_{r,a}$"), 
       fill = "")



s2 <- cowplot::plot_grid(p1, p2)

ggsave(here::here("figures/", "s_parameterspace.png"), width = 8.5*1.5, height = 8.5/2.5*1.5)




















  

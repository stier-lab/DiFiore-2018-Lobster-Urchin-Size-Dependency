library(here)
library(tidyverse)

#-------------------------------------------------------------------------
## Post-hoc relationship between FR and body size ratios
#-------------------------------------------------------------------------

# read in data 

df <- read.csv(here("data", "paramest-all.csv")) %>% filter(package == "JAGS")

# get lobster metadata

md <- read.csv(here("data", "lob-metadata.csv"))
md$id.n <- as.numeric(as.factor(md$id))

md <- md %>% arrange(id.n) %>% 
  mutate(name = id, 
         id  = NULL)

df <- df %>% left_join(md, by = c("ind" = "id.n"))
df$ratio <- df$mc / df$mr

ggplot(df, aes(x = ratio, y = a))+
  geom_point()+
  # scale_y_continuous(limits = c(0,0.10))+
  facet_wrap(~package, scales = "free")

ggplot(df, aes(x = ratio, y = log(h)))+
  geom_point()+
  # scale_y_continuous(limits = c(0,0.10))+
  facet_wrap(~package, scales = "free")


ggplot(df, aes(x = log(mc), y = log10(a+1)))+
  geom_point(aes(color = treatment))+
  geom_smooth(aes(color = treatment), method = "lm")+
  # scale_y_continuous(limits = c(0,0.10))+
  facet_wrap(~package, scales = "free")

ggplot(df, aes(x = log(mc), y = log10(I(1/h)+1)))+
  geom_point(aes(color = treatment))+
  geom_smooth(aes(color = treatment), method = "lm")+
  # scale_y_continuous(limits = c(0,0.10))+
  facet_wrap(~package, scales = "free")

ggplot(df, aes(x = ratio, y = log10(a+1)))+
  geom_point()


plot(log10(a+1) ~ ratio, data = df[df$package == "JAGS", ], ylab = "log(a + 1)")

lm1 <- lm(log10(a+1) ~ ratio, data = df[df$package == "JAGS", ])
summary(lm1)

#--------------------------------------------------------------------------
## Barrios-Oneill h ~ R 
#--------------------------------------------------------------------------

bo.h <- function(ratio, c1,c2){
  c1*exp(c2*log10((ratio)+1))
}

mod <- nls(log10(h+1) ~ c1*exp(c2*log10((ratio)+1)),
  start = c(c1 = 0.5, c2 = -1), data = df[df$package == "JAGS", ])
summary(mod)

newdat <- data.frame(ratio = seq(min(df$ratio, na.rm = T), max(df$ratio, na.rm = T), length.out = 1000))

newdat$predicted <- predict(mod, newdata = newdat)


png(here("figures", "ha-bodysize-fit.png"), width = 1200, height = 600, res = 200)
d <- par(mfrow = c(1,2), mar = c(4,4,1,1))

plot(log10(a+1) ~ ratio, data = df[df$package == "JAGS", ], ylab = "log(a + 1)", xlab = "Predator:prey body size ratio")

plot(log10(h+1) ~ ratio, data = df[df$package == "JAGS", ], ylab = "log(h +1)", xlab = "Predator:prey body size ratio")
lines(predicted ~ ratio, data = newdat)

par(d)
dev.off()










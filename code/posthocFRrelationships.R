
#-------------------------------------------------------------------------
## Post-hoc relationship between FR and body size ratios
#-------------------------------------------------------------------------

# read in data 

df <- read.csv(here("data", "paramest-all.csv"))

# get lobster metadata

md <- read.csv(here("data", "lob-metadata.csv"))
md$id.n <- as.numeric(md$id)

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

plot(log10(h+1) ~ ratio, data = df)
lines(predicted ~ ratio, data = newdat)





















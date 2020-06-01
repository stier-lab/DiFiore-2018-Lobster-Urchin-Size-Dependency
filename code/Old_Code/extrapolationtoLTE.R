#-------------------------------------------------------------------------------------------
## Predicted size-dependent consumption rates
#-------------------------------------------------------------------------------------------
library(here)
source(here("code", "Functions.R"))

# get the lte urchin size data

urc.s <- read.csv(here("data/LTER", "LTE_Urchin_All_Years.csv"), header = T) %>% # get urchin size data for SBC LTER LTE transects.
  filter(TREATMENT == "CONTROL", COMMON_NAME == "Purple Urchin", YEAR > 2011) %>%
  select(YEAR, MONTH, DATE, SITE, TRANSECT, SP_CODE, SIZE, COUNT) %>% 
  na_if(-99999) %>%
  rename_all(tolower) %>%
  mutate(id = paste(site, transect, sep = "")) %>%
  filter(month == 8) 

#extract filter to prune lobster survey data to only lte control transects

lte_controls <- distinct(urc.s, site, transect) %>% mutate(id = paste(site, transect, sep = ""))

lte_controls <- as.vector(lte_controls$id)


urc.s <- urc.s %>% 
  group_by_at(vars(-count)) %>%
  summarize(count.id = sum(count)) %>%
  complete(count.id = full_seq(1:count.id, 1))%>%
  ungroup() %>%
  mutate(size = as.numeric(size), 
         id = paste(site, transect, sep = ""), 
         count.id = NULL, 
         biomass =  0.000592598*(size*10)^2.872636198*1.01) %>%
  drop_na(size)



# get the community data

urc.a <- read.csv(here("data/LTER", "LTE_All_Species_Biomass_at_transect_20200108.csv"), stringsAsFactors = F,na.strings ="-99999") %>%
  rename_all(tolower) %>%
  dplyr::select(year, month, site, transect, treatment, sp_code, density, dry_gm2) %>%
  mutate(id = paste(site, transect, sep = "")) %>%
  filter(treatment == "CONTROL", sp_code %in% c("SPL"), id %in% lte_controls, year > 2011, month == 8) %>% #only the urchin biomass density data for now..
  mutate(treatment = NULL,
         lid = paste(site, year, sep = "-"), 
         biomass = dry_gm2,
         dry_gm2 = NULL) %>%
  arrange(lid)


# Organize and clean lobster data

lob <- read.csv(here("data/LTER", "Lobster_Abundance_All_Years.csv"), header = T) %>% # get lobster abundance and size data from LTER
  na_if(-99999) %>% 
  mutate(DATE = as.character(DATE),
         SIZE_MM = SIZE_MM) %>%
  rename_all(tolower) %>%
  mutate(sp_code = "PAIN", # add species identifier
         id = paste(site, transect, sep = "")) %>%
  filter(id %in% lte_controls)


lob.a <- lob %>%
  as_tibble() %>%
  group_by(year, month, site, transect, area, sp_code, id) %>%
  mutate(biomass = (0.001352821*(size_mm)^2.913963113)*count) %>%
  summarize(density = sum(count, na.rm = T)/1200, 
            biomass = sum(biomass, na.rm = T)/1200) %>% # estimate biomass density (i.e. per m2 biomass)
  ungroup() %>%
  mutate(area = NULL, 
         lid = paste(site, year, sep = "-"))

lob.s <- lob %>% 
  dplyr::select(-num_ao) %>%
  group_by_at(vars(-replicate,-count)) %>%
  summarize(count.id = sum(count)) %>%
  complete(count.id = full_seq(1:count.id, 1))%>%
  ungroup() %>%
  mutate(size_mm = as.numeric(size_mm), 
         id = paste(site, transect, sep = ""), 
         area = NULL, 
         count.id = NULL, 
         biomass = 0.001352821*(size_mm)^2.913963113) %>%
  drop_na(size_mm)
  


df.a <- bind_rows(urc.a, lob.a)


# Plan: make a unique identifier for each year/site combination. Build a matrix of lobster sizes and urchin sizes where each row is the size distribution for each year/site combination. Run a for loop that applies the boot strap function 1000 times, thereby generating a distribution (n = 1000) of consumption estimates.

l.s <- lob.s %>% 
  mutate(lid = paste(site, year, sep = "-")) %>%
  dplyr::select(lid, biomass) %>%
  group_by(lid) %>% 
  arrange(lid)

fix <- data.frame(lid = "NAPL-2019", biomass = NA)

l.s <- bind_rows(l.s, fix)

l.s <- as.data.frame(l.s)

l.s <- split(l.s$biomass, l.s$lid)

u.s <- urc.s  %>% 
  mutate(lid = paste(site, year, sep = "-")) %>%
  dplyr::select(lid, biomass) %>%
  group_by(lid) %>% 
  arrange(lid)

u.s <- as.data.frame(u.s)
u.s <- split(u.s$biomass, u.s$lid)



bd.FR <- function(n, a, p, t, mc, mr){
  
  h.log <- coef(mte1.h)[1] + coef(mte1.h)[2]*log(mc) + coef(mte1.h)[3]*log(mr)
  h <- exp(h.log)
  
  #predicted consumption rate = 
  a*n*p*t/(1+a*h*n)
}

out <- list()

for(i in 1:length(urc.a$density)){
  
  mc.draw <- sample(l.s[[i]], size = 1000, replace = T)
  mr.draw <- sample(u.s[[i]], size = 1000, replace = T)
  
  n <- as.vector(urc.a$density)[i]
  p <- as.vector(lob.a$density)[i]
  
 out[[i]] <-  bd.FR(n = n, a = 0.02, p = p, t = 1, mc = mc.draw, mr = mr.draw)
  
}

mean.df <- data.frame(lid = urc.a$lid, mean.cr =  unlist(lapply(out, FUN = mean)), median.cr =  unlist(lapply(out, FUN = median)))
hist(mean.df$mean.cr)
hist(mean.df$median.cr)

d <- par(mfrow = c(4,4), mar = c(2,2,2,2))
for(i in 1:39){
  hist(out[[i]], breaks = 30)
}
par(d)


# bring in the LTE community data
cd <- read.csv(here("data/LTER", "LTE_All_Species_Biomass_at_transect_20200108.csv"), stringsAsFactors = F,na.strings ="-99999") %>%
  as_tibble() %>%
  rename_all(tolower) %>%
  dplyr::select(year, month, date, site, transect, treatment, sp_code, dry_gm2, group) %>%
  mutate(id = paste(site, transect, sep = ""), 
         lid = paste(site, year, sep = "-")) %>%
  filter(treatment == "CONTROL", id %in% lte_controls, year > 2011, month == 8)

cd$group2 <- ifelse(cd$sp_code == "MAPY", "KELP", 
                   ifelse(cd$group == "ALGAE" & cd$sp_code != "MAPY", "UA", 
                          ifelse(cd$sp_code == "SPL", "URC", cd$group)))

cd <- cd %>% filter(group2 %in% c("UA", "KELP", "URC")) %>%
  group_by(lid, group2) %>%
  summarize(biomass = sum(dry_gm2)) %>%
  spread(group2, biomass) %>%
  arrange(lid) %>%
  left_join(mean.df) %>%
  separate(lid, into = c("site", "year"), sep = "[-]")

cd <- cbind(cd, urc.d = urc.a$density, lob.d = lob.a$density, lob.b = lob.a$biomass)


psych::pairs.panels(cd[,-c(1,2)], 
                    ellipses = F)


plot(cd[,-c(1,2)])


cd %>% gather(fg, estimate, -c(site, year)) %>%
  filter(fg %in% c("KELP", "UA", "URC", "median.cr", "lob.b")) %>%
  ggplot(aes(x = as.numeric(year), y = estimate))+
  geom_line(aes(color = site))+
  facet_wrap(~fg, ncol = 1, scales = "free")


# Food web biomass based models

  # resource-resource competition
  mod0 <- lmer(UA ~ KELP + (1|site), cd)
  summary(mod0) # not much evidence at least in the LTE control transects summer sampling from 2012-2019
  library(car)
  modelassump(mod0)
  
      mod0.lm <- lm(UA ~ KELP, cd)
      summary(mod0.lm)

  # consumer-resource biomass link 

  mod1 <- lmer(KELP ~ URC + (1|site), cd)
  summary(mod1) # evidence of a negative correlation between urchin biomass and kelp biomass
  modelassump(mod1)
  
  mod2 <- lmer(UA ~ URC + (1|site), cd)
  summary(mod2)
  modelassump(mod2)
  
  # predator-consumer biomass link 
  
  mod3 <- lmer(URC ~ lob.b * as.numeric(year) + (lob.b|site), cd)
  summary(mod3)
  modelassump(mod3)
  
  ggplot(cd, aes(x = lob.b, y = URC))+
    geom_point(aes(color = site))
  
  ggplot(cd)+
    geom_line(aes(x = as.numeric(year), y = scale(lob.b)), col = "red")+
    geom_line(aes(x = as.numeric(year), y = scale(URC)), col = "purple")+
    geom_line(aes(x = as.numeric(year), y = scale(KELP)), col = "forestgreen")+
    facet_wrap(~site, scales = "free")
      
  mod3.0 <- glmer(URC ~ scale(lob.b) + scale(as.numeric(year)) + (1|site), cd, family = gaussian(link = "inverse"))
      summary(mod3.0)
      modelassump(mod3.0)
  
      mod3.lm <- lm(URC ~ lob.b, cd)
      summary(mod3.lm)
      
  # predator- resource biomass indirect link
  
  mod4 <- lmer(KELP ~ lob.b + (1|site), cd)
  summary(mod4)
  modelassump(mod4)
  
  mod5 <- lmer(UA ~ lob.b + (1|site), cd)
  summary(mod5)
  modelassump(mod5)
  
      # mod5.0 <- glmer(UA ~ lob.b + (1|site), cd, family = gaussian(link = "log"))
      # summary(mod5.0)
      # modelassump(mod5)


# ---------------------------------------------
## Stability/Variability analysis
# ---------------------------------------------
  
cd %>% group_by(site) %>%
    summarize(cv.kelp = sd(KELP)/mean(KELP), 
              cv.ua = sd(UA)/mean(UA), 
              mean.mean.cr = mean(mean.cr, na.rm = T)*24) %>%
    ggplot(aes(y = cv.kelp, x = mean.mean.cr))+
    geom_point()
  
cd %>% group_by(site) %>%
  summarize(cv.kelp = sd(KELP)/mean(KELP), 
            cv.ua = sd(UA)/mean(UA), 
            mean.mean.cr = mean(mean.cr, na.rm = T)*24) %>%
  ggplot(aes(y = cv.ua, x = mean.mean.cr))+
  geom_point()

  
# ---------------------------------------------
## Change in kelp and UA analysis 
# ---------------------------------------------


  
  
  
  
  
  
d <- par(mfrow = c(4,4), mar = c(2,2,2,2))
for(i in 1:40){
hist(log(u.s[[i]]))
  
  }
  
par(d)  
  

temp <- lmer(sqrt(KELP) ~ scale(URC) * scale(median.cr) + (1|site), cd)
summary(temp)

temp <- glmer((KELP+0.1) ~ scale(URC) * scale(median.cr) + (1|site), cd, family = Gamma(link = "log"))
summary(temp)


temp <- lmer(UA ~ scale(URC) * scale(median.cr) + (1|site), cd)
summary(temp)

















temp <- lmer(median.cr ~ urc.d + lob.d + (1|site), cd)
summary(temp)

ggplot(cd, aes( x = urc.d, y = median.cr))+
  geom_point(aes(color = site))

ggplot(cd, aes( x = lob.d, y = median.cr))+
  geom_point(aes(color = site))


ggplot(cd, aes(x = log(median.cr), y = KELP))+
  geom_point(aes(color = site))

ggplot(cd, aes(x = median.cr, y = URC))+
  geom_point(aes(color = site))

ggplot(cd, aes(x = URC, y = median.cr))+
  geom_point(aes(color = site))

ggplot(cd, aes(x = log(median.cr), y = UA))+
  geom_point(aes(color = site))


library(lmer)
library(lme4)
library(lmerTest)

lmer1 <- lmer(KELP ~ log(median.cr+1) + (1|site), data = cd)
summary(lmer1)


lmer1 <- lmer(sqrt(UA) ~ log(scale(median.cr)+1) + scale(KELP) + (1|site), data = cd)
summary(lmer1)

lmer1 <- lmer(URC ~ log(median.cr+1) + (1|site), data = cd)
summary(lmer1)















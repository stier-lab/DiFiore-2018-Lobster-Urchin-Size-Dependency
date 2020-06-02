#-------------------------------------------------------------------------------------------
## Predicted size-dependent consumption rates
#-------------------------------------------------------------------------------------------
library(here)
source(here("code", "Base_functions/Functions.R"))
source(here("code", "1_setup.R"))

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
         mass =  0.000592598*(size*10)^2.872636198*1.01) %>%
  drop_na(size) %>% 
  group_by(year, site) %>%
  select(-c(id, size, month, date, transect, sp_code)) %>%
  nest(urc.mass = mass)



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
         area = NULL, 
         count.id = NULL, 
         mass = 0.001352821*(size_mm)^2.913963113) %>%
  drop_na(size_mm) %>% 
  group_by(year, site) %>%
  select(-c(id, size_mm, month, date, transect, sp_code)) %>%
  nest(lob.mass = mass)


# lob.s <- lob.s %>% group_by(year, site, transect, sp_code) %>%
#   select(-c(month, date)) %>%
#   nest()
# 
# urc.s <- urc.s %>% group_by(year, site, transect, sp_code) %>%
#   select(-c(month, date)) %>%
#   nest()



df.a <- bind_rows(urc.a, lob.a) %>% select(year, site, sp_code, biomass) %>%
  spread(sp_code, biomass) %>%
  rename(lob.a = PAIN, urc.a = SPL)

s <- left_join(lob.s, urc.s) %>%
  left_join(df.a) %>%
  mutate(id = paste(year, site, sep = "-"))
#-----------------------------------------------------------------
## 
#-----------------------------------------------------------------

post.a <- read.csv(here::here("data/cleaned", "posteriors_posthoc_a.csv")) %>%
  spread(parameter,estimate) %>% filter(model == "model_wpuncert")
post.h <- read.csv(here::here("data/cleaned", "posteriors_posthoc_h.csv")) %>%
  spread(parameter,estimate) %>% filter(model == "model_wpuncert")

beta1a = sample(post.a$beta1, 10000, replace = T)
beta2a = sample(post.a$beta2, 10000, replace = T)
beta1h = sample(post.h$beta1, 10000, replace = T)
beta2h = sample(post.h$beta2, 10000, replace = T)
h0 = exp(sample(post.h$alpha, 10000, replace = T))
a0 = exp(sample(post.a$alpha, 10000, replace = T))

predict.fun <- function(lob.samples, urc.samples, urc.a, lob.a, beta1a. = beta1a, beta2a. = beta2a, beta1h. = beta1h, beta2h. = beta2h, h0. = h0, a0. = a0, T = 1, ...){
  
  a = a0.*lob.samples^beta1a.*urc.samples^beta2a.
  h = h0.*lob.samples^beta1h.*urc.samples^beta2h.
  a*urc.a*lob.a*T / (1 + a*h*urc.a)
  
}



# k <- predict.fun(mc = sample(lob.size, 10000, replace = T), 
#                  mr = sample(urc.size, 10000, replace = T), 
#                  N = urc.den, 
#                  P = lob.den, 
#                  beta1a = sample(post.a$beta1, 10000, replace = T), 
#                  beta2a = sample(post.a$beta2, 10000, replace = T), 
#                  a0 = exp(sample(post.a$alpha, 10000, replace = T)),
#                  beta1h = sample(post.h$beta1, 10000, replace = T), 
#                  beta2h = sample(post.h$beta2, 10000, replace = T), 
#                  h0 = exp(sample(post.h$alpha, 10000, replace = T)))


names <- s %>% mutate(id = paste(year, site, sep = "-")) %>%
  ungroup() %>%
  select(id) %>%
  distinct(id)

names <- as.vector(names$id)

t <- s %>%
  group_by(year, site) %>%
  mutate(urc.samples = purrr::map(urc.mass, sample_n, 10000, replace = T), 
         lob.samples = purrr::map(lob.mass, sample_n, 10000, replace = T)) %>%
  purrr::pmap(predict.fun) %>% 
  purrr::flatten() %>%
  set_names(names) %>%
  as_tibble() %>%
  gather(id, prediction)
  
t %>%
  ggplot()+
  geom_histogram(aes(x = prediction))+
  facet_wrap(~id, scales = "free")

t %>%
  separate(id, into = c("year", "site"), sep = "-") %>%
  group_by(year, site) %>%
  summarize(mean.prediction = mean(prediction)) %>%
  ggplot(aes(x = as.numeric(year), y = mean.prediction))+
  geom_line(aes(color = site))

t %>%
  separate(id, into = c("year", "site"), sep = "-") %>%
  group_by(year, site) %>%
  median_qi() %>%
  ggplot(aes(x = as.numeric(year), y = prediction))+
  geom_line(aes(color = site))+
  geom_point(aes(color = site))

#-----------------------------------------------------------------------------------------------
## Prediction analysis
#-----------------------------------------------------------------------------------------------

# Ok, we know that interaction strengths can vary, and that body size is one of the reasons that interaction strengths vary. Furthermore, we know that most discussion of interaction strengths is focused at the interspecfic variation in body size, but there is considerable intra-specific variation of body size that can alter how strongly two species interact. I'm interested in how different our estimates of interaction strength are if you ignore intraspecific variation in body size, compared to if you incorporate intraspecific variation. Furthermore, most body size distributions are skewed (i.e. non-normal), so how would our predictions differ for real populations?

# No body size dependence (need to fit this model to compare with)
# Body size dependence, population mean sizes
# Body size dependence, normally-distributed body size distributions
# Body size dependence, skewed-body size distributions
# Body size dependence, observed body size distributions.


# Predicted consumption based only on the mean body size at each site/year

mu <- s %>%
  group_by(year, site) %>%
  mutate(lob.samples = purrr::map_dbl(lob.mass, ~ mean(.x$mass)), 
         urc.samples = purrr::map_dbl(urc.mass, ~ mean(.x$mass))) %>%
  purrr::pmap(predict.fun)%>% 
  set_names(names) %>%
  as_tibble() %>%
  gather(id, predict.mu)
           
mu %>%
  ggplot()+
  geom_histogram(aes(x = predict.mu))+
  facet_wrap(~id, scales = "free")         
         
mu %>% 
  mutate(estimate = "predict.mu") %>%
  rename(value = predict.mu) %>%
  bind_rows(t %>% mutate(estimate = "prediction") %>% rename(value = prediction)) %>%
  group_by(id) %>%
  spread(value, estimate) %>%
  summarize(predict.mu = median(predict.mu), 
            prediction = median(prediction)) %>%
  ggplot(aes(x = predict.mu, y = prediction))+
  geom_point(aes(shape = site, color = as.numeric(year) ))+
  geom_abline(slope =1, intercept = 0, linetype = "dashed", color = "gray")
         
mu %>% arrange(id) %>% bind_cols(select(t, prediction)) %>%
  group_by(id) %>%
  gather(estimate, value, -id) %>%
  filter(id %in% c("2012-IVEE", "2013-IVEE")) %>%
  ggplot()+
  geom_density(aes(x = value, fill = estimate), adjust = 10)+
  facet_wrap(~id)

mu %>% arrange(id) %>% bind_cols(select(t, prediction)) %>%
  group_by(id) %>%
  gather(estimate, value, -id) %>%
  separate(id, into = c("year", "site"), sep = "-") %>%
  ggplot()+
  ggridges::geom_density_ridges(aes(x = value, y = as.factor(year), fill = as.factor(estimate)), rel_min_height = 0.001)+
  facet_wrap(~site, scales = "free")+
  coord_cartesian(xlim = c(0, 1))
  
mu %>% 
  mutate(estimate = "predict.mu") %>%
  rename(value = predict.mu) %>%
  bind_rows(t %>% mutate(estimate = "prediction") %>% rename(value = prediction)) %>%
  separate(id, into = c("year", "site"), sep = "-") %>%
  group_by(site, year, estimate) %>%
  median_qi() %>%
  ggplot(aes(x = value, y = year))+
  geom_pointintervalh(aes(color = estimate), position = position_dodge(width = 0.5))+
  facet_wrap(~site, scales = "free")
  
  geom_line(aes(color = estimate))+
  geom_point(aes(color = site))



df.mod <- mu %>% arrange(id) %>% bind_cols(select(t, prediction)) %>%
  group_by(id) %>%
  summarize(predict.mu = median(predict.mu), 
            prediction = median(prediction))
plot(prediction ~ predict.mu, df.mod)
lm1 <- lm(prediction ~ predict.mu, df.mod)
summary(lm1)
plot(resid(lm1) ~ predict.mu, df.mod)













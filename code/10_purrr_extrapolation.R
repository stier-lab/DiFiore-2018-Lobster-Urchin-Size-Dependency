#-------------------------------------------------------------------------------------------
## Predicted size-dependent consumption rates
#-------------------------------------------------------------------------------------------
library(here)
source(here("code", "Base_functions/Functions.R"))
source(here("code", "1_setup.R"))

# Tank size

# Full tank is 137 x 76. I'm going to stick with areal measures (rather than volumetric) because LTER works off of an areal basis. Therefore experimental tanks were 68.5 x 76. 

tsize <- 137/2 * 76 /10000

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



df.a <- bind_rows(urc.a, lob.a) %>% select(year, site, sp_code, density) %>%
  spread(sp_code, density) %>%
  rename(lob.a = PAIN, urc.a = SPL)

s <- left_join(lob.s, urc.s) %>%
  left_join(df.a) %>%
  mutate(id = paste(year, site, sep = "-"))

#-----------------------------------------------------------------
## Simulation function
#-----------------------------------------------------------------

post <- read.csv(here::here("data/cleaned/posteriors", "allometric_population.csv"
)) %>% as_tibble() %>% sample_draws(10000)

# beta1a = sample(post$beta1, 10000, replace = T)
# beta2a = sample(post.a$beta2, 10000, replace = T)
# beta1h = sample(post.h$beta1, 10000, replace = T)
# beta2h = sample(post.h$beta2, 10000, replace = T)
# h0 = exp(sample(post.h$alpha, 10000, replace = T))
# a0 = exp(sample(post.a$alpha, 10000, replace = T))

predict.fun <- function(lob.samples, urc.samples, urc_density, lob_density, beta1a. = post$beta1.a, beta2a. = post$beta2.a, beta1h. = post$beta1.h, beta2h. = post$beta2.h, h0. = post$mu.alpha.h, a0. = post$mu.alpha.a, T = 1, ...){
  
  loga <- a0. + beta1a.*log(lob.samples) + beta2a.*log(urc.samples)
  logh <- h0. + beta1h.*log(lob.samples) + beta2h.*log(urc.samples)
  a <- exp(loga)/tsize
  h <- exp(logh)
  a*urc_density*lob_density*T / (1 + a*h*urc_density)
  
}

names <- s %>% mutate(id = paste(year, site, sep = "-")) %>%
  ungroup() %>%
  select(id) %>%
  distinct(id)

names <- as.vector(names$id)

# t <- s %>%
#   group_by(year, site) %>%
#   mutate(urc.samples = purrr::map(urc.mass, sample_n, 10000, replace = T), 
#          lob.samples = purrr::map(lob.mass, sample_n, 10000, replace = T)) %>%
#   purrr::pmap(predict.fun) %>% 
#   purrr::flatten() %>%
#   set_names(names) %>%
#   as_tibble() %>%
#   gather(id, prediction)
#   
# # t %>%
# #   ggplot()+
# #   geom_histogram(aes(x = prediction))+
# #   facet_wrap(~id, scales = "free")
# 
# t %>%
#   separate(id, into = c("year", "site"), sep = "-") %>%
#   group_by(year, site) %>%
#   summarize(mean.prediction = mean(prediction)) %>%
#   ggplot(aes(x = as.numeric(year), y = mean.prediction))+
#   geom_line(aes(color = site))
# 
# t %>%
#   separate(id, into = c("year", "site"), sep = "-") %>%
#   group_by(year, site) %>%
#   median_qi() %>%
#   ggplot(aes(x = as.numeric(year), y = prediction))+
#   geom_line(aes(color = site))+
#   geom_point(aes(color = site))

#-----------------------------------------------------------------------------------------------
## Prediction analysis
#-----------------------------------------------------------------------------------------------

# Ok, we know that interaction strengths can vary, and that body size is one of the reasons that interaction strengths vary. Furthermore, we know that most discussion of interaction strengths is focused at the interspecfic variation in body size, but there is considerable intra-specific variation of body size that can alter how strongly two species interact. I'm interested in how different our estimates of interaction strength are if you ignore intraspecific variation in body size, compared to if you incorporate intraspecific variation. Furthermore, most body size distributions are skewed (i.e. non-normal), so how would our predictions differ for real populations?

# No body size dependence
# Body size dependence, population mean sizes
# Body size dependence, observed body size distributions.
# Body size dependence - according to MTE

#-------------------------------
# No body size dependence
#-------------------------------

df.null <- read.csv(here::here("data/cleaned/posteriors", "posteriors_null.csv")) %>% as_tibble() %>% sample_draws(10000)

predict.FR <- function(urc_density, lob_density, a = df.null$a, h = df.null$h, ...){
  (a/tsize)*urc_density*lob_density/(1 + (a/tsize)*h*urc_density)
}

null <- s %>%
  group_by(year, site) %>%
  purrr::pmap(predict.FR)%>% 
  set_names(names) %>%
  as_tibble() %>%
  gather(id, prediction) %>%
  separate(id, into = c("year", "site"), sep = "[-]") %>%
  mutate(estimate = "null")

ggplot(null)+
  ggridges::geom_density_ridges(aes(x = prediction, y = as.factor(year)), rel_min_height = 0.001)+
  facet_wrap(~site, scales = "free")
  coord_cartesian(xlim = c(0, 1))

#--------------------------------
## Mean allometric scaling
#--------------------------------

mu <- s %>%
  group_by(year, site) %>%
  mutate(lob.samples = purrr::map_dbl(lob.mass, ~ mean(.x$mass)), 
         urc.samples = purrr::map_dbl(urc.mass, ~ mean(.x$mass))) %>%
  purrr::pmap(predict.fun)%>% 
  set_names(names) %>%
  as_tibble() %>%
  gather(id, prediction) %>%
  separate(id, into = c("year", "site"), sep = "[-]") %>%
    mutate(estimate = "mu")
  
           
ggplot(mu)+
  ggridges::geom_density_ridges(aes(x = prediction, y = as.factor(year)), rel_min_height = 0.001)+
  facet_wrap(~site, scales = "free")
         


#------------------------------------------------------------------------
## Allometric scaling conditional on observed body size distributions
#------------------------------------------------------------------------


full <- s %>%
  group_by(year, site) %>%
  mutate(urc.samples = purrr::map(urc.mass, sample_n, 10000, replace = T), 
         lob.samples = purrr::map(lob.mass, sample_n, 10000, replace = T)) %>%
  purrr::pmap(predict.fun) %>% 
  purrr::flatten() %>%
  set_names(names) %>%
  as_tibble() %>%
  gather(id, prediction) %>%
  separate(id, into = c("year", "site"), sep = "[-]") %>%
  mutate(estimate = "full")

ggplot(full)+
  ggridges::geom_density_ridges(aes(x = prediction, y = as.factor(year)), rel_min_height = 0.001)+
  facet_wrap(~site, scales = "free")

#------------------------------------------------------------------------
# Empirical expectations of allometric scaling from Rall et al.
#------------------------------------------------------------------------

#df.mte <- read.csv(here::here("data/cleaned/posteriors", "mte_population.csv")) %>% as_tibble() %>% sample_draws(10000)

# Predict solely based on the exponents for marine invertebrates from Rall et al. 2012. Need to deal with units. Rall uses mg, m2, s. The relationship is: 


predict.RALL <- function(lob.samples, urc.samples, urc_density, lob_density, beta1a. = 0.85, beta2a. = 0.09, beta1h. = -0.76, beta2h. = 0.76, h0. = 10.38, a0. = -21.23, T = 1, ...){
  
  loga <- a0. + beta1a.*log(lob.samples) + beta2a.*log(urc.samples)
  logh <- h0. + beta1h.*log(lob.samples) + beta2h.*log(urc.samples)
  a <- exp(loga)*(60*60)
  h <- exp(logh)/(60*60)
  a*urc_density*lob_density*T / (1 + a*h*urc_density)
  
}

RALL <- s %>%
  group_by(year, site) %>%
  mutate(lob.samples = purrr::map_dbl(lob.mass, ~ mean(.x$mass)), 
         urc.samples = purrr::map_dbl(urc.mass, ~ mean(.x$mass))) %>%
  purrr::pmap(predict.RALL) %>% 
  #purrr::flatten() %>%
  set_names(names) %>%
  as_tibble() %>%
  gather(id, prediction) %>%
  separate(id, into = c("year", "site"), sep = "[-]") %>%
  mutate(estimate = "Rall")

#------------------------------------------------------------------------
# Theoretical expectations
#------------------------------------------------------------------------

# Theory doesn't provide a method of estimating the intercepts, therefore, we used our model to generate an estimate of h0 and a0 given the theoretical expectations of the scaling exponents by forcing the variance in the priors of the parameter estimates to ~0. 
df.mte <- read.csv(here::here("data/cleaned/posteriors", "mte_population.csv")) %>% as_tibble() %>% sample_draws(10000)

# Alternatively, we could borrow from the Rall paper to get the intercepts and fix the exponents.

predict.MTE <- function(lob.samples, urc.samples, urc_density, lob_density, beta1a. = 0.58, beta2a. = 0.33, beta1h. = -0.75, beta2h. = 0.5, h0. = df.mte$mu.alpha.h, a0. = df.mte$mu.alpha.a, T = 1, ...){
  
  loga <- a0. + beta1a.*log(lob.samples) + beta2a.*log(urc.samples)
  logh <- h0. + beta1h.*log(lob.samples) + beta2h.*log(urc.samples)
  a <- exp(loga)
  h <- exp(logh)
  a*urc_density*lob_density*T / (1 + a*h*urc_density)
  
}

MTE <- s %>%
  group_by(year, site) %>%
  mutate(lob.samples = purrr::map_dbl(lob.mass, ~ mean(.x$mass)), 
         urc.samples = purrr::map_dbl(urc.mass, ~ mean(.x$mass))) %>%
  purrr::pmap(predict.MTE) %>% 
  #purrr::flatten() %>%
  set_names(names) %>%
  as_tibble() %>%
  gather(id, prediction) %>%
  separate(id, into = c("year", "site"), sep = "[-]") %>%
  mutate(estimate = "MTE")


MTE.fixed <- s %>%
  group_by(year, site) %>%
  mutate(lob.samples = purrr::map_dbl(lob.mass, ~ mean(.x$mass)), 
         urc.samples = purrr::map_dbl(urc.mass, ~ mean(.x$mass))) %>%
  purrr::pmap(predict.MTE, h0. = 10.38, a0. = -21.23) %>% 
  #purrr::flatten() %>%
  set_names(names) %>%
  as_tibble() %>%
  gather(id, prediction) %>%
  separate(id, into = c("year", "site"), sep = "[-]") %>%
  mutate(estimate = "MTE.fixed")

#------------------------------------------------------------------------
# Allometric scaling of prey size only!
#------------------------------------------------------------------------

# df.preyonly <- read.csv(here::here("data/cleaned/posteriors", "preyonly_population.csv")) %>% as_tibble() %>% sample_draws(10000)
# 
# predict.prey <- function(lob.samples, urc.samples, urc.a, lob.a, beta2a. = df.preyonly$beta2.a, beta2h. = df.preyonly$beta2.h, h0. = df.preyonly$mu.alpha.h, a0. = df.preyonly$mu.alpha.a, T = 1, ...){
#   
#   loga <- a0. + beta2a.*log(urc.samples)
#   logh <- h0. + beta2h.*log(urc.samples)
#   a <- exp(loga)
#   h <- exp(logh)
#   a*urc.a*lob.a*T / (1 + a*h*urc.a)
#   
# }
# 
# prey <- s %>%
#   group_by(year, site) %>%
#   mutate(lob.samples = purrr::map_dbl(lob.mass, ~ mean(.x$mass)), 
#          urc.samples = purrr::map_dbl(urc.mass, ~ mean(.x$mass))) %>%
#   purrr::pmap(predict.prey) %>% 
#   #purrr::flatten() %>%
#   set_names(names) %>%
#   as_tibble() %>%
#   gather(id, prediction) %>%
#   separate(id, into = c("year", "site"), sep = "[-]") %>%
#   mutate(estimate = "prey")

#--------------------------------------------------------
# Combine for plotting
#--------------------------------------------------------

df <- rbind(null, mu, full, RALL, MTE.fixed)
write.csv(df, here::here("data/cleaned/posteriors", "observational_predictions.csv"), row.names = F)
df <- read.csv(here::here("data/cleaned/posteriors/", "observational_predictions.csv"))

ggplot(df)+
  ggridges::geom_density_ridges(aes(x = prediction, y = as.factor(year), fill = estimate), rel_min_height = 0.001)+
  facet_wrap(~site, scales = "free")

# fully facetted plot for supplement (?)
p6 <- df %>%
  group_by(year, site, estimate) %>%
  median_qi(.width = c(.95,.75)) %>%
  ungroup() %>%
  mutate(estimate = recode(estimate, full = "5. Fully integrated allometric model", mu = "4. Allometric model w/ mean body size", null = "1. Unstructured model", Rall = "3. Empirical expectations", MTE = "2. Theoretical expectations")) %>% 
  ggplot(aes(x = prediction, y = forcats::fct_rev(year)))+
  geom_pointintervalh(aes(color = estimate), position = position_dodge(width = 0.5))+
  scale_color_manual(values = c('#AF8DC3','#C3AF8D', '#c3958d', '#8DC3AF', 
                                #'#c38db5'
                                "black"))+
  facet_wrap(~site, scales = "free")+
  labs(x = expression(paste("Predicted consumption rate (ind. m"^-2,"h"^-1,")")), y = "", color = "")+
  theme(legend.position = c(0.7, 0.3))

ggsave(here::here("figures/", "observational.png"), p6, width = 3*4, height = 3*2)



# plot for main text that shows variation but points to differences in mean estimates from models
groups <- df %>% 
  mutate(id = paste(year, site, sep = "")) %>%
  group_by(year, site, estimate, id) %>%
  mean_qi(prediction)

# mean estimate averaged across sites and years. Variance represents spatio-temporal variance NOT parameter or size distribution variance. 
mean <- groups %>% 
  group_by(estimate) %>%
  mean_qi(prediction)
  
coef <- ggplot(groups, aes(x = prediction, y = estimate))+
  geom_pointintervalh(aes(group = id, color = estimate), position = position_dodge(0.25), alpha = 0.25)+
  scale_color_manual(values = c('#AF8DC3','#C3AF8D', '#c3958d', '#8DC3AF', '#c38db5'))+
  geom_pointintervalh(data = mean, aes(x = prediction, y = estimate)) +
  labs(y = "", x = expression(paste("Predicted consumption rate (ind. m"^-2,"h"^-1,")")))


# time series plot
timeseries <- groups %>% filter(estimate == "full") %>%
ggplot(aes(x = as.numeric(year), y = prediction))+
  geom_line(aes(color = site))+
  geom_point(aes(color = site))+
  scale_color_manual(values = c('#AF8DC3','#C3AF8D', '#c3958d', '#8DC3AF', '#c38db5'))+
  labs(x = "", y = expression(paste("Predicted consumption rate (ind. m"^-2,"h"^-1,")")), color = "Site")+
  theme(legend.position = c(0.1, 0.8))
  

p6b <- plot_grid(timeseries, coef, nrow = 1)
ggsave(here::here("figures/", "observational_alt.png"), p6b, width = 3*4, height = 3*2)


#----------------------------------------------------------------------------
## Scrap
#----------------------------------------------------------------------------

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



df.mod <- mu %>% 
  mutate(estimate = "predict.mu") %>%
  rename(value = predict.mu) %>%
  bind_rows(t %>% mutate(estimate = "prediction") %>% rename(value = prediction)) %>%
  separate(id, into = c("year", "site"), sep = "-") %>%
  group_by(year, site, estimate) %>%
  spread(estimate, value)

ggplot(df.mod, aes(x = ))

plot(prediction ~ predict.mu, df.mod)
lm1 <- lm(prediction ~ predict.mu, df.mod)
summary(lm1)
plot(resid(lm1) ~ predict.mu, df.mod)


temp <- groups %>% select(year, site, prediction, estimate) %>% pivot_wider(values_from = prediction, names_from = estimate)

plot( mu ~ null, temp)
abline(0,1)

plot(mu ~ Rall, temp)
abline(0, 1)

plot(mu ~ MTE, temp)
abline(0, 1)

mean(temp$null)/sd(temp$null)
mean(temp$mu)/sd(temp$mu)
var(df$prediction[df$estimate == "mu"])
var(df$prediction[df$estimate == "null"])



ggplot(df, aes(x = ))









temp <- do.call("rbind", replicate(3, mu, simplify = FALSE)) %>% select(-estimate) %>% rename(mu = "prediction")
temp2 <- rbind(null[, c("prediction", "estimate")], full[, c("prediction", "estimate")], MTE[, c("prediction", "estimate")])

temp$alt <- temp2$prediction
temp$estimate <- temp2$estimate


ggplot(temp, aes(x = mu, y = alt))+
  geom_point(aes(color = estimate), shape = 21, alpha = 0.5)+
  geom_abline(slope = 1, intercept = 0)

#------------------------------------------------------------------------
## 1:1 plot comparing MU estiamtes to other estimates
#------------------------------------------------------------------------

sum <- df %>% group_by(year, site, estimate) %>%
  summarize(mean = mean(prediction)) %>%
  pivot_wider(names_from = estimate, values_from = mean) %>%
  pivot_longer(cols = c(null, full, MTE, Rall, MTE.fixed), names_to = "estimate", values_to = "prediction") %>%
  filter(estimate != "MTE")


fig6.p3 <- df %>% group_by(year, site, estimate) %>% mutate(.draw = 1:n()) %>% sample_draws(500) %>%
  pivot_wider(names_from = estimate, values_from = prediction) %>%
  pivot_longer(cols = c(null, full, MTE), names_to = "estimate", values_to = "prediction") %>%
  filter(estimate != "MTE") %>%
  ggplot(aes(x = mu, y = prediction))+
  geom_point(aes(fill = estimate), shape = 21, alpha = 0.25, color = "white")+
  geom_point(data = sum, aes(x  = mu, y = prediction, fill = estimate), size = 2, shape = 21, color = "black")+
  geom_abline(slope = 1, intercept = 0, linetype = "dashed")+
  ylim(0, 0.0035)+ # this cuts off a whole bunch of points!
  xlim(0, 0.0025)+
  labs(y = expression(paste("Predicted consumption rate (ind. m"^-2,"h"^-1,")")), x = expression(paste("Predicted consumption rate for mean experimental model (ind. m"^-2,"h"^-1,")")), fill = "Model")+
  theme(legend.position = c(0.75, 0.2))
  # scale_x_log10()+
  # scale_y_log10()+

ggsave(here::here("figures", "fig6panel3.png"), fig6.p3, device = "png")

p6_bc <- plot_grid(timeseries, fig6.p3, align = "h")

ggsave(here::here("figures/", "fig6panelbc.png"), p6_bc, width = 8.25*2, height = 8.25/3*2)






df %>% group_by(estimate) %>%
  summarize(mean = mean(prediction), 
            ci.low = PI(prediction)[1], 
            ci.high = PI(prediction)[2]) %>%
  ggplot(aes(x = estimate, y = mean))+
    geom_point()+
    geom_errorbar(aes(ymin = ci.low, ymax = ci.high))+
    labs(x = "Model type", y = "Predicted consumption")



#--------------------------------------------------------------------------------------
## Extremes plots
#-------------------------------------------------------------------------------------
l.sum <- lob.s %>% unnest(cols = c(lob.mass)) %>% group_by(year, site) %>%
  summarize(mean = mean(mass))

summary(l.sum$mean)
hist(l.sum$mean)


lob.s %>% unnest(cols = c(lob.mass)) %>% group_by(year, site) %>%
  summarize(mean = mean(mass)) %>%
  ungroup() %>%
  drop_na(mean) %>%
  filter(mean == max(mean) | 
           mean == min(mean) |
           abs(mean - mean(mean)) == min(abs(mean - mean(mean))))

filter(abs(b - 1.43) == min(abs(b - 1.43)))




u.sum <- urc.s %>% unnest(cols = c(urc.mass)) %>% group_by(year, site) %>%
  summarize(mean = mean(mass))
summary(u.sum$mean)
hist(u.sum$mean)

urc.s %>% unnest(cols = c(urc.mass)) %>% group_by(year, site) %>%
  summarize(mean = mean(mass)) %>%
  ungroup() %>%
  filter(mean == max(mean) | 
         mean == min(mean) |
           abs(mean - mean(mean)) == min(abs(mean - mean(mean))))


ex <- s %>% 
  mutate(avg_lsize = map_dbl(lob.mass, ~{mean(.x$mass)}), 
         avg_usize = map_dbl(urc.mass, ~{mean(.x$mass)}), 
         diff = avg_lsize - avg_usize) %>%
  ungroup() %>%
  filter(diff == max(diff) | 
           diff == min(diff) |
           abs(diff -mean(diff)) == min(abs(diff - mean(diff)))) %>%
  select(year, site, lob.mass, urc.mass, urc.a, id)

names <- as.vector(ex$id)
urc.a.mean <- mean(s$urc.a)
  

null <- ex %>%
  group_by(year, site) %>%
  purrr::pmap(predict.FR, lob.a = 1)%>% 
  set_names(names) %>%
  as_tibble() %>%
  gather(id, prediction) %>%
  separate(id, into = c("year", "site"), sep = "[-]") %>%
  mutate(estimate = "null")


mu <- ex %>%
  group_by(year, site) %>%
  mutate(lob.samples = purrr::map_dbl(lob.mass, ~ mean(.x$mass)), 
         urc.samples = purrr::map_dbl(urc.mass, ~ mean(.x$mass))) %>%
  purrr::pmap(predict.fun, lob.a = 1)%>% 
  set_names(names) %>%
  as_tibble() %>%
  gather(id, prediction) %>%
  separate(id, into = c("year", "site"), sep = "[-]") %>%
  mutate(estimate = "mu")


full <- ex %>%
  group_by(year, site) %>%
  mutate(urc.samples = purrr::map(urc.mass, sample_n, 10000, replace = T), 
         lob.samples = purrr::map(lob.mass, sample_n, 10000, replace = T)) %>%
  purrr::pmap(predict.fun, lob.a = 1) %>% 
  purrr::flatten() %>%
  set_names(names) %>%
  as_tibble() %>%
  gather(id, prediction) %>%
  separate(id, into = c("year", "site"), sep = "[-]") %>%
  mutate(estimate = "full")

RALL <- ex %>%
  group_by(year, site) %>%
  mutate(lob.samples = purrr::map_dbl(lob.mass, ~ mean(.x$mass)), 
         urc.samples = purrr::map_dbl(urc.mass, ~ mean(.x$mass))) %>%
  purrr::pmap(predict.RALL, lob.a = 1) %>% 
  #purrr::flatten() %>%
  set_names(names) %>%
  as_tibble() %>%
  gather(id, prediction) %>%
  separate(id, into = c("year", "site"), sep = "[-]") %>%
  mutate(estimate = "Rall")


MTE <- ex %>%
  group_by(year, site) %>%
  mutate(lob.samples = purrr::map_dbl(lob.mass, ~ mean(.x$mass)), 
         urc.samples = purrr::map_dbl(urc.mass, ~ mean(.x$mass))) %>%
  purrr::pmap(predict.MTE, lob.a = 1) %>% 
  #purrr::flatten() %>%
  set_names(names) %>%
  as_tibble() %>%
  gather(id, prediction) %>%
  separate(id, into = c("year", "site"), sep = "[-]") %>%
  mutate(estimate = "MTE")


MTE.fixed <- ex %>%
  group_by(year, site) %>%
  mutate(lob.samples = purrr::map_dbl(lob.mass, ~ mean(.x$mass)), 
         urc.samples = purrr::map_dbl(urc.mass, ~ mean(.x$mass))) %>%
  purrr::pmap(predict.MTE, h0. = 10.38, a0. = -21.23, lob.a = 1) %>% 
  #purrr::flatten() %>%
  set_names(names) %>%
  as_tibble() %>%
  gather(id, prediction) %>%
  separate(id, into = c("year", "site"), sep = "[-]") %>%
  mutate(estimate = "MTE.fixed")


df.x <- rbind(null, full, mu, MTE, MTE.fixed, RALL) %>%
  mutate(id = paste(year, site, sep = "-")) %>%
  filter(estimate != "MTE", estimate != "full") %>%
  mutate(names = as.factor(recode(id, "2013-CARP" = "Large predators, small prey", "2014-AQUE" = "Mean predators, mean prey", "2015-MOHK" = "Small predators, large prey" ))) %>%
  mutate(estimate = forcats::fct_relevel(as.factor(estimate), "full", "mu", "null", "MTE.fixed", "Rall"), 
         estimate = recode(estimate, full = "Fully-integrated\nsize-structured", null = "Unstructured", MTE.fixed = "First principles\nexpectations", Rall = "Metanalysis\nexpectations", mu = "Mean\nsize-structured"))

df.x.sum <- df.x %>% group_by(estimate, names) %>%
  filter(estimate != "Metanalysis\nexpectations", estimate != "First principles\nexpectations") %>%
  mean_qi(prediction, .width = c(0.75, 0.95))

extremes <- ggplot(df.x, aes(x = estimate, y = prediction))+
  geom_jitter(aes(color = estimate), show.legend = F)+
  geom_pointinterval(data = df.x.sum)+
  scale_color_manual(values = c('#AF8DC3','#C3AF8D', '#c3958d', '#8DC3AF', '#c38db5'))+
  facet_wrap(~names)+
  theme(axis.text.x = element_text(angle = 33, vjust = 0.9, hjust=1))+
  labs(x = "", y = expression(paste("Predicted consumption (ind. m"^-2,"h"^-1,")", sep = "")))


ggsave(here::here("figures/", "extremes.png"), extremes, device = "png", width = 8.5*1.5, height = 8.5*1.5/3)



#-----------------------------------------------------------------------------
## Residual based analysis
#-----------------------------------------------------------------------------

temp <- df %>%
  filter(estimate == "null" | estimate == "mu") %>%
  group_by(estimate) %>%
  mutate(id.vec = 1:n()) %>%
  group_by(year, site, id.vec, id) %>%
  pivot_wider(names_from = estimate, values_from = prediction) %>%
  ungroup() %>%
  mutate(diff = null-mu, 
         year = as.integer(year))

# lm1 <- lm(mu ~ null, data = temp)
# summary(lm1)


formerge <- s %>% 
  mutate(avg_lsize = map_dbl(lob.mass, ~{mean(.x$mass)}), 
         avg_usize = map_dbl(urc.mass, ~{mean(.x$mass)}), 
         diff = avg_lsize - avg_usize) %>%
  select(year, site, lob.a, urc.a, id, avg_lsize, avg_usize)

# 
# temp$res <- residuals(lm1)

temp <- temp %>% left_join(formerge, by = c("id", "year", "site"))

lm2 <- lm(diff ~ lob.a + urc.a + avg_lsize + avg_usize, temp)
summary(lm2)
 
lmer <- lmer(diff ~ scale(urc.a) + scale(avg_lsize) * scale(avg_usize) + (1|site) + (1|year), temp)
summary(lmer)

plot(diff ~ avg_lsize, temp)

plot(mu ~ null, data = temp)
abline(a = 0, b = 1)


temp <- s %>% ungroup() %>%
mutate(lob.samples = purrr::map_dbl(lob.mass, ~ mean(.x$mass)), 
       urc.samples = purrr::map_dbl(urc.mass, ~ mean(.x$mass)), 
       lob.max = purrr::map_dbl(lob.mass, ~max(.x$mass)))
summary(temp[,c("lob.samples", "urc.samples")])










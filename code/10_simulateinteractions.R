#-----------------------------------------------------------------
## Simulation function
#-----------------------------------------------------------------


post.a <- read.csv(here::here("data/cleaned/posteriors", "posteriors_posthoc_a.csv"))
post.h <- read.csv(here::here("data/cleaned/posteriors", "posteriors_posthoc_h.csv"))

meta <- read.csv(here::here("data/", "lob-metadata.csv"))


predict.fun <- function(lob.samples, urc.samples, urc_density, lob_density, beta1a. = post.a$beta1, beta2a. = post.a$beta2, beta1h. = post.h$beta1, beta2h. = post.h$beta2, h0. = post.h$alpha, a0. = post.a$alpha, T = 1, ...){
  
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


#-----------------------------------------------------------------------------------------------
## Simulate interactions
#-----------------------------------------------------------------------------------------------

## Mean allometric scaling

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


## Allometric scaling conditional on observed body size distributions

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
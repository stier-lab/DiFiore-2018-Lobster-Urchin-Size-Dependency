#----------------------------------------------------------------------------
## Set up and functions
#----------------------------------------------------------------------------
source(here::here("code", "1_setup.R"))
source(here::here("code", "10a_clean-obsdata.R"))
source(here::here("code", "theme.R"))


# Add line for merge conflict resolution.

# Misc. numbers needed later in code
      # Tank size
      
      # Full tank is 137 x 76. I'm going to stick with areal measures (rather than volumetric) because LTER works off of an areal basis. Therefore experimental tanks were 68.5 x 76. 
      
      tsize <- 137/2 * 76 /10000
      
      names <- s %>% mutate(id = paste(year, site, sep = "-")) %>%
        ungroup() %>%
        select(id) %>%
        distinct(id)
      
      names <- as.vector(names$id)
      

      
      plain <- function(x,...) {
        format(x, ..., scientific = FALSE, trim = TRUE, drop0trailing = T)
      }
      
#---------------------------------
## Data 
#---------------------------------
      
      ndraws = 10000

# Bring in the posteriors from the post hoc analysis of a and h ~ body size
# post.a <- read.csv(here::here("data/cleaned/posteriors", "posteriors_posthoc_a.csv")) %>% sample_draws(10000)
# post.h <- read.csv(here::here("data/cleaned/posteriors", "posteriors_posthoc_h.csv")) %>% sample_draws(10000)
      
      # # eliminate parameter uncertainty
      post.a <- read.csv(here::here("data/cleaned/posteriors", "allometric_populationSTAN.csv")) %>%
        as_tibble() %>% 
        pivot_wider(names_from = .variable, values_from = .value) %>%
        summarize(alpha = median(alphaa), beta1 = median(beta1a), beta2 = median(beta2a))
      
      post.h <- read.csv(here::here("data/cleaned/posteriors", "allometric_populationSTAN.csv")) %>%
        as_tibble() %>% 
        pivot_wider(names_from = .variable, values_from = .value) %>%
        summarize(alpha = median(alphah), beta1 = median(beta1h), beta2 = median(beta2h))

s # This dataframe is a tibble with nested body size distributions for lobsters and urchins.

set.seed(112615)
r.s <- s %>% # reproducible random draws from the size frequency distribution
  group_by(year, site) %>%
  mutate(urc_mass = purrr::map(urc.mass, sample_n, ndraws, replace = T), 
         lob_mass = purrr::map(lob.mass, sample_n, ndraws, replace = T)) %>% 
  dplyr::select(-urc.mass, -lob.mass)

#----------------------------------------------------------------------------
## Simuations based on experimental data
#----------------------------------------------------------------------------

# This is a function for the allometric functional response:

allometricFR <- function(lob_mass, urc_mass, urc_density, lob_density, beta1a., beta2a., beta1h., beta2h., h0., a0., T = 1, ...){
  
  loga <- a0. + beta1a.*log(lob_mass) + beta2a.*log(urc_mass)
  logh <- h0. + beta1h.*log(lob_mass) + beta2h.*log(urc_mass)
  a <- exp(loga) * tsize * 2 # convert units from per arena per 2 days to per m2 per day 
  h <- exp(logh) / 2 # convert units from 2 day trials to 1 day. 
  a*urc_density*lob_density*T / (1 + a*h*urc_density)
  
}

# Simulate interactions assuming all sources of uncertainty
full <- r.s %>%
  purrr::pmap(allometricFR,
              a0. = post.a$alpha,
              h0. = post.h$alpha,
              beta1a. = post.a$beta1,
              beta2a. = post.a$beta2,
              beta1h. = post.h$beta1,
              beta2h. = post.h$beta2) %>%
  purrr::flatten() %>%
  set_names(names) %>%
  as_tibble() %>%
  gather(id, prediction) %>%
  separate(id, into = c("year", "site"), sep = "[-]") %>%
  mutate(prediction = prediction,
         estimate = "full")

        # Summary stats for paper
        
        summary(full$prediction)
        rethinking::PI(full$prediction, 0.95)
        
        # CV spatial
        full %>%
          group_by(year) %>%
          summarize(cv_Xsite = sd(prediction)/mean(prediction)) %>%
          summarize(mean_cvXsite = mean(cv_Xsite), 
                    sd_cvXsite = sd(cv_Xsite))
        
        # CV temporal
        full %>%
          group_by(site) %>%
          summarize(cv_Xyear = sd(prediction)/mean(prediction)) %>%
          summarize(mean_cvXyear = mean(cv_Xyear), 
                    sd_cvXyear = sd(cv_Xyear))
        
        # CV within site/years
        full %>%
          group_by(year, site) %>%
          summarize(cv_Xboth = sd(prediction)/mean(prediction)) %>%
          ungroup() %>%
          summarize(mean_cvXboth = mean(cv_Xboth), 
                    sd_cvXboth = sd(cv_Xboth))


#-------------------------------------------------------------------------
## Estimates from the literature
#-------------------------------------------------------------------------

#-----------------------------        
# Rall
    # Rall et al. 2012 incorporate temp as a covariate. The basic regression equations are: 
        # log(a) = a0 + beta1*log(Mc) + beta2*log(Mr) + beta3*(T-T0)/k*T*T0
        # log(h) = h0 + beta1*log(Mc) + beta2*log(Mr) + beta3*(T-T0)/k*T*T0
    # where T is temp in kelvin, T0 = 293.15, mass are mg, and time is seconds, k is the boltzman's constant
        
bz_ev_k = 8.61733326 * 10^-5
        
allometricRall <- function(lob_mass, urc_mass, urc_density, lob_density, beta1a., beta2a., beta1h., beta2h., h0., a0., beta3a., beta3h., T = 1, ...){
  
  urc_mass_mg = urc_mass*1000
  lob_mass_mg = lob_mass*1000
  temp = (289.25 - 293.15)/(bz_ev_k*289.25*293.15)
  loga <- a0. + beta1a.*log(lob_mass_mg) + beta2a.*log(urc_mass_mg) + beta3a.*temp
  logh <- h0. + beta1h.*log(lob_mass_mg) + beta2h.*log(urc_mass_mg) + beta3h.*temp
  a <- exp(loga) 
  h <- exp(logh) 
  consumption_mg_m2_s = a*urc_density*lob_density*T / (1 + a*h*urc_density)
  consumption_g_m2_day = consumption_mg_m2_s * (60*60*24) / 1000 # convert mg per m2 per second to g per m2 per day.
  consumption_g_m2_day
  
}

rall <- r.s %>%
  purrr::pmap(allometricRall, 
              beta1a. = 0.85, 
              beta2a. = 0.09, 
              beta1h. = -0.76, 
              beta2h. = 0.76, 
              h0. = 10.38, 
              a0. = -21.23, 
              beta3a. = 0.42, 
              beta3h. = -0.30) %>% 
  purrr::flatten() %>%
  set_names(names) %>%
  as_tibble() %>%
  gather(id, prediction) %>%
  separate(id, into = c("year", "site"), sep = "[-]") %>%
  mutate(prediction = prediction,
         estimate = "rall")

summary(rall$prediction)        
hist(rall$prediction)        

#-----------------------------
        
        
        
        
#-----------------------------
# Barrios-Oneil
    # BO et al. 2019 use a different model structure with a temperature covariate
    # basic model structure: 
        # log(a) = a_0 + beta1*I_ACstat + beta2*I_F + beta3*log(Mc) + beta4*log(Mr) + beta5*log(T) + beta6*log(Mr)*log(T)
        # log(cmax) = C_0 + beta1*I_ACstat + beta2*I_F + beta3*log(Mc) + beta4*log(Mr) + beta5*log(T) + beta6*log(Mc)*log(T)
    # I assume that I_F = 0 and I_ACstat = 1, such that we have the regression equation for the active predators preying on static prey.
    # All beta values from Fig. 4 in text and Table S5. Intercepts estimated as the population level intercept + adjustment for active predators foraging on static prey (AC_stat) + the random intercept deviation from the mean for crustaceans. 
    # Units in the paper at g, days, and m^2. So no need for unit conversions.
        
allometricBO <- function(lob_mass, urc_mass, urc_density, lob_density, beta1a., beta2a., beta1h., beta2h., h0., a0., T = 1, temp, beta5a., beta6a., beta5h., beta6h., ...){
  
  loga <- a0. + beta1a.*log(lob_mass) + beta2a.*log(urc_mass) + beta5a.*log(temp) + beta6a.*log(urc_mass)*log(temp)
  logC <- h0. + beta1h.*log(lob_mass) + beta2h.*log(urc_mass) + beta5h.*log(temp) + beta6h.*log(lob_mass)*log(temp)
  a <- exp(loga)
  C <- exp(logC)
  h <- 1/C
  a*urc_density*lob_density*T / (1 + a*h*urc_density)
  
}
        
BO <- r.s %>%
  purrr::pmap(allometricBO, 
              beta1a. = 0.58, 
              beta2a. = 0.59, 
              beta1h. = 1.44, 
              beta2h. = 0.27,
              temp = 16.1, 
              beta5a. = 1.84,
              beta6a. = 0.12,
              beta5h. = 1.78, 
              beta6h. = -0.39,
              h0. = -6.07 + 0.55 + 0.48, 
              a0. = -8.08 + -1.07 + 0.48) %>% 
  purrr::flatten() %>%
  set_names(names) %>%
  as_tibble() %>%
  gather(id, prediction) %>%
  separate(id, into = c("year", "site"), sep = "[-]") %>%
  mutate(estimate = "BO")
#-----------------------------------------


#-----------------------------------------
# Uiterwall and Delong
  # UD 2020 account for arena size, dimensionality, and temp. The model structures are:
    # log(a) = a0 + beta1*T + beta2*T^2 + beta3*Mc + beta4*Mr + beta5*arena
    # log(h) = h0 + beta1*T + beta2*T^2 + beta3*Mc + beta4*Mr + beta5*Dim + beta6*arena
  # Values for all beta values and intercepts are from Table 1 and Table 2 in the main text of the manuscript.
  # I assume that temperature = the mean temperature over our experiment and arena size = our arena size.
  # Units in the manuscript are mg, cm^2, and days

allometricUD <- function(lob_mass, urc_mass, urc_density, lob_density, beta1a., beta2a., beta1h., beta2h., h0., a0., beta3a., beta3h., beta4a., beta4h., beta5a., beta5h., beta6h., temp, arena_m2, dim, T = 1, ...){
  lob_mass_mg = lob_mass*1000 # convert g to mg
  urc_mass_mg = urc_mass*1000 # convert g to mg
  arena_cm2 = arena_m2 * 10^4 # convert m2 to cm2
  urc_density_cm2 = urc_density / 10^4 # convert ind./ m2 to ind./ cm2
  lob_density_cm2 = lob_density / 10^4 # convert ind./m2 to ind./cm2
  loga <- a0. + beta1a.*log(lob_mass_mg) + beta2a.*log(urc_mass_mg) + beta3a.*temp + beta4a.*temp^2 + beta5a.*log(arena_cm2)
  logh <- h0. + beta1h.*log(lob_mass_mg) + beta2h.*log(urc_mass_mg) + beta3h.*temp + beta4h.*temp^2 + beta5h.*log(arena_cm2) + beta6h.*dim 
  a <- exp(loga) # per cm2 per day
  h <- exp(logh) # days
  consumption_mg_cm2_d = a*urc_density_cm2*lob_density_cm2*T / (1 + a*h*urc_density_cm2)
  consumption_g_cm2_d = consumption_mg_cm2_d / 1000 # convert mg to g
  consumption_g_m2_d = consumption_g_cm2_d * 10^4 # convert g per cm2 per day to g per m2 per day
  consumption_g_m2_d
}


UD <- r.s %>%
  purrr::pmap(allometricUD, 
              beta1a. = 0.05, 
              beta2a. = -0.0005, 
              beta1h. = -0.25, 
              beta2h. = 0.34, 
              h0. = 0.83, 
              a0. = -8.45, 
              beta3a. = 0.1, 
              beta3h. = -0.24, 
              beta4a. = -0.003,
              beta4h. = 0.005,
              beta5a. = 0.98, 
              beta5h. = -0.01, 
              beta6h. = -0.64, 
              temp = 16.1, 
              arena_m2 = tsize, 
              dim = 2) %>% 
  purrr::flatten() %>%
  set_names(names) %>%
  as_tibble() %>%
  gather(id, prediction) %>%
  separate(id, into = c("year", "site"), sep = "[-]") %>%
  mutate(estimate = "UD")

summary(UD$prediction)
hist(UD$prediction)

#-------------------------------------------------------------------------
## Comparison of variation due to body size and density.
#-------------------------------------------------------------------------

# Here I fix density to the regional averages. Therefore, the variance in IS estimated by this code will represent the variance due to variation in BODY SIZE not density.
full_meandensity <- r.s %>%
  ungroup() %>%
  mutate(urc_density = mean(urc_density),
         lob_density = mean(lob_density)) %>%
  # mutate(urc_density = 1, 
  #        lob_density = 1) %>%
  purrr::pmap(allometricFR, 
              a0. = post.a$alpha, 
              h0. = post.h$alpha, 
              beta1a. = post.a$beta1, 
              beta2a. = post.a$beta2, 
              beta1h. = post.h$beta1, 
              beta2h. = post.h$beta2) %>% 
  purrr::flatten() %>%
  set_names(names) %>%
  as_tibble() %>%
  gather(id, prediction) %>%
  separate(id, into = c("year", "site"), sep = "[-]") %>%
  mutate(prediction = prediction/2/tsize,
         estimate = "full_meandensity")

# Fix body size to the regional averages. Therefore, the variance in IS represents the variance due to variation in DENSITY not body size.

full_meanbodysize <- r.s %>%
  ungroup() %>%
  select(-urc_mass, -lob_mass) %>% 
  mutate(urc_mass = purrr::map(mean_urc_mass, rep, ndraws), 
         lob_mass = purrr::map(mean_lob_mass, rep, ndraws))%>%
  purrr::pmap(allometricFR, 
              a0. = post.a$alpha, 
              h0. = post.h$alpha, 
              beta1a. = post.a$beta1, 
              beta2a. = post.a$beta2, 
              beta1h. = post.h$beta1, 
              beta2h. = post.h$beta2) %>% 
  set_names(names) %>%
  as_tibble() %>%
  gather(id, prediction) %>%
  separate(id, into = c("year", "site"), sep = "[-]") %>%
  mutate(prediction = prediction/2/tsize,
         estimate = "full_meanbodysize")


df2 <- rbind(full, full_meandensity, full_meanbodysize)

df2 %>% filter(estimate != "full_meandensity" ) %>% ggplot()+
  geom_histogram(aes(x = prediction, ..ncount.., fill = estimate), position = "dodge")+
  labs(x = expression(paste("Interaction strength (ind. m"^-2,"d"^-1,")")), y = "Count")+
  scale_x_log10(labels = plain)+
  theme_bd()


ggplot(full)+
  geom_histogram(aes(x = prediction, ..ncount..), fill = "white", color = "gray50", bins = 30)+
  labs(x = expression(paste("Interaction strength (ind. m"^-2,"d"^-1,")")), y = "Count")+
  scale_x_log10(labels = plain)+
  theme_bd()+
  theme(legend.position = c(0.25, 0.75), text = element_text(size = 18, color = "lightgray"), axis.line = element_line(color = "lightgray"), axis.ticks = element_line(color = "lightgray"))+
  theme(
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill='transparent'),
    legend.box.background = element_rect(fill='transparent')
  )

ggsave("figures/esahisto1.png", width = 6, height = 10/3, bg = "transparent")

ggplot(full)+
  geom_histogram(aes(x = prediction, ..ncount..), fill = "white", color = "gray50", bins = 30)+
  geom_histogram(data = full_meanbodysize, aes(x = prediction, ..ncount..), fill = "#4472C4", color = "gray50", alpha = 0.5, bins = 30)+
  labs(x = expression(paste("Interaction strength (ind. m"^-2,"d"^-1,")")), y = "Count")+
  scale_x_log10(labels = plain)+
  theme_bd()+
  theme(legend.position = c(0.25, 0.75), text = element_text(size = 18, color = "lightgray"), axis.line = element_line(color = "lightgray"), axis.ticks = element_line(color = "lightgray"))+
  theme(
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill='transparent'),
    legend.box.background = element_rect(fill='transparent')
  )

ggsave("figures/esahisto2.png", width = 6, height = 10/3, bg = "transparent")


p2 <- ggplot(full)+
  geom_density(aes(x = prediction, ..scaled.., fill = "Total variation"), alpha = 0.5, adjust = 1.25)+
  geom_histogram(data = full_meanbodysize, aes(x = prediction, ..ncount.., fill = "Variation due\nto density"), position = "identity", color = "black")+
  scale_fill_manual(values = c("#AEBFA8","gray90"))+
  labs(x = expression(paste("Interaction strength (ind. m"^-2,"d"^-1,")")), y = "Scaled density/count", fill = "")+
  scale_x_log10(labels = plain)+
  coord_cartesian(ylim = c(0, 1))+
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1))+
  theme_bd()+
  theme(legend.position = c(0.2, 0.8))


#----------------------
## Combine and plot
#----------------------


sum <- data.frame(estimate = c("mean_full", "median_full", "Rall", "UD", "median_BO", "mean_BO", "mean_density", "mean_bodysize"), 
                  prediction = c(mean(full$prediction), 
                                 median(full$prediction),
                                 mean(rall$prediction), 
                                 mean(UD$prediction), 
                                 median(BO$prediction),
                                 mean(BO$prediction),
                                 mean(full_meandensity$prediction), 
                                 mean(full_meanbodysize$prediction)))

sum2 <- filter(sum, estimate %in% c("mean_full", "Rall", "UD", "mean_BO"))

# Summary stats for paper

# based on means
# initial - final / initial
((1.149810e-02 - 3.283366e-03) / 3.283366e-03)*100

1.149810e-02 / 3.283366e-03


((3.283366e-03 - 1.149810e-02) / 3.283366e-03)*100

# based on medians
((0.005907806 - 0.001862737) /  0.001862737) *100

0.005907806 / 0.001862737

library(calecopal)
#chaparal1
cols <- cal_palette("chaparral3", n = 4)
calecopal::chaparal1

histo <- ggplot(full)+
  geom_histogram(aes(x = prediction, ..ncount..), alpha = 0.1, color = "black")+
  geom_vline(data = sum2, aes(xintercept = prediction, color = estimate, linetype = estimate), show.legend = F, lwd = 2, linetype = c(1,2,3,4), color = c(cols[4], cols[2], cols[1], cols[3]))+
  scale_x_log10(labels = plain)+
  annotate('text', x = c(0.0375, 0.00025, 0.00001, 0.001),
           y = c(1.05, 0.6, 0.25, 1), label = c(
    "bar(italic(f(N,P,m[c],m[r])))",
    "bold(Marine)", 
    "bold(Cross~taxa)", 
    "bold(Mobile)"), parse = T, color = c(cols[4], cols[2], cols[1], cols[3]))+
  annotate('text', x = c(0.00025, 0.001), y = c(0.6-0.05, 1-0.05), label = c("bold(invertebrates)", "bold(crustaceans)"), color = c(cols[2], cols[3]), parse = T)+
  labs(x = expression(paste("Interaction strength (ind. m"^-2,"d"^-1,")")), y = "Scaled count")+
  coord_cartesian(ylim = c(0, 1.1))+
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1))+
  theme_bd()


#--------------------------------------------------------------
## Rank order figure
#--------------------------------------------------------------

# The idea here is to prove that estimates from the literature provide qualitatively accurate predictions of which site will have high or low interaction strength. 

s1 <- full %>% 
  mutate(id = paste(site, year, sep = "-")) %>%
  bind_rows(rbind(BO, rall, UD) %>% mutate(id = paste(site, year, sep = "-")) ) %>% 
  group_by(estimate, id) %>% 
  tidybayes::median_qi(prediction) %>%
  mutate(id = forcats::fct_reorder(id, filter(., estimate == "full") %>% pull(prediction))) %>%
  separate(id, into = c("site", "year"), sep = "-", remove = F) %>%
  mutate(character = case_when(estimate == "full" ~ "Experimental prediction", 
                               estimate == "BO" ~ "Marine crustaceans", 
                               estimate == "UD" ~ "Cross taxa", 
                               estimate == "rall" ~ "Marine invertebrates")) %>%
  ggplot(aes(x = prediction, y = id))+
  tidybayes::geom_pointinterval(aes(x = prediction, xmin = .lower, xmax =.upper, color = character), position = "dodge")+
  scale_x_log10()+
  facet_wrap(~site, scales = "free_y")+
  labs(y = "", x = expression(paste("Interaction strength (ind. m"^-2,"h"^-1,")")), color = "")+
  theme_bd()+
  theme(legend.position = c(0.8,0.2))

ggsave("figures/sup_fig-qualitative.png", plot = s1)

temp <- full %>% 
  mutate(id = paste(site, year, sep = "-")) %>%
  bind_rows(rbind(BO, rall, UD) %>% mutate(id = paste(site, year, sep = "-")) ) %>% 
  group_by(site, year, estimate, id) %>% 
  tidybayes::median_qi(prediction)%>%
  select(site, year, prediction, estimate) %>%
  pivot_wider(id_cols = c(site, year), names_from = estimate, values_from = prediction)


cor.test(temp$full, temp$BO, method = "spearman")
cor.test(temp$full, temp$UD, method = "spearman")
cor.test(temp$full, temp$rall, method = "spearman")


#--------------------------------------------------------------
## site to site comparison
#--------------------------------------------------------------

# Here I imaging a two panel scatter plot. To top is IS ~ mean pred size / mean prey size, the bottom is IS/pred density ~ urchin density. So what I need is to build a summary style data frame with one row per site/year.

df.sum <- s %>% 
  group_by(year, site) %>%
  left_join(full %>% mutate(year = as.integer(year))) %>% 
  nest(prediction = prediction) %>% 
  mutate(median_urc.mass = map_dbl(urc.mass, ~median(.x$mass)), 
         median_lob.mass = map_dbl(lob.mass, ~median(.x$mass)), 
         median_prediction = map_dbl(prediction, ~median(.x$prediction))) %>% 
  select(-c(urc.mass, lob.mass, estimate, prediction))

p1.r <- cor.test(df.sum$median_lob.mass/df.sum$median_urc.mass, df.sum$median_prediction/df.sum$lob_density)
p2.r <- cor.test(df.sum$urc_density, df.sum$median_prediction/df.sum$lob_density)

eq <- function(x) {
  value = as.numeric(x)
  as.character(
    as.expression(
      substitute(italic(r) == value)
    )
  )
}

p1 <- ggplot(df.sum, aes(y = median_prediction/lob_density, x = median_lob.mass/median_urc.mass))+
  geom_point(aes(size = urc_density))+
  labs(x = "Median predator:prey body mass ratio", y = expression(paste("Interaction strength (ind. m"^-2,"h"^-1,"P"^-1,")")), size = expression(paste("Urchin density (ind. m"^-2,")")))+
  annotate(x = 30, y = 0.25, geom="text", 
           label = eq(round(p1.r$estimate, 2)), parse = T, size = 5)+
  theme_bd()+
  theme(legend.position = c(0.2,0.8), text = element_text(size = 14))

p2 <- ggplot(df.sum, aes(y = median_prediction/lob_density, x = urc_density))+
  geom_point(aes(size = median_lob.mass/median_urc.mass ))+
  labs(size = "Median predator:prey\nbody mass ratio", y = expression(paste("Interaction strength (ind. m"^-2,"h"^-1,"P"^-1,")")), x = expression(paste("Urchin density (ind. m"^-2,")")))+
  annotate(x = 30, y = 0.25, geom="text", 
           label = eq(round(p2.r$estimate, 2)), parse = T, size = 5)+
  theme_bd()+
  theme(legend.position = c(0.2,0.8), text = element_text(size = 14))


s2 <- cowplot::plot_grid(p1, p2, nrow = 2, align = "v")
ggsave("figures/sup_fig-sitebysite.png", s2, width = 10*0.75, height = 10)








#-------------------------------------------------------------

# for SBC talk


# June 22, 2022: Idea, generate a time series plot showing AVERAGE/MEDIAN urchin density, lobster density, urchin size, lobster size, interaction strength, and kelp density. 


ggplot(full)+
  geom_histogram(aes(x = prediction), alpha = 0.1, color = "black")+
  scale_x_log10(labels = plain)+
  labs(x = expression(paste("Interaction strength (ind. m"^-2,"d"^-1,")")), y = "Count")+
  theme(text = element_text(size = 24), axis.text.x = element_text(size = 16))

ggsave(filename = "figures/histo-bland.png", device = "png")





full %>% 
  mutate(mpa = ifelse(site %in% c("IVEE", "NAPL"), "MPA", "Fished")) %>%
ggplot()+
  geom_histogram(aes(x = prediction, fill = mpa), alpha = 0.75)+
  scale_fill_manual(values = c("gray70", "gray10"))+
  scale_x_log10(labels = plain)+
  labs(x = expression(paste("Interaction strength (ind. m"^-2,"d"^-1,")")), y = "Count")+
  facet_wrap(~year, nrow = 2)+
  theme(text = element_text(size = 24))

ggsave(filename = "figures/hist-mpaXyear.png", device = "png", width = 13.33, height = 6 )


full %>% 
  mutate(mpa = ifelse(site %in% c("IVEE", "NAPL"), "MPA", "Fished")) %>%
  ggplot()+
  geom_histogram(aes(x = prediction, fill = mpa), color = "black", alpha = 0.5)+
  scale_fill_manual(values = c("gray70", "gray10"))+
  scale_x_log10(labels = plain)+
  labs(x = expression(paste("Interaction strength (ind. m"^-2,"d"^-1,")")), y = "Count")+
  theme(text = element_text(size = 24))

ggsave(filename = "figures/hist-mpa.png", device = "png")

temp <- full %>% 
  mutate(mpa = ifelse(site %in% c("IVEE", "NAPL"), "MPA", "Fished")) %>% 
  group_by(mpa, year) %>% 
  tidybayes::mean_qi(prediction)

temp <- temp %>% 
  group_by(mpa) %>%
  mutate(delta = prediction - first(prediction), 
         delta.low = .lower - first(.lower), 
         delta.high = .upper - first(.upper))


ggplot(temp, aes(x = as.numeric(year), y = prediction))+
  geom_line(aes(color = mpa))


ggplot(temp, aes(x = as.numeric(year), y = delta))+
  geom_line(aes(color = mpa))+
  geom_errorbar(aes(color = mpa, ymin = delta.low, ymax = delta.high))




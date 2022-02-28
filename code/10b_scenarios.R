#----------------------------------------------------------------------------
## Set up and functions
#----------------------------------------------------------------------------
source(here::here("code", "1_setup.R"))
source(here::here("code", "10a_clean-obsdata.R"))


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
      
      # This is a function for the allometric functional response:
      
      allometricFR <- function(lob_mass, urc_mass, urc_density, lob_density, beta1a., beta2a., beta1h., beta2h., h0., a0., T = 1, ...){
        
        loga <- a0. + beta1a.*log(lob_mass) + beta2a.*log(urc_mass)
        logh <- h0. + beta1h.*log(lob_mass) + beta2h.*log(urc_mass)
        a <- exp(loga)
        h <- exp(logh)
        a*urc_density*lob_density*T / (1 + a*h*urc_density)
        
      }
      
      # This is for Barrios-Oneil. Predcitions will be in units of # of prey eaten per m2 per day BUT they make predictions based on maximum consumption rate so we need to modify the formula
      allometricBO <- function(lob_mass, urc_mass, urc_density, lob_density, beta1a., beta2a., beta1h., beta2h., h0., a0., T = 1, ...){
        
        loga <- a0. + beta1a.*log(lob_mass) + beta2a.*log(urc_mass)
        logC <- h0. + beta1h.*log(lob_mass) + beta2h.*log(urc_mass)
        a <- exp(loga)
        C <- exp(logC)
        h <- 1/C
        a*urc_density*lob_density*T / (1 + a*h*urc_density)
        
      }
      
      # UNITS !!!!! All predictions of consumption rate should be in units of # of individual prey consumed per m2 per day!!!! 
      #- here we will make predictions based on the unit scale of the original data and then convert to consistent units afterwards
      
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
## Simuations
#----------------------------------------------------------------------------

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
  mutate(prediction = prediction/2/tsize, # get into units of # prey consumed per day per unit squared. The original time units are for 48 hours (or 2 days) therefore 1 day is /2
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
        

# Rall
rall <- r.s %>%
  purrr::pmap(allometricFR, 
              beta1a. = 0.85, 
              beta2a. = 0.09, 
              beta1h. = -0.76, 
              beta2h. = 0.76, 
              h0. = 10.38, 
              a0. = -21.23) %>% 
  purrr::flatten() %>%
  set_names(names) %>%
  as_tibble() %>%
  gather(id, prediction) %>%
  separate(id, into = c("year", "site"), sep = "[-]") %>%
  mutate(prediction = prediction*60*60*24,
         estimate = "rall")

# Barrios-Oneil
BO <- r.s %>%
  purrr::pmap(allometricBO, 
              beta1a. = 0.58, 
              beta2a. = 0.59, 
              beta1h. = 1.44, 
              beta2h. = 0.27, 
              h0. = -6.07 + 0.55 + 0.48, 
              a0. = -8.08 + -1.07) %>% 
  purrr::flatten() %>%
  set_names(names) %>%
  as_tibble() %>%
  gather(id, prediction) %>%
  separate(id, into = c("year", "site"), sep = "[-]") %>%
  mutate(estimate = "BO")

# Uiterwall and Delong
UD <- r.s %>%
  purrr::pmap(allometricFR, 
              beta1a. = 0.05, 
              beta2a. = -0.0005, 
              beta1h. = -0.25, 
              beta2h. = 0.34, 
              h0. = 0.83, 
              a0. = -8.45) %>% 
  purrr::flatten() %>%
  set_names(names) %>%
  as_tibble() %>%
  gather(id, prediction) %>%
  separate(id, into = c("year", "site"), sep = "[-]") %>%
  mutate(estimate = "UD")

# Here I fix density to the regional averages. Therefore, the variance in IS estimated by this code will represent the variance due to variation in BODY SIZE not density.
full_meandensity <- r.s %>%
  ungroup() %>%
  mutate(urc_density = max(urc_density),
         lob_density = max(lob_density)) %>%
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
  geom_histogram(aes(x = prediction, ..ncount..), alpha = 0.1, color = "gray50", bins = 100)+
  geom_histogram(data = full_meanbodysize, aes(x = prediction, ..ncount..), fill = "#F19B34", bins = 100)+
  labs(x = expression(paste("Interaction strength (ind. m"^-2,"d"^-1,")")), y = "Count")+
  scale_x_log10(labels = plain)+
  theme_bd()+
  theme(legend.position = c(0.25, 0.75), text = element_text(size = 18))


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


sum <- data.frame(estimate = c("mean_full", "median_full", "Rall", "UD", "BO", "mean_density", "mean_bodysize"), 
                  prediction = c(mean(full$prediction), 
                                 median(full$prediction),
                                 mean(rall$prediction), 
                                 mean(UD$prediction), 
                                 mean(BO$prediction),
                                 mean(full_meandensity$prediction), 
                                 mean(full_meanbodysize$prediction)))

sum2 <- filter(sum, estimate %in% c("mean_full", "Rall", "UD", "BO"))

# Summary stats for paper

sum2 %>% group_by(estimate) %>% 
  summarize(mean = mean(prediction))

# initial - final / initial
(0.0137 - 0.00328) / 0.00328 

0.0137 / 0.00328

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


timeseries <- full %>% 
  group_by(year, site) %>%
  mean_qi(prediction) %>% 
  ggplot(aes(x = as.numeric(year), y = prediction))+
  geom_line(aes(color = site))+
  geom_point(aes(color = site))+
  scale_color_manual(values = c('#AF8DC3','#C3AF8D', '#c3958d', '#8DC3AF', '#c38db5'))+
  labs(x = "", y = expression(paste("Interaction strength (ind. m"^-2,"h"^-1,")")), color = "Site")+
  theme_bd()+
  theme(legend.position = c(0.2, 0.75), text = element_text(size = 18), axis.text.y = element_text(size = 14))

p6b <- plot_grid(timeseries, histo, nrow = 1, align = "h", rel_widths = c(1,1.5))
ggsave(here::here("figures/", "observational_alt2.png"), p6b, width = 8.21429*1.5, height = 3*1.5)

fig5 <- plot_grid(histo, p2, nrow =1)

fig5


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
  pivot_wider(id_cols = c(site, year), names_from = estimate, values_from = prediction) %>%
  ungroup() %>%
  group_by(site) %>%
  purrr::map(~ cor.test(full, BO, method = "spearman"))


cor.test(temp$full, temp$BO, method = "spearman")
cor.test(temp$full, temp$UD, method = "spearman")
cor.test(temp$full, temp$rall, method = "spearman")


temp2 <- full %>% 
  mutate(id = paste(site, year, sep = "-")) %>%
  bind_rows(rbind(BO, rall, UD) %>% mutate(id = paste(site, year, sep = "-")) ) %>% 
  group_by(estimate) %>% 
  tidybayes::median_qi(prediction)

temp2$prediction[temp2$estimate == "full"] / temp2$prediction[temp2$estimate == "BO"]

3.003534*temp2$prediction[temp2$estimate == "BO"]

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




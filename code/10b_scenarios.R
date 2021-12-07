#----------------------------------------------------------------------------
## Set up and functions
#----------------------------------------------------------------------------
source(here::here("code", "1_setup.R"))
source(here::here("code", "10a_clean-obsdata.R"))


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

# Bring in the posteriors from the post hoc analysis of a and h ~ body size
# post.a <- read.csv(here::here("data/cleaned/posteriors", "posteriors_posthoc_a.csv")) %>% sample_draws(10000)
# post.h <- read.csv(here::here("data/cleaned/posteriors", "posteriors_posthoc_h.csv")) %>% sample_draws(10000)
      
      
      #TROUBLE SHOOTING1!!! 
      # # eliminate parameter uncertainty
      post.a <- read.csv(here::here("data/cleaned/posteriors", "posteriors_posthoc_a.csv")) %>% summarize(alpha = median(alpha), beta1 = median(beta1), beta2 = median(beta2))
      post.h <- read.csv(here::here("data/cleaned/posteriors", "posteriors_posthoc_h.csv"))%>% summarize(alpha = median(alpha), beta1 = median(beta1), beta2 = median(beta2))

s # This dataframe is a tibble with nested body size distributions for lobsters and urchins.

set.seed(112615)
r.s <- s %>% # reproducible random draws from the size frequency distribution
  group_by(year, site) %>%
  mutate(urc_mass = purrr::map(urc.mass, sample_n, 10000, replace = T), 
         lob_mass = purrr::map(lob.mass, sample_n, 10000, replace = T)) %>% 
  select(-urc.mass, -lob.mass)

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
  mutate(prediction = prediction*24/tsize,
         estimate = "full")


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
  mutate(prediction = prediction*24/tsize,
         estimate = "full_meandensity")

# Fix body size to the regional averages. Therefore, the variance in IS represents the variance due to variation in DENSITY not body size.
full_meanbodysize <- r.s %>%
  mutate(lob_mass = mean(s.avg$mean_lob_mass),
         urc_mass = mean(s.avg$mean_urc_mass)) %>%
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
  mutate(prediction = prediction*24/tsize,
         estimate = "full_meanbodysize")


df2 <- rbind(full, full_meandensity, full_meanbodysize)

histo2 <- df2 %>% filter(estimate != "full_meanbodysize" ) %>% ggplot()+
  geom_histogram(aes(x = prediction, fill = estimate))+
  labs(x = expression(paste("Interaction strength (ind. m"^-2,"d"^-1,")")), y = "Count")+
  scale_x_log10(labels = plain)+
  theme_classic()


df2 %>% filter(estimate != "full_meanbodysize" ) %>% ggplot()+
  geom_histogram(aes(x = prediction, fill = estimate))+
  labs(x = expression(paste("Interaction strength (ind. m"^-2,"d"^-1,")")), y = "Count")+
  scale_x_log10(labels = plain)+
  facet_wrap(~estimate)+
  theme_classic()

# calculate proportion of variance explained by body size
var_bodysize <- var(log(full_meandensity$prediction))
var_total <- var(log(full$prediction))
  
var_bodysize / var_total # this equal the proportion of the total variance due to body size!!!

# calcuate the proportion of variance explained by density
var_density <- var(log(full_meanbodysize$prediction)) 

var_density/var_total

1 - (var_bodysize / var_total) # this equals the proportion of the total variance due to density!!!










proposal.plot <- df2 %>% ggplot()+
  geom_density(aes(x = prediction, fill = estimate), alpha = 0.5)+
  labs(x = expression(paste("Interaction strength (ind. m"^-2,"d"^-1,")")), y = "Frequency", fill = "")+
  scale_x_log10(labels = plain)+
  scale_fill_discrete(labels = c("Variation in\nbody size and density", "Variation in\nbody size only"))+
  theme_classic()+
  theme(legend.position = c(0.25, 0.75), text = element_text(size = 18))

ggsave(here::here("figures/", "poposal_histogram.png"), proposal.plot)

proposal <- cowplot::plot_grid(NULL, timeseries, proposal.plot, histo, labels = "AUTO")
ggsave(here::here("figures", "proposal_ch1.png"), proposal, width = 8.5*1.75, height = 8.5*0.65*1.75)

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

library(calecopal)
#chaparal1
cal_palette("chaparral1")
calecopal::chaparal1

histo <- ggplot(full)+
  geom_histogram(aes(x = prediction), alpha = 0.1, color = "black")+
  geom_vline(data = sum2, aes(xintercept = prediction, color = estimate), show.legend = F, lwd = 1.5, linetype = 2, color = cal_palette("chaparral1", n = 4))+
  scale_x_log10(labels = plain)+
  annotate('text', x = sum2$prediction, y = 10000*c(4.75, 3.5, 3.5, 2), label = c(
    "bar(italic(f(mc,mr)))",
    "Rall", 
    "Uiterwaal", 
    "BarriosONeill"), parse=T)+
  labs(x = expression(paste("Interaction strength (ind. m"^-2,"d"^-1,")")), y = "Count")+
  theme_classic()+
  theme(text = element_text(size = 18))


timeseries <- df %>% 
  mutate(id = paste(year, site, sep = "")) %>%
  group_by(year, site, estimate, id) %>%
  mean_qi(prediction) %>% 
  filter(estimate == "full_wparameter") %>%
  ggplot(aes(x = as.numeric(year), y = prediction))+
  geom_line(aes(color = site))+
  geom_point(aes(color = site))+
  scale_color_manual(values = c('#AF8DC3','#C3AF8D', '#c3958d', '#8DC3AF', '#c38db5'))+
  labs(x = "", y = expression(paste("Interaction strength (ind. m"^-2,"h"^-1,")")), color = "Site")+
  theme_classic()+
  theme(legend.position = c(0.2, 0.75), text = element_text(size = 18), axis.text.y = element_text(size = 14))

p6b <- plot_grid(timeseries, histo, nrow = 1, align = "h", rel_widths = c(1,1.5))
ggsave(here::here("figures/", "observational_alt2.png"), p6b, width = 8.21429*1.5, height = 3*1.5)





















#-------------------------------------------------------------

# for SBC talk

ggplot(s.2a_full_wparameter)+
  geom_histogram(aes(x = prediction), alpha = 0.1, color = "black")+
  scale_x_log10(labels = plain)+
  labs(x = expression(paste("Interaction strength (ind. m"^-2,"d"^-1,")")), y = "Count")+
  theme(text = element_text(size = 24), axis.text.x = element_text(size = 16))

ggsave(filename = "figures/histo-bland.png", device = "png")





s.2a_full_wparameter %>% 
  mutate(mpa = ifelse(site %in% c("IVEE", "NAPL"), "MPA", "Fished")) %>%
ggplot()+
  geom_histogram(aes(x = prediction, fill = mpa), alpha = 0.75)+
  scale_fill_manual(values = c("gray70", "gray10"))+
  scale_x_log10(labels = plain)+
  labs(x = expression(paste("Interaction strength (ind. m"^-2,"d"^-1,")")), y = "Count")+
  facet_wrap(~year, nrow = 2)+
  theme(text = element_text(size = 24))

ggsave(filename = "figures/hist-mpaXyear.png", device = "png", width = 13.33, height = 6 )


s.2a_full_wparameter %>% 
  mutate(mpa = ifelse(site %in% c("IVEE", "NAPL"), "MPA", "Fished")) %>%
  ggplot()+
  geom_histogram(aes(x = prediction, fill = mpa), color = "black", alpha = 0.5)+
  scale_fill_manual(values = c("gray70", "gray10"))+
  scale_x_log10(labels = plain)+
  labs(x = expression(paste("Interaction strength (ind. m"^-2,"d"^-1,")")), y = "Count")+
  theme(text = element_text(size = 24))

ggsave(filename = "figures/hist-mpa.png", device = "png")

temp <- s.2a_full_wparameter %>% 
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






#---------------------------------------------------------------------------------------
## Scrap
#---------------------------------------------------------------------------------------

ggplot(df)+
  geom_density(aes(x = prediction, fill = estimate))+
  scale_x_log10()+
  geom_vline(xintercept = c(s.2c))+
  facet_wrap(~estimate, ncol = 1)

df %>% group_by(estimate) %>%
  mutate(.draw = 1:n()) %>%
  sample_draws(n = 10000) %>%
  ggplot()+
  stat_sample_slabinterval(aes(x = estimate, y = prediction), .width = c(0.9,0.05), slab_type = "histogram")+
  scale_y_log10()




ggplot(df[df$estimate == "full",])+
  geom_density(aes(x = prediction))+
  geom_vline(xintercept = c(s.2c, mean(df$prediction[df$estimate == "full"])))+
  scale_x_log10()


sum <- data.frame(estimate = c("mean_s.2a", "median_s.2a", "s.2c", "s.2d_Rall", "s.2d_UD", "s.2d_BO", "s.2e_unstructured"), 
                  prediction = c(mean(s.2a_full$prediction), 
                                 median(s.2a_full$prediction),
                                 s.2c, 
                                 mean(s.2d_rall$prediction), 
                                 mean(s.2d_UD$prediction), 
                                 mean(s.2d_BO$prediction), 
                                 mean(s.2e_unstructured$prediction)))


ggplot(df[df$estimate == "full",])+
  geom_density(aes(x = prediction))+
  geom_vline(data = filter(sum, estimate %in% c("mean_s.2a", "median_s.2a", "s.2c")), aes(xintercept = prediction, color = estimate))+
  scale_x_log10()



df %>% group_by(estimate) %>%
  mutate(.draw = 1:n()) %>%
  sample_draws(n = 10000) %>%
  ggplot()+
  #geom_histogram(data = df[df$estimate == "full",], aes(x = prediction))+
  stat_pointintervalh(aes(y = forcats::fct_relevel(estimate, "full", "BO", "rall", "UD", "unstructured"), x = prediction), width = 10000)+
  #stat_sample_slabinterval(data = df[df$estimate == "full"], aes(x = estimate, y = prediction), .width = c(0.9,0.05), slab_type = "histogram")+
  scale_x_log10()

df %>% group_by(estimate) %>%
  summarize(mean = mean(prediction), 
            median = median(prediction))


p1 <- df %>% filter(estimate != "full")
p2 <- df %>% filter(estimate == "full")

ggplot()+
  stat_sample_slabinterval(data = p2, aes(x= estimate, y = prediction))+
  stat_pointinterval(data = p1, aes(x = estimate, y = prediction))+
  scale_y_log10()




out1 <- allometricFR(lob_mass = 1000, 
                     urc_mass = 20, 
                     lob_density = 0.02, 
                     urc_density = 1:100,
                     beta1a. = 0.75, 
                     beta2a. = 0.33, 
                     beta1h. = -0.75, 
                     beta2h. = 0.5, 
                     h0. = -10, 
                     a0. = 4)
plot(1:100, out1)




out2 <- allometricFR(lob_mass = 1000, 
                     urc_mass = 20, 
                     lob_density = 0.02, 
                     urc_density = 1:100,
                     beta1a. = 0.75, 
                     beta2a. = 0.33, 
                     beta1h. = -0.75, 
                     beta2h. = 0.5, 
                     h0. = -5, 
                     a0. = 8)



plot(1:100, out1, ylim = c(0, 15000))
points(1:100, out2, col = "red")


plot(1:100, out2)


x <- 1:100
y <- 10 + 2*rnorm(length(x), x, 10)
plot(scale(y) ~ x)


summary(lm(y~x))
summary(lm(scale(y)~x))

lm1 <- lm(scale(y)~x)

plot(scale(y) ~ x)
abline(lm1)

plot(y~x)
abline(lm1)

plot(y ~x)
abline(lm(y~scale(x)))







# Ok, so if you fail to account for spatio-temporal variation in body size and density, you overestimate interaction strengths ~6 fold! See comparison of the mean of the fully integrated allometric model predictions and the prediciton of interaction strenght based on regional averages of body size and density!!!

# Interaction strength distributions at any site or time point are extremely left skewed. This means that there are lots of weak interactions and a few very strong interactions. Or in other words the majority of individual predators interact with the majority of prey extremely weakly! and only a few predator individuals interact with their prey strongly! This will have consequences on how we expect food webs to function! 

# On average, predictions of interaction strength based on published sources in the literature underestimate interaction strengths by an order of magnitude! But with higher taxonomic specificity the predictions are closer to our species-specific estimates.

# If we estimate interaction strength without accounting for body size (i.e. no allometric scaling of the FR), we will on average overestiamte interaction strength. 













































#----------------------------------------------------------------------------
## Contrast 1: How do predictions of interaction strength that incorporate body size differ from predictions that ignore body size? AND How accurate are predictions of interaction strength based on published allometric scaling relationships relative to our experimental results? 
#----------------------------------------------------------------------------

# For this contrast we are NOT interested in spatio-temporal variation in density or body size. Therefore, we will use the average densities and body sizes across the SBC dataset.

s.avg




# Scenario 1a. Allometric scaling based on experimental results

post <- read.csv(here::here("data/cleaned/posteriors", "allometric_population.csv"
)) %>% as_tibble() %>% sample_draws(10000) # predictions of consumption rate based on these posteriors will have units of # of prey consumed per hour per arena area (137/2 * 76 /10000)

s.1a <- allometricFR(lob_mass = s.avg$mean_lob_mass, 
                     urc_mass = s.avg$mean_urc_density, 
                     lob_density = s.avg$mean_lob_density, 
                     urc_density = s.avg$mean_urc_density,
                     beta1a. = post$beta1.a, 
                     beta2a. = post$beta2.a, 
                     beta1h. = post$beta1.h, 
                     beta2h. = post$beta2.h, 
                     h0. = post$mu.alpha.h, 
                     a0. = post$mu.alpha.a)

s.1a <- s.1a*24/tsize # convert estimates to units of # of prey consumed per day per m2
#s.1a <- s.1a/tsize

# Scenario 1b. Allometric scaling based on literature review

# Construct distributions based on Rall et al. 2012, Uitterwaal and DeLong 2020, and first principles

# This is for Rall et al. 2012. Precictions will be in units of # of prey eaten per m2 per second
s.1b.rall <- allometricFR(lob_mass = s.avg$mean_lob_mass, 
                          urc_mass = s.avg$mean_urc_density, 
                          lob_density = s.avg$mean_lob_density, 
                          urc_density = s.avg$mean_urc_density,
                          beta1a. = 0.85, 
                          beta2a. = 0.09, 
                          beta1h. = -0.76, 
                          beta2h. = 0.76, 
                          h0. = 10.38, 
                          a0. = -21.23)

s.1b.rall <- s.1b.rall*60*60*24

# This is for Barrios-Oneil. Predcitions will be in units of # of prey eaten per m2 per day BUT they make predictions based on maximum consumption rate so we need to modify the formula
allometricBO <- function(lob_mass, urc_mass, urc_density, lob_density, beta1a., beta2a., beta1h., beta2h., h0., a0., T = 1, ...){
  
  loga <- a0. + beta1a.*log(lob_mass) + beta2a.*log(urc_mass)
  logC <- h0. + beta1h.*log(lob_mass) + beta2h.*log(urc_mass)
  a <- exp(loga)
  C <- exp(logC)
  h <- 1/C
  a*urc_density*lob_density*T / (1 + a*h*urc_density)
  
}


s.1b.BO <- allometricBO(lob_mass = s.avg$mean_lob_mass, 
                        urc_mass = s.avg$mean_urc_density, 
                        lob_density = s.avg$mean_lob_density, 
                        urc_density = s.avg$mean_urc_density,
                        beta1a. = 0.58, 
                        beta2a. = 0.59, 
                        beta1h. = 1.44, 
                        beta2h. = 0.27, 
                        h0. = -6.07 + 0.55 + 0.48, 
                        a0. = -8.08 + -1.07)

# This is for Uiterwaal and DeLong 2020. Precictions will be in units of # of prey eaten per m2 per second
s.1b.UD <- allometricFR(lob_mass = s.avg$mean_lob_mass, 
                        urc_mass = s.avg$mean_urc_density, 
                        lob_density = s.avg$mean_lob_density, 
                        urc_density = s.avg$mean_urc_density,
                        beta1a. = 0.05, 
                        beta2a. = -0.0005, 
                        beta1h. = -0.25, 
                        beta2h. = 0.34, 
                        h0. = 0.83, 
                        a0. = -8.45)



# Scenario 1c. Ignore body size, a and h estimated from experimetnal data

df.null <- read.csv(here::here("data/cleaned/posteriors", "posteriors_null.csv")) %>% as_tibble() %>% sample_draws(10000) # units here will be in # of prey eaten per hour per experimental unit

s.1c <- allometricFR(lob_mass = s.avg$mean_lob_mass, 
                     urc_mass = s.avg$mean_urc_density, 
                     lob_density = s.avg$mean_lob_density, 
                     urc_density = s.avg$mean_urc_density,
                     beta1a. = 0, 
                     beta2a. = 0, 
                     beta1h. = 0, 
                     beta2h. = 0, 
                     h0. = log(df.null$h), 
                     a0. = log(df.null$a))

s.1c <- s.1c*24/tsize
#s.1c <- s.1c/tsize


plot(density(s.1a), type = "n", xlim = c(0, 0.1))
polygon(density(s.1a), col = "wheat")
#polygon(density(s.1c), col = "blue")
abline(v = c(s.1b.rall, s.1b.BO, s.1b.UD))
# Based on what we see here, there doesn't seem to be much of a difference between the unstructured model and the size structured model. I think the reason for this is that the unstructured model isn't truely unstructured because it is fit to data that was generated by the foraging of a gradient of predator sizes and prey sizes. So this sort of makes sense that we don't see major differences. What is worrying me is that 3 different estimates from the literature all give comparible estimates of interaction strength BUT our results are orders of magnitude higher. I'm worried that this has something to do with my unit conversions. 



plot(density(s.1c), type = "n", xlim = c(0, 0.004))
polygon(density(s.1a), col = "wheat")
polygon(density(s.1c), col = "blue")
abline(v = c(s.1b.rall, s.1b.BO, s.1b.UD), col = c("red", "pink", "green"))














#----------------------
## Scenario 2b

# A number for each site/year
#----------------------

s.2b_mu <- s %>%
  group_by(year, site) %>%
  mutate(lob_mass = purrr::map_dbl(lob.mass, ~ mean(.x$mass)), 
         urc_mass = purrr::map_dbl(urc.mass, ~ mean(.x$mass))) %>%
  purrr::pmap(allometricFR, 
              a0. = median(post$mu.alpha.a), 
              h0. = median(post$mu.alpha.h), 
              beta1a. = median(post$beta1.a), 
              beta2a. = median(post$beta2.a), 
              beta1h. = median(post$beta1.h), 
              beta2h. = median(post$beta2.h))%>% 
  set_names(names) %>%
  as_tibble() %>%
  gather(id, prediction) %>%
  separate(id, into = c("year", "site"), sep = "[-]") %>%
  mutate(estimate = "mu")




#----------------------
## Scenario 2e
#----------------------

s.2e_unstructured <- s %>%
  group_by(year, site) %>%
  mutate(urc_mass = purrr::map(urc.mass, sample_n, 10000, replace = T), 
         lob_mass = purrr::map(lob.mass, sample_n, 10000, replace = T)) %>%
  purrr::pmap(allometricFR, 
              a0. = log(median(df.null$a)), 
              h0. = log(median(df.null$h)), 
              beta1a. = 0, 
              beta2a. = 0, 
              beta1h. = 0, 
              beta2h. = 0) %>% 
  purrr::flatten() %>%
  set_names(names) %>%
  as_tibble() %>%
  gather(id, prediction) %>%
  separate(id, into = c("year", "site"), sep = "[-]") %>%
  mutate(prediction = prediction*24/tsize,
         estimate = "unstructured")

post.pop <- read.csv(here::here("data/cleaned/posteriors", "posteriors_population.csv"
)) %>% as_tibble() %>% sample_draws(10000)

FR <- function(a, h, urc_density, lob_density, T = 1, ...){
  a*urc_density*lob_density*T / (1 + a*h*urc_density)
  
}


s.nobodysize <- s %>%
  select(year, site, urc_density, lob_density) %>%
  group_by(year, site) %>%
  purrr::pmap(FR, 
              a = median(post.pop$mu.a), 
              h = median(post.pop$mu.h)) %>% 
  set_names(names) %>%
  as_tibble() %>%
  gather(id, prediction) %>%
  separate(id, into = c("year", "site"), sep = "[-]") %>%
  mutate(prediction = prediction*24/tsize,
         estimate = "nobodysize")



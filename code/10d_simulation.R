source("code/theme.R")
source("code/10b_scenarios.R")

set.seed(100024)

all_urcmass = unnest(s, cols = c(urc.mass))$mass
all_lobmass = unnest(s, cols = c(lob.mass))$mass

params <- list(lobdensity_mean = mean(r.s$lob_density), 
               lobdensity_25 = quantile(r.s$lob_density, 0.25), 
               urcdensity_mean = mean(r.s$urc_density), 
               urcdensity_25 = quantile(r.s$urc_density, 0.25), 
               lobmass_mean = mean(all_lobmass), 
               lobmass_25 = quantile(all_lobmass, 0.25), 
               urcmass_mean = mean(all_urcmass), 
               urcmass_25 = quantile(all_urcmass, 0.25))

#Attempted based just on the mean 

sim <- data.frame(site = c("initial", "initial", "final", "final"), 
                  group = c("A", "B", "A", "B"),
                  lob_density = rep(params$lobdensity_25, 4),
                  urc_density = c(params$urcdensity_25, 
                                  params$urcdensity_25,
                                  params$urcdensity_25*10, 
                                  params$urcdensity_25), 
                  lob_mass = c(params$lobmass_25, 
                               params$lobmass_25, 
                               params$lobmass_25,
                               params$lobmass_25*10), 
                  urc_mass = rep(params$urcmass_25, 4))


sim$IS <- allometricFR(lob_mass = sim$lob_mass, 
                       urc_mass = sim$urc_mass, 
                       urc_density = sim$urc_density, 
                       lob_density = sim$lob_density, 
                       a0. = post.a$alpha,
                       h0. = post.h$alpha,
                       beta1a. = post.a$beta1,
                       beta2a. = post.a$beta2,
                       beta1h. = post.h$beta1,
                       beta2h. = post.h$beta2)


sim %>% ggplot(aes(x = site, y = IS))+
  geom_point(aes(color =  site))+
  geom_line(aes(group = group))


# Attempted based on size frequency distributions


sim <- data.frame(site = c("Initial", "Initial", "Final", "Final"), 
                  group = c("A", "B", "A", "B"),
                  lob_density = rep(params$lobdensity_25, 4),
                  urc_density = c(params$urcdensity_25, 
                                  params$urcdensity_25,
                                  params$urcdensity_25*10, 
                                  params$urcdensity_25))

fitdist1 <- fitdistrplus::fitdist(all_urcmass, "gamma", method = "mle")

urcmass <- rgamma(n = sim$urc_density[1]*10000, shape = fitdist1$estimate[1], rate = fitdist1$estimate[2])
hist(urcmass)
summary(urcmass)
summary(all_urcmass)
start_urcmass <- expand.grid(group = c("A", "B"), urc_mass = urcmass) %>% mutate(site = "Initial")
end_urcmass <- data.frame(group = rep(c("A", "B"), each = length(urcmass)), 
                          site = rep("Final", 2*length(urcmass)), 
                          urc_mass = c(urcmass, urcmass))

forjoin_urc <- rbind(start_urcmass, end_urcmass)


fitdist2 <- fitdistrplus::fitdist(all_lobmass, "gamma", method = "mle")
lobmass <- rgamma(n = sim$lob_density[1]*10000, shape = fitdist2$estimate[1], rate = fitdist2$estimate[2])
hist(lobmass)
summary(lobmass)
summary(all_lobmass)
hist(all_lobmass)


start_lobmass <- expand.grid(group = c("A", "B"), lob_mass = lobmass) %>% mutate(site = "Initial")
end_lobmass <- data.frame(group = rep(c("A", "B"), each = length(lobmass)), 
                          site = rep("Final", 2*length(lobmass)), 
                          lob_mass = c(lobmass, lobmass*10))

forjoin_lob <- rbind(start_lobmass, end_lobmass)

df.sim <- sim %>% left_join(forjoin_urc) %>% nest(urc_mass = urc_mass) %>% 
  left_join(forjoin_lob) %>% nest(lob_mass = lob_mass)

names <- paste(sim$site, sim$group, sep = "-")

output <- df.sim %>% #reproducible random draws from the size frequency distribution
  group_by(site, group) %>%
  mutate(urc_mass = purrr::map(urc_mass, sample_n, ndraws, replace = T), 
         lob_mass = purrr::map(lob_mass, sample_n, ndraws, replace = T)) %>%
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
  separate(id, into = c("site", "group"), sep = "[-]") %>%
  mutate(prediction = prediction/2/tsize)

output %>% group_by(site, group) %>% 
  summarize(mean = mean(prediction), 
            median = median(prediction))


#Final - Initial / Initial

(0.00316 - 0.00312) / 0.00312
(0.01416 - 0.00389) / 0.00389

cal_palette("chaparral3", n = 5)

p3 <- output %>% 
  mutate(group = case_when(group == "A" ~ "10x increase\nin density", 
                           group == "B" ~ "10x increase\nin predator size")) %>%
  ggplot(aes(x = prediction, y = forcats::fct_rev(site)))+
  tidybayes::stat_slab(aes(fill = group), alpha = 0.5)+
  scale_fill_manual(values = c("#BED6B3", "#304818"))+
  tidybayes::stat_pointinterval(aes(group = group), .width = c(0))+
  scale_y_discrete(labels = c("", ""))+
  stat_summary(fun = "median", aes(linetype = group, group = group), geom = "line")+
  scale_linetype_manual(values = c(4,1))+
  coord_cartesian(xlim = c(0, 0.02))+
  labs(x = expression(paste("Interaction strength (ind. m"^-2,"d"^-1,")")), y = "Simulated change in predator biomass", fill = "", linetype = "")+
  annotate(geom = "text", y = c(1.5, 1.5), x = c(0.0025, 0.0070), label = c("+1.3%", "+264%"))+
  theme_bd()+
  theme(legend.position = c(0.8, 0.5))

fig5 <- plot_grid(p2, p3, histo, nrow = 3, align = "v", axis = "l")

ggsave("figures/fig5_histos.svg", fig5, device = "svg", width = 6, height = 10)
  



















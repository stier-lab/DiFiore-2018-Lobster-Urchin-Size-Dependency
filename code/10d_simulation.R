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


sim <- data.frame(site = c("initial", "initial", "final", "final"), 
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
start_urcmass <- expand.grid(group = c("A", "B"), urc_mass = urcmass) %>% mutate(site = "initial")
end_urcmass <- data.frame(group = rep(c("A", "B"), each = length(urcmass)), 
                          site = rep("final", 2*length(urcmass)), 
                          urc_mass = c(urcmass, urcmass))

forjoin_urc <- rbind(start_urcmass, end_urcmass)


fitdist2 <- fitdistrplus::fitdist(all_lobmass, "gamma", method = "mle")
lobmass <- rgamma(n = sim$lob_density[1]*10000, shape = fitdist2$estimate[1], rate = fitdist2$estimate[2])
hist(lobmass)
summary(lobmass)
summary(all_lobmass)
hist(all_lobmass)


start_lobmass <- expand.grid(group = c("A", "B"), lob_mass = lobmass) %>% mutate(site = "initial")
end_lobmass <- data.frame(group = rep(c("A", "B"), each = length(lobmass)), 
                          site = rep("final", 2*length(lobmass)), 
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


#final - initial / initial

(0.00316 - 0.00312) / 0.00312
(0.01416 - 0.00389) / 0.00389

cal_palette("chaparral2", n = 6)

p3 <- output %>% 
  mutate(group = case_when(group == "A" ~ "10x increase\nin density", 
                           group == "B" ~ "10x increase\nin predator size")) %>%
  ggplot(aes(x = forcats::fct_rev(site), y = prediction))+
  tidybayes::stat_slab(aes(fill = group, group = group), alpha = 0.5)+
  scale_fill_manual(values = c("#D98A63", "#A7C2CD"))+
  tidybayes::stat_pointinterval(aes(group = group), .width = c(0))+
  stat_summary(fun = "median", aes(linetype = group, group = group), geom = "line")+
  scale_linetype_manual(values = c(4,1))+
  coord_cartesian(ylim = c(0, 0.02))+
  labs(y = expression(paste("Interaction strength (ind. m"^-2,"d"^-1,")")), x = "", fill = "", linetype = "")+
  annotate(geom = "text", x = c(1.5, 1.5), y = c(0.0030, 0.0070), label = c("+1.3%", "+264%"))+
  theme_classic()+
  theme(text = element_text(size = 18), legend.position = c(0.2, 0.8))

toprow = plot_grid(p2, p3, nrow = 1)
bottomrow = plot_grid(NULL, histo, NULL, nrow = 1, rel_widths = c(0.15, 0.7, 0.15))



fig5 <- cowplot::plot_grid(toprow, bottomrow, nrow = 2)
ggsave(filename = "figures/figure5_histos.png", plot = fig5, width = 8.5*1.5, height = 8.5*0.65*1.5)















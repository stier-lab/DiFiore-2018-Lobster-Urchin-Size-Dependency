source(here::here("code/", "1_setup.R"))
source(here::here("code", "10a_clean-obsdata.R"))

df <- s %>% select(-c(urc_density, lob_density)) %>% pivot_longer(names_to = "species", values_to = "mass", -c(year, site)) %>% unnest(cols = c(mass))

scale2 <- function(x, na.rm = FALSE) (x - mean(x, na.rm = na.rm)) / sd(x, na.rm)

df %>% filter(site %in% c("MOHK", "NAPL")) %>%
  group_by(species) %>%
  mutate(mass = scale2(mass))%>%
  ggplot(aes(x = mass, y = as.factor(year), fill = species))+
  ggridges::stat_density_ridges(geom = "density_ridges_gradient", bandwidth = 0.35)+
  # geom_vline(aes(xintercept=median, col = reserve))+
  scale_fill_manual(values = c("#FEE08B", "#af8dc3"), labels = c("Lobster", "Urchin"))+
  geom_vline(xintercept = 0, col = "gray", lty = 4)+
  labs(x = "Standardized body mass", y = "Density", fill = "", y = "")+
  facet_wrap(~site)+
  theme_pubclean()

ggsave(filename = here::here("figures", "normalized_hist.png"), device = "png")


p1 <- df %>%
  ggplot(aes(x = mass)) +
  geom_density(aes(fill = species), adjust = 2, show.legend = F) +
  scale_fill_manual(values = c("#a13617", "#1782a1"))+
  facet_wrap(~species, scales = "free")+
  labs(x = "Body mass (g)", y = "Density")+ 
  theme(strip.background = element_blank(), strip.text = element_blank(), text = element_text(size = 24))

p2 <- s %>%
  select(year, site, urc_density, lob_density) %>%
  pivot_longer(cols = c(urc_density, lob_density), names_to = "species", values_to = "density") %>%
  ggplot(aes(x = density)) +
  geom_density(aes(fill = species), adjust = 1, show.legend = F) +
  scale_fill_manual(values = c("#a13617", "#1782a1"))+
  facet_wrap(~species, scales = "free")+
  labs(x = expression(paste("Abundance (ind. m"^-2,")")), y = "Density")+ 
  theme(strip.background = element_blank(), strip.text = element_blank(), text = element_text(size = 24))

cowplot::plot_grid(p1, p2, align = "hv", nrow = 2)
ggsave(filename = "figures/histo-byspecies.png", device = "png")

fct_relevel(temp$species, "lob.mass", "urc.mass")


barnes <- read.csv(file = "data/Predator_and_prey_body_sizes_in_marine_food_webs_vsn4.csv") %>% select(Predator.mass, Prey.mass, Predator, Prey)




allometricFR <- function(lob_mass, urc_mass, urc_density, lob_density, beta1a., beta2a., beta1h., beta2h., h0., a0., T = 1, ...){
  
  loga <- a0. + beta1a.*log(lob_mass) + beta2a.*log(urc_mass)
  logh <- h0. + beta1h.*log(lob_mass) + beta2h.*log(urc_mass)
  a <- exp(loga)
  h <- exp(logh)
  a*urc_density*lob_density*T / (1 + a*h*urc_density)
  
}

allometricMC <- function(lob_mass, urc_mass,beta1a., beta2a., beta1h., beta2h., h0., a0., T = 1, ...){
  
  logh <- h0. + beta1h.*log(lob_mass) + beta2h.*log(urc_mass)
  h <- exp(logh)
  1/h
  
}

barnes$maxC <- allometricMC(lob_mass = barnes$Predator.mass, urc_mass = barnes$Prey.mass, beta1a. = 0.05, beta2a. = -0.0005, beta1h. = -0.25, beta2h. = 0.34, h0. = 0.83, a0. =  -8.45)

plot(log(maxC) ~ log(Predator.mass), barnes)

summary(lmer(log(maxC)~ log(scale(Predator.mass)) + log(scale(Prey.mass)) + (1|Predator), barnes))

ggplot(barnes, aes(x = Predator.mass, y = maxC))+
  geom_point()+
  geom_smooth(method = "lm")+
  scale_y_log10()+
  scale_x_log10()


plot(y~x)

x <- c(1:30, seq(30, 80, length.out = 1000), 80:140)
y <- rnorm(n = length(x), mean = x, sd = 10)

plot(y ~ x)







forage <- read.csv("~/Downloads/resource_map_doi_10_5063_F17H1GTQ/data/FoRAGE_db_12_19_18_data_set.csv") %>% janitor::clean_names() %>%
  mutate(maxC = 1/fittted_h_day) %>%
  drop_na(predator_mass_mg, prey_mass_mg, temp_c)

lm1 <- lm(log(maxC) ~ log(predator_mass_mg) + log(prey_mass_mg) + temp_c + I(temp_c^2), forage)

car::crPlot(model = lm1, variable = "log(predator_mass_mg)")

lm1.1 <- lm(log(maxC) ~ log(predator_mass_mg) + temp_c + I(temp_c^2), forage)
lm1.2 <- lm(log(predator_mass_mg) ~ log(prey_mass_mg) + temp_c + I(temp_c^2), forage)

plot(residuals(lm1.1) ~ residuals(lm1.2) , forage)
summary(lm1)

forage$reslm1.1 <- residuals(lm1.1)
forage$reslm1.2 <- residuals(lm1.2)


col <- colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))

p1 <- forage %>% filter(exp(reslm1.2) > 10e-6, major_grouping_1 != "Rotifer") %>%
ggplot(aes(y = reslm1.1, x = exp(reslm1.2)))+
  geom_point(alpha = 0.75, size = 1.5)+
  scale_x_log10(breaks = c(0.0001, 1, 100, 1000, 10000), labels = c(0.0001, 1, 100, 1000, 10000))+
  labs(x = "Predator size (mg)", y = "Consumption rate")

ggsave(filename = "figures/pres-p1.png", plot = p1, device = "png")

p2 <- forage %>% filter(exp(reslm1.2) > 10e-6, major_grouping_1 != "Rotifer") %>%
  ggplot(aes(y = reslm1.1, x = exp(reslm1.2)))+
  geom_point(alpha = 0.75, size = 1.5)+
  geom_smooth(method = "lm", col = "darkred")+
  scale_x_log10(breaks = c(0.0001, 1, 100, 1000, 10000), labels = c(0.0001, 1, 100, 1000, 10000))+
  coord_cartesian(ylim = c(-10, 25))+
  labs(x = "Predator size (mg)", y = "Consumption rate")

ggsave(filename = "figures/pres-p2.png", plot = p2, device = "png")

p3 <- forage %>% filter(exp(reslm1.2) > 10e-6, major_grouping_1 != "Rotifer") %>%
  ggplot(aes(y = reslm1.1, x = exp(reslm1.2)))+
  geom_point(aes(color = major_grouping_1), show.legend = F, alpha = 0.75, size = 1.5)+
  scale_color_manual(values = col(19))+
  geom_smooth(method = "lm", aes(col = major_grouping_1), alpha = 0.25, show.legend = F)+
  geom_smooth(method = "lm", col = "darkred", se = F, lwd = 1.5)+
  scale_x_log10(breaks = c(0.0001, 1, 100, 1000, 10000), labels = c(0.0001, 1, 100, 1000, 10000))+
  coord_cartesian(ylim = c(-10, 25))+
  labs(x = "Predator size (mg)", y = "Consumption rate")

ggsave(filename = "figures/pres-p3.png", plot = p3, device = "png")


forgg <- forage %>% filter(exp(reslm1.2) > 10e-6, major_grouping_1 != "Rotifer")

bo <- data.frame(mc = forgg$predator_mass_mg, mr = forgg$prey_mass_mg)
bo$MaxC <- exp(-5.04) + 1.44*log(bo$mc) + 0.27*log(bo$mr)
bo$MaxC <- exp(bo$MaxC)

bo1.1 <- lm(log(MaxC) ~ log(mc), bo)
bo1.2 <- lm(log(mc) ~ log(mr), bo)

bo$bo1.1 <- residuals(bo1.1)
bo$bo1.2 <- residuals(bo1.2)

plot(log(MaxC) ~ log(mc), bo)
plot(residuals(bo1.1) ~ log(residuals(bo1.2)))
  
ggplot(forgg, aes(y = reslm1.1, x = exp(reslm1.2)))+
  geom_point(show.legend = F, alpha = 0.75, size = 1.5)+
  scale_color_manual(values = col(19))+
  geom_smooth(method = "lm", alpha = 0.25, show.legend = F, col = "darkred", lwd = 1.5)+
  geom_smooth(data = filter(forgg, major_grouping_1 == "Crustacean"), method = "lm", col = "darkred", se = F, lwd = 1.5)+
  geom_smooth(data = bo, aes(y = log(MaxC), x = mc), method  = "lm")+
  scale_x_log10(breaks = c(0.0001, 1, 100, 1000, 10000), labels = c(0.0001, 1, 100, 1000, 10000))+
  coord_cartesian(ylim = c(-10, 25))+
  labs(x = "Predator size (mg)", y = "Consumption rate")



lm1.1 <- lm(log(maxC) ~ log(predator_mass_mg)*major_grouping_1 + temp_c + I(temp_c^2), forage)
lm1.2 <- lm(log(predator_mass_mg) ~ log(prey_mass_mg) + temp_c + I(temp_c^2), forage)



x <- rnorm(10000, mean = 3, sd = 1)
y <- rnorm(10000, mean = 0, sd = 1)

gg <- data.frame(name = rep(c("x", "y"), each = 10000), value = c(x, y))

ggplot(gg)+
  geom_density(aes(x = value, fill = name), show.legend = F, alpha = 0.75)+
  scale_fill_manual(values = c("#7030A0", "#8B0000"))+
  labs(x = "", y = "")+
  theme(axis.text = element_blank(), 
        axis.line.y = element_blank(), 
        axis.ticks = element_blank())

ggsave(filename = "figures/pres-p4.png", device = "png", width = 3, height = 2)


#-------------------------------------------------------------------------------------------
## Scrap
#-------------------------------------------------------------------------------------------

df %>% filter(site %in% c("MOHK", "NAPL")) %>%
  # group_by(species) %>%
  # mutate(mass = scale(mass)) %>%
  ggplot(aes(x = mass, y = as.factor(year)))+
  ggridges::geom_density_ridges(aes(fill = species), scale = 2, alpha = 0.9)+
  # geom_vline(aes(xintercept=median, col = reserve))+
  scale_fill_manual(values = c("#FEE08B", "#af8dc3"))+
  # scale_x_log10()+
  # scale_color_manual(values = c("#7fbf7b", "#af8dc3"))+
  # labs(x = "Size (cm)", y = "Density", fill = "")+
  facet_wrap(~site)+
  theme_pubclean()


df %>% filter(site %in% c("MOHK", "NAPL")) %>%
  group_by(species, site) %>%
  mutate(mass = scale(mass)) %>%
  ggplot(aes(x = mass, y = as.factor(year), fill = species, height = ..density..))+
  ggridges::geom_density_ridges(scale = 1.5, alpha = 0.9, stat = "density", bandwidth = 0.01)+
  # geom_vline(aes(xintercept=median, col = reserve))+
  scale_fill_manual(values = c("#FEE08B", "#af8dc3"))+
  # scale_color_manual(values = c("#7fbf7b", "#af8dc3"))+
  # labs(x = "Size (cm)", y = "Density", fill = "")+
  facet_wrap(~site)+
  theme_pubclean()


df %>% filter(site %in% c("MOHK", "NAPL")) %>%
  # group_by(species, site) %>%
  # mutate(mass = scale(mass)) %>%
  ggplot(aes(x = mass))+
  geom_histogram(aes(fill = species))+
  # geom_vline(aes(xintercept=median, col = reserve))+
  scale_fill_manual(values = c("#FEE08B", "#af8dc3"))+
  # scale_color_manual(values = c("#7fbf7b", "#af8dc3"))+
  # labs(x = "Size (cm)", y = "Density", fill = "")+
  facet_wrap(~site)+
  theme_pubclean()


df %>% ggplot(aes(x = mass)) + 
  geom_histogram()+
  scale_x_log10()+
  facet_wrap(~species, scales = "free")

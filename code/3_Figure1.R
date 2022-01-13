source("code/1_setup.R")

#------------------------------------------------------------
## Clean and combine SBC LTER data
#------------------------------------------------------------

# Biomass conversion from: Nelson, C, D. Reed, S. Harrer, R. Miller. 2021. SBC LTER: Reef: Coefficients for estimating biomass from body size or percent cover for kelp forest species ver 3. Environmental Data Initiative. https://doi.org/10.6073/pasta/0fe9233dabe35df5d61fb3b07f8fb51e. Accessed 2022-01-04.



# get the lte urchin size data

# data citation: Reed, D, R. Miller. 2021. SBC LTER: Reef: Long-term experiment: Kelp removal: Urchin size frequency distribution ver 21. Environmental Data Initiative. https://doi.org/10.6073/pasta/fd564dddfe7b77fe9e4bd8417f166057. Accessed 2022-01-03.

urc.s <- read.csv(here("data/LTER", "LTE_Urchin_All_Years_20210209.csv"), header = T) %>% # get urchin size data for SBC LTER LTE transects.
  filter(TREATMENT == "CONTROL", COMMON_NAME == "Purple Urchin", YEAR > 2011) %>%
  dplyr::select(YEAR, MONTH, DATE, SITE, TRANSECT, SIZE, COUNT) %>% 
  na_if(-99999) %>%
  rename_all(tolower) %>%
  filter(month == 8) 

sites <- distinct(urc.s, site)$site #save sites to filter other data

urc.s <- urc.s %>% 
  group_by_at(vars(-count)) %>% # group by all variables other than count
  complete(count = full_seq(1:count, 1))%>% # add dummy variable that is a running count within a size bin
  ungroup() %>%
  mutate(size = as.numeric(size) + 0.25, # add adjustment so that the size is the center of the size bin
         count = NULL, 
         mass =  0.000592598*(size*10)^2.872636198*1.01) %>% # estimate mass based on SBC LTER test diameter - weigh relationship.  IS THIS WET OR DRY MASS!!!!!!
  drop_na(size) %>% 
  group_by(year, site) %>%
  dplyr::select(-c(size, month, date, transect)) 


# Organize and clean lobster data

# Data citation: Reed, D, R. Miller. 2021. SBC LTER: Reef: Abundance, size and fishing effort for California Spiny Lobster (Panulirus interruptus), ongoing since 2012 ver 6. Environmental Data Initiative. https://doi.org/10.6073/pasta/0bcdc7e8b22b8f2c1801085e8ca24d59. Accessed 2022-01-03.

lob <- read.csv(here("data/LTER", "Lobster_Abundance_All_Years_20210412.csv"), header = T, na.strings = c(-99999), stringsAsFactors = F) %>% # get lobster abundance and size data from LTER
  rename_all(tolower) %>%
  filter(site %in% sites) %>% 
  filter(size_mm != 600) %>% # get rid of the size error for one lobster
  dplyr::select(year, month, site, transect, replicate, size_mm, count)

lob.s <- lob %>% 
  group_by(year, site, size_mm) %>% 
  summarize(count = sum(count))%>%
  group_by_at(vars(-count)) %>% 
  complete(count = full_seq(1:count, 1)) %>%
  ungroup() %>%
  mutate(size_mm = as.numeric(size_mm), 
         count = NULL, 
         mass = 0.001352821*(size_mm)^2.913963113) %>% #IS THIS WET OR DRY MASS!!!!!!!!!
  drop_na(size_mm) %>% 
  group_by(year, site) %>%
  dplyr::select(-size_mm)


urc.s %>% mutate(species = "Urchin") %>%
  bind_rows(lob.s %>% mutate(species = "Lobster")) %>%
  group_by(year, site, species) %>%
  filter(site %in% c("MOHK", "NAPL")) %>%
  ggplot(aes(x = mass, y = as.factor(year)))+
  geom_density_ridges(aes(fill = species), alpha = 0.9)+
  # geom_vline(aes(xintercept=median, col = reserve))+
  scale_fill_manual(values = c("#FEE08B", "#af8dc3"))+
  # scale_color_manual(values = c("#7fbf7b", "#af8dc3"))+
  labs(x = "Size (cm)", y = "Density", fill = "")+
  facet_wrap(~site)+
  theme_pubclean()

urc.s %>% mutate(species = "Urchin") %>%
  bind_rows(lob.s %>% mutate(species = "Lobster")) %>%
  group_by(year, site, species) %>%
  filter(site %in% c("MOHK", "NAPL")) %>%
  median_qi(mass) %>%
  ggplot(aes(x = mass, y = as.factor(year)))+
  geom_pointinterval(aes(color = species, xmin = .lower, xmax = .upper ), position = position_dodge(width = 0.25)) + 
  facet_wrap(~site)
  

urc.s %>% mutate(species = "Urchin") %>%
  bind_rows(lob.s %>% mutate(species = "Lobster")) %>%
  group_by(species) %>%
  mutate(mass.s = scale(mass)) %>%
  group_by(year, site, species) %>%
  filter(site %in% c("MOHK", "NAPL")) %>%
  median_qi(mass.s) %>%
  ggplot(aes(x = mass.s, y = as.factor(year)))+
  geom_pointinterval(aes(color = species, xmin = .lower, xmax = .upper ), position = position_dodge(width = 0.25)) + 
  facet_wrap(~site)+
  theme_pubclean()


p1 <- urc.s %>% mutate(species = "Urchin") %>%
  bind_rows(lob.s %>% mutate(species = "Lobster")) %>%
  group_by(species) %>%
  mutate(mass.s = scale(mass)) %>%
  group_by(year, site, species) %>%
  filter(site %in% c("MOHK", "NAPL")) %>%
  mutate(site = case_when(site == "MOHK" ~ "Site 1", 
                          site == "NAPL" ~ "Site 2")) %>%
  ggplot(aes(x = mass.s, y = as.factor(year)))+
  geom_density_ridges(aes(fill = species), scale = 1.5, alpha = 0.9)+
  scale_fill_manual(values = c("#FEE08B", "#af8dc3"))+
  labs(x = "Scaled body mass", y = "Density", fill = "")+
  facet_wrap(~site)+
  theme_pubclean()+
  theme(strip.background = element_rect(fill= "transparent"), 
        strip.text = element_text(size = 18), 
        legend.position = "right",
        legend.text = element_text(size = 14), 
        axis.text.y = element_text(size = 12), 
        axis.title.y = element_text(size = 14,margin = margin(t = 0, r = 20, b = 0, l = 0)), 
        axis.title.x = element_text(size = 14,margin = margin(t = 20, r = 0, b = 0, l = 0)))


ggsave("figures/fig1_2sitecompare.png", p1, device = "png")




#------------------------------------
## Summary stats for P1 of results
#------------------------------------

summary(lob.s$mass)
rethinking::PI(lob.s$mass, prob = .95)
# 3%       98% 
#   88.84656 897.78328 

summary(urc.s$mass)
rethinking::PI(urc.s$mass, prob = .95)
# 3%        98% 
#   8.161326 132.176989 





























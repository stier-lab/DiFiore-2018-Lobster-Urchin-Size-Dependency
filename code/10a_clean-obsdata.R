source("code/1_setup.R")


#------------------------------------------------------------
## Clean and combine SBC LTER data
#------------------------------------------------------------

# Notes: The SBC LTER collects urchin size data at one transect at each site. Lobster size and abundance data is collected at >3 transects at each of the five sites. Urchin abundance is estimated at >3 transects per site. Here we will use the urchin size data from the transect on which it is recorded. We will use individual observations of lobster at each size for the size estimation, but will average lobster abundance across transects to have a single estimate of lobster abundance at the site level. Similarly, we will estimate the mean urchin abundance across transects. 

# get the lte urchin size data

urc.s <- read.csv(here("data/LTER", "LTE_Urchin_All_Years.csv"), header = T) %>% # get urchin size data for SBC LTER LTE transects.
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

mean_urc_mass <- mean(urc.s$mass)


urc.s <- urc.s %>%
  nest(urc.mass = mass)



# get the community data

urc.a <- read.csv(here("data/LTER", "Annual_All_Species_Biomass_at_transect.csv"), stringsAsFactors = F,na.strings ="-99999") %>%
  rename_all(tolower) %>%
  dplyr::select(year, month, site, transect, sp_code, density) %>%
  mutate(id = paste(site, transect, sep = "")) %>%
  filter(sp_code %in% c("SPL"), year > 2011, site %in% sites) %>% #only the urchin biomass density data for now..
  group_by(year, site) %>%
  summarize(urc_density = mean(density, na.rm = T))

mean_urc_density <- mean(urc.a$urc_density)

# Organize and clean lobster data

lob <- read.csv(here("data/LTER", "Lobster_Abundance_All_Years.csv"), header = T) %>% # get lobster abundance and size data from LTER
  na_if(-99999) %>% 
  rename_all(tolower) %>%
  filter(site %in% sites) %>% 
  dplyr::select(year, month, site, transect, replicate, size_mm, count)

lob.a <- lob %>%
  group_by(year, site, transect) %>%
  summarize(density = sum(count, na.rm = T)/1200) %>%
  group_by(year, site) %>%
  summarize(lob_density = mean(density, na.rm = T))

mean_lob_density <- mean(lob.a$lob_density)

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

mean_lob_mass <- mean(lob.s$mass)


lob.s <- lob.s %>%
  nest(lob.mass = mass)

df <- left_join(urc.a, lob.a)
df <- left_join(df, urc.s)
df <- left_join(df, lob.s)
s <- df



s.avg <- list(mean_lob_density = mean_lob_density,
              mean_urc_density = mean_urc_density,
              mean_lob_mass = mean_lob_mass, 
              mean_urc_mass = mean_urc_mass)











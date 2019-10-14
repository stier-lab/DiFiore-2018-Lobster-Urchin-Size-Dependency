#----------------------------------------------------------------------------------
#The goal of this code is to compile data on the size distribtions and abundances of lobsters and purple urchins at all sites. Each data set has its own ideosyncracies in how lobsters and urchins are monitored and sized, with methods shifting through time. The goal is to creat two data sets 1) a data set of lobster and urchin density at all monitored sites through time, and 2) a data set of lobster and urchin sizes at all monitored sites through time. As a caution not all sites were monitored in all years, and lobsters and urchins were not sized in each year. Furthermore, PISCO and LTER survey each species differently, i.e. different survey areas, and different sizing techniques.

# LTER data description
# The SBC LTER has monitored lobster and urchin abundance at ~11 sites since 2002. Only urchins > 2.5 cm are recorded in surveys. At the core sites lobster and urchin average size is estimated, however we will ignore the average sizes at a particularly transect at a particualar time because were interested in size frequency distributions NOT averages. Since 2012, the SBC LTER has monitored lobsters in targeted surveys conducted annually at 9 sites. Each lobster is sized visually by trained observers. We will focus on this data for lobster size distrubtions for the SBC LTER sites. Since 2008 the SBC LTER has estimated the sizes of purple urchins at 5 LTE (long term experiment) sites quarterly. We will focus on this data for urchin size distributions for the SBC LTER. 

# PISCO data description
# PISCO has monitored lobster and urchin abundance at ~30 sites since 2003. Only urchins > 2.5 cm are recorded in surveys. Since 2010 PISCO has estimated the size of lobsters encountered along transects. Lobster sized was visually estimated by trained observers. From 2003-2014 PISCO has also quantified urchin sizes by manually removing 150 purple urchins from each site (note: this is NOT at each transect, but only at the site level). Urchins were manually sized with calipers onboard. Our estimes of lobster and urchin density at the transect level will extend from 2003 to 2018, while our estimates of lobster size will cover 2010-2018 and urchin size will cover 2003-2014. 

# Note on units 
# We standarized all estimates of density as ind. per m2. Size for lobsters was carapace length (cm), and test diameter for urchins (cm). We estimated mass for lobsters and urchins based on published carapace/test size ~ weight relationships. All experimental units are prey consumed per m2 per predator per hour. 

#----------------------------------------------------------------------------------


#----------------------------------------------------------------------------------
## Get data
# Organize two data files one for size and one for abundance, each including both lobsters and urchins

library(here)
#pisco 

meta <- read.csv(here("data/PISCO", "Swath_Transects_2003-2018.csv")) %>%
  dplyr::select(site, MPA_STATUS, RESERVE, YEAR_MPA) %>%
  distinct() # this creates a metadata df with info on MPA status and time of MPA implementation. NOTE NOT ALL SMCA PROHIBIT TAKE OF LOBSTERS, ex. Anacapa west SMCA allows recreational and commercial lobster take


pus <- read.csv(here("data/PISCO", "MESFRA_STRPUR_size_frequency_2003-2014.csv")) %>% dplyr::select(-campus, -method, -observer, -depth) %>% filter(classcode != "MESFRAAD") %>% # this is the urchin size frequency data
  left_join(meta)

pf <- read.csv(here("data/PISCO", "PANINT_MESFRA_STRPUR_counts_from_swath_2003-2018.csv")) %>% filter(classcode != "MESFRAAD")%>% dplyr::select(-campus, -method, -observer, -heading, -depth) %>% # this is all the data for the swaths, which includes abundance of urchins and lobsters plus size data for lobsters after 2010
  left_join(meta) 

sz <- filter(pf, classcode == "PANINT", year >=2010) %>% # this is the size data for lobsters (and the frequency data for after 2010)
  select(year, month, day, site, side, zone, transect, classcode, size, count, MPA_STATUS, RESERVE, YEAR_MPA) %>%
  bind_rows(pus) %>% # combine with the urchin size data
  mutate(area = NA)

ab <- pf %>% group_by(year, month, day, site, side, zone, transect, classcode, MPA_STATUS, RESERVE, YEAR_MPA) %>%
  summarize(abundance = sum(count)) # calculate lobster (and urchin) abundance regardless of size


# SBC LTER

  # Put together size data
lob <- read.csv(here("data/LTER", "Lobster_Abundance_All_Years.csv"), header = T) %>% # get lobster abundance and size data from LTER
  na_if(-99999) %>% 
  mutate(DATE = as.character(DATE),
         SIZE_MM = SIZE_MM*0.1,
         MPA_STATUS = ifelse(SITE == "IVEE" | SITE == "NAPL", "SMCA", "REFERENCE"),
         RESERVE = ifelse(SITE == "IVEE" | SITE == "NAPL", "IN", "OUT"), 
         YEAR_MPA = ifelse(SITE == "IVEE" | SITE == "NAPL", 2012, 0)
         ) %>%
  separate(DATE, into = c("j1", "day", "j2"), sep = "[/]") %>%
  mutate(j1 = NULL, 
         j2 = NULL, 
         transect = NA, 
         NUM_AO = NULL, 
         classcode = "PANINT")

# NOTE: because of the differences in sampling styles, I am calling LTER transects "side" and LTER lobster replicates "zone"

names(lob) <- c("year", "month", "day", "site", "side-LTERtransect", "zone-LTERreplicate", "size", "count", "area", "mpa_status", "reserve", "year_mpa", "transect", "classcode")

lob <- lob[, c("year", "month", "day", "site", "side-LTERtransect", "zone-LTERreplicate", "transect", "classcode", "size", "count", "mpa_status", "reserve", "year_mpa", "area" )] # reorder the columns for match with PISCO


urc <- read.csv(here("data/LTER", "LTE_Urchin_All_Years_20190611.csv"), header = T) %>% # get urchin abundance and size data for SBC LTER LTE transects.
  filter(TREATMENT == "CONTROL", COMMON_NAME == "Purple Urchin") %>%
  select(YEAR, MONTH, DATE, SITE, TRANSECT, SIZE, COUNT) %>% 
  na_if(-99999) %>% 
  mutate(DATE = as.character(DATE),
         MPA_STATUS = ifelse(SITE == "IVEE" | SITE == "NAPL", "SMCA", "REFERENCE"),
         RESERVE = ifelse(SITE == "IVEE" | SITE == "NAPL", "IN", "OUT"), 
         YEAR_MPA = ifelse(SITE == "IVEE" | SITE == "NAPL", 2012, 0)) %>%
  separate(DATE, into = c("year", "month", "day"), sep = "[-]") %>%
  mutate(year = NULL, 
         month = NULL, 
         transect = NA, 
         NUM_AO = NULL, 
         classcode = "STRPURAD", 
         `zone-LTERreplicate` = NA, 
         area = NA)

names(urc) <- c("year", "month", "day", "site", "side-LTERtransect", "size", "count", "mpa_status", "reserve", "year_mpa", "transect", "classcode", "zone-LTERreplicate", "area")

urc <- urc[, c("year", "month", "day", "site", "side-LTERtransect", "zone-LTERreplicate", "transect", "classcode", "size", "count", "mpa_status", "reserve", "year_mpa", "area" )] # reorder the columns for match with PISCO

names(sz) <- c("year", "month", "day", "site", "side-LTERtransect", "zone-LTERreplicate", "transect", "classcode", "size", "count", "mpa_status", "reserve", "year_mpa", "area" )

urc$`side-LTERtransect` <- as.character(urc$`side-LTERtransect`)
lob$`side-LTERtransect` <- as.character(lob$`side-LTERtransect`)

size <- sz %>% # combine LTER size data with PISCO size data
  mutate(day = as.character(day)) %>%
  bind_rows(urc, lob) %>% 
  arrange(year, month, day, site, `side-LTERtransect`, `zone-LTERreplicate`, transect)

size <- size %>% # convert size df to a df with a row for each individual (rather than counts of individuals nested within sizes)
  group_by_at(vars(-count)) %>%
  summarize(count = sum(count)) %>%
  drop_na(count)%>%
  complete(count = full_seq(1:count, 1))%>%
  ungroup() %>%
  mutate(size = as.numeric(size))

write.csv(size, here("data/cleaned", "size_frequencies.csv"), row.names = F) # this is a cleaned data frame to examine patterns in size frequencies...

  # Organize an abundance df

ltlobab <- lob %>% group_by_at(vars(-size, -count)) %>%
  summarize(abundance = sum(count))

lturcab <- urc %>% group_by_at(vars(-size, -count)) %>%
  summarize(abundance = sum(count))



#----------------------------------------------------------------------------------

ks.test(size$size[size$classcode == "PANINT" & size$reserve == "IN"], 
        size$size[size$classcode == "PANINT" & size$reserve == "OUT"])

ks.test(size$size[size$classcode == "STRPURAD" & size$reserve == "IN"], 
        size$size[size$classcode == "STRPURAD" & size$reserve == "OUT"])

size %>%
  filter(year >= 2012) %>%
  group_by(reserve, classcode) %>%
  mutate(median = median(size, na.rm = T), 
            mean = mean(size, na.rm = T)) %>%
ggplot(aes(x = size, group = reserve))+
  geom_density(aes(fill = reserve), alpha = 0.75, adjust = 1.75)+
  geom_vline(aes(xintercept=median, col = reserve))+
  scale_fill_manual(values = c("#762a83", "#d73027"))+
  labs(x = "Size (cm)", y = "Density", fill = "")+
  facet_wrap(~classcode, scales = "free")+
  theme_pubclean()

size %>% group_by(reserve, classcode) %>%
  summarize(median = median(size, na.rm = T), 
            mean = mean(size, na.rm = T))


# Notes: 
#the lter data, supplemented by the data that we collected, is the only data that could be useful for modeling size shift in urchin populations due to lobsters. The pisco data is great for getting abundance, and size distributions for more sites... but really only has 4 years of overlapping data where both lobsters and urchins were sized. Furthermore, the pisco sampling takes sizes of urchins only at the anchor site but not along the transect.

# model shifts in median lobster size inside and outside of reserves

df <- size %>% filter(classcode == "PANINT") %>%
  group_by_at(vars(-size)) %>%
  summarize(size.median = median(size, na.rm = T))

lme <- lmer(size.median ~ reserve + (1|year) + (1|site), data = df)
summary(lme)
boxplot(size.median ~ reserve, data = df)

df <- size %>% filter(classcode == "STRPURAD") %>%
  group_by_at(vars(-size)) %>%
  summarize(size.median = median(size, na.rm = T))

lme <- lmer(size.median ~ reserve + (1|year) + (1|site), data = df)
summary(lme)
boxplot(size.median ~ reserve, data = df)


# some random plots to explore the data
ab %>% spread(classcode, abundance) %>% ggplot(aes(x = PANINT, y = STRPURAD))+
  geom_jitter()




ggplot(pf, aes(x = abundance))+
  geom_histogram()+
  facet_wrap(~classcode, scales = "free")














urc <- read.csv("data/LTER/LTE_Urchin_All_Years_20190611.csv", header = T) %>% 
  filter(TREATMENT == "CONTROL") %>% select(YEAR, MONTH, DATE, SITE, TRANSECT, SIZE, COUNT, COMMON_NAME) 

names(urc) <- tolower(names(urc))

urc %>% group_by(year, month, date, site, transect, common_name, size) %>%
  summarize(count = sum(count)) %>%
  group_by(year, month, date, site, transect, common_name, size) %>%
  drop_na(count)%>%
  complete(count = full_seq(1:count, 1))%>%
  # select(-count) %>%
  ungroup() %>%
  mutate(size = as.numeric(size)) %>%
  mutate(protection = ifelse(site == "IVEE" | site == "NAPL", "MPA", "FISHED")) %>%
  View()













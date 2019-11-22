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
source(here("code/", "setup.R"))
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
  summarize(abundance = sum(count)) %>%# calculate lobster (and urchin) abundance regardless of size
  mutate(area = 30*2)



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
  mutate(size = as.numeric(size) 
         # ,mass = ifelse(classcode == "PANINT", 0.001352821*(size*10)^2.913963113, 
         #               ifelse(classcode == "STRPURAD", 0.000592598*(size*10)^2.872636198*1.01, NA))
         )


write.csv(size, here("data/cleaned", "size_frequencies.csv"), row.names = F) # this is a cleaned data frame to examine patterns in size frequencies...

  # Organize an abundance df

ltlobab <- lob %>% group_by_at(vars(-size, -count)) %>%
  summarize(abundance = sum(count)) %>%
  group_by(year, site, classcode, `side-LTERtransect`) %>%
  summarize(abundance = sum(abundance), 
            area = 4*60*5,
            density = sum(abundance)/(4*60*5))


lturc <- read.csv(here("data/LTER", "Annual_All_Species_Biomass_at_transect.csv"), stringsAsFactors = F,na.strings ="-99999") %>% # this is the annual data
  dplyr::select("YEAR", "SITE", "TRANSECT", "SP_CODE", "DENSITY") %>%
  filter(SP_CODE == "SPL") %>%
  mutate(abundance = NA, 
         area = 40*2, 
         den = DENSITY, 
         DENSITY = NULL, 
         classcode = ifelse(SP_CODE == "SPL", "STRPURAD", NA), 
         SP_CODE = NULL)

names(lturc) <- tolower(names(lturc))
lturc <- lturc[, c("year", "site", "transect", "classcode", "den", "abundance", "area")]

ltesites <- c("AQUE", "CARP", "IVEE", "MOHK", "NAPL")

ltlob <- read.csv(here("data/LTER", "Annual_All_Species_Biomass_at_transect.csv"), stringsAsFactors = F,na.strings ="-99999") %>% # this is the annual data
  dplyr::select("YEAR", "SITE", "TRANSECT", "SP_CODE", "DENSITY") %>%
  filter(SP_CODE == "PAIN")

allpre2012 <- ltlob[ltlob$YEAR < 2012, ]
nonlteafter2012 <- ltlob[(!ltlob$SITE %in% ltesites) & ltlob$YEAR >=2012, ]

lobster <- rbind(allpre2012, nonlteafter2012) %>% mutate(abundance = NA, 
                                                         area = 40*2, 
                                                         classcode = "PANINT", 
                                                         SP_CODE =NULL)
names(lobster) <- c("year", "site", "transect", "den", "abundance", "area", "classcode")
lobster <- lobster[, c("year", "site", "transect","classcode", "den", "abundance","area")]


names(ltlobab) <- c("year", "site", "classcode", "transect", "abundance", "area", "den")

ltlobab <- ltlobab[, c("year", "site", "transect","classcode", "den", "abundance","area")]

lt <- ltlobab %>%
  mutate(transect = as.numeric(transect)) %>%
  bind_rows(lobster) %>%
  bind_rows(lturc) %>%
  arrange(year, site, classcode)%>%
  mutate(abundance = NULL, 
         area = NULL) %>%
  group_by(year, site, transect, classcode) %>%
  spread(classcode, den, fill = 0) %>%
  mutate(MPA_STATUS = ifelse(site == "IVEE" | site == "NAPL", "SMCA", "REFERENCE"),
         RESERVE = ifelse(site == "IVEE" | site == "NAPL", "IN", "OUT"), 
         YEAR_MPA = ifelse(site == "IVEE" | site == "NAPL", 2012, 0))
  
pisab <- ab %>%
  mutate(den = abundance/area, 
         abundance = NULL, 
         area = NULL) %>%
  group_by(year, site, side, zone, transect, classcode) %>%
  spread(classcode, den, fill = 0) %>%
  ungroup() %>%
  mutate(transectid = paste(side, zone, transect, sep = "-")) %>%
  dplyr::select(-month, -day, -side, -zone, -transect)

names(pisab) <- c("year", "site", "mpa_status", "reserve", "year_mpa", "PANINT", "STRPURAD", "transect")
pisab <- pisab[, c("year", "site", "transect", "PANINT", "STRPURAD","mpa_status", "reserve", "year_mpa")]

names(lt) <- c("year", "site", "transect", "PANINT", "STRPURAD","mpa_status", "reserve", "year_mpa")


den <- lt %>%
  ungroup() %>%
  mutate(transect = as.character(transect)) %>%
  bind_rows(pisab) %>%
  arrange(year, site, transect)

write.csv(abun, here("data/cleaned", "all_densities.csv"), row.names = F)

  

#----------------------------------------------------------------------------------

  # Analyze size data
ks.test(size$size[size$classcode == "PANINT" & size$reserve == "IN"], 
        size$size[size$classcode == "PANINT" & size$reserve == "OUT"])

ks.test(size$size[size$classcode == "STRPURAD" & size$reserve == "IN"], 
        size$size[size$classcode == "STRPURAD" & size$reserve == "OUT"])

size.plot <- size %>%
  filter(year >= 2011) %>%
  group_by(reserve, classcode) %>%
  # mutate(median = median(size, na.rm = T), 
  #           mean = mean(size, na.rm = T)) %>%
ggplot(aes(x = size, group = reserve))+
  geom_density(aes(fill = reserve), alpha = 0.75, adjust = 1.75)+
  # geom_vline(aes(xintercept=median, col = reserve))+
  scale_fill_manual(values = c("#7fbf7b", "#af8dc3"))+
  # scale_color_manual(values = c("#7fbf7b", "#af8dc3"))+
  labs(x = "Size (cm)", y = "Density", fill = "")+
  facet_wrap(~classcode+year, scales = "free")+
  theme_pubclean()

size %>%
  filter(year >= 2011) %>%
  group_by(reserve, classcode) %>%
  # mutate(median = median(size, na.rm = T), 
  #           mean = mean(size, na.rm = T)) %>%
  mutate(size = as.numeric(size) 
         ,mass = ifelse(classcode == "PANINT", 0.001352821*(size*10)^2.913963113,
                       ifelse(classcode == "STRPURAD", 0.000592598*(size*10)^2.872636198*1.01, NA))
  )%>%
  ggplot(aes(x = mass, group = reserve))+
  geom_density(aes(fill = reserve), alpha = 0.75)+
  # geom_vline(aes(xintercept=median, col = reserve))+
  scale_fill_manual(values = c("#7fbf7b", "#af8dc3"))+
  # scale_color_manual(values = c("#7fbf7b", "#af8dc3"))+
  labs(x = "Size (cm)", y = "Density", fill = "")+
  facet_wrap(~classcode+year, scales = "free")+
  theme_pubclean()

#this is only the LTER data
size %>%
  filter(year > 2011, classcode == "PANINT", site %in% c("NAPL", "IVEE", "CARP", "MOHK", "AQUE")) %>%
  group_by(reserve, classcode) %>%
  # mutate(median = median(size, na.rm = T), 
  #           mean = mean(size, na.rm = T)) %>%
  ggplot(aes(x = size, group = reserve))+
  geom_density(aes(fill = reserve), alpha = 0.75, adjust = 1)+
  coord_flip()+
  # geom_vline(aes(xintercept=median, col = reserve))+
  scale_fill_manual(values = c("#7fbf7b", "#af8dc3"))+
  # scale_color_manual(values = c("#7fbf7b", "#af8dc3"))+
  labs(x = "Size (cm)", y = "Density", fill = "")+
  facet_wrap(~year, nrow = 1)+
  geom_vline(xintercept = 8.25, linetype = "dashed")+
  theme_pubclean()


# so this is only the reserves that were implemented since 2012 and we can see significant shifts in the size distributions before and after implementation. Before implementation there were more legal sized lobsters outside of reserves. After implementation there were more legal sized lobsters inside reserves.
size %>%
  filter(classcode == "PANINT", year_mpa == 2012 | year_mpa == 0) %>%
  group_by(reserve, classcode) %>%
  # mutate(median = median(size, na.rm = T), 
  #           mean = mean(size, na.rm = T)) %>%
  ggplot(aes(x = size, group = reserve))+
  geom_density(aes(fill = reserve), alpha = 0.75, adjust = 1)+
  
  coord_flip()+
  # geom_vline(aes(xintercept=median, col = reserve))+
  scale_fill_manual(values = c("#7fbf7b", "#af8dc3"))+
  # scale_color_manual(values = c("#7fbf7b", "#af8dc3"))+
  labs(x = "Size (cm)", y = "Density", fill = "")+
  facet_wrap(~year, nrow = 1)+
  geom_vline(xintercept = 8.25, linetype = "dashed")+
  theme_pubclean()

# so these are only the older MPA's
size %>%
  filter(classcode == "PANINT", year_mpa < 2012) %>%
  group_by(reserve, classcode) %>%
  # mutate(median = median(size, na.rm = T), 
  #           mean = mean(size, na.rm = T)) %>%
  ggplot(aes(x = size, group = reserve))+
  geom_density(aes(fill = reserve), alpha = 0.75, adjust = 1)+
  
  coord_flip()+
  # geom_vline(aes(xintercept=median, col = reserve))+
  scale_fill_manual(values = c("#7fbf7b", "#af8dc3"))+
  # scale_color_manual(values = c("#7fbf7b", "#af8dc3"))+
  labs(x = "Size (cm)", y = "Density", fill = "")+
  facet_wrap(~year, nrow = 1)+
  geom_vline(xintercept = 8.25, linetype = "dashed")+
  theme_pubclean()

#this is everything
size %>%
  filter(classcode == "PANINT") %>%
  group_by(reserve, classcode) %>%
  # mutate(median = median(size, na.rm = T), 
  #           mean = mean(size, na.rm = T)) %>%
  ggplot(aes(x = size, group = reserve))+
  geom_density(aes(fill = reserve), alpha = 0.75, adjust = 1)+
  coord_flip()+
  # geom_vline(aes(xintercept=median, col = reserve))+
  scale_fill_manual(values = c("#7fbf7b", "#af8dc3"))+
  scale_x_log10()+
  # scale_color_manual(values = c("#7fbf7b", "#af8dc3"))+
  labs(x = "Size (cm)", y = "Density", fill = "")+
  facet_wrap(~year, nrow = 1)+
  geom_vline(xintercept = 8.25, linetype = "dashed")+
  theme_pubclean()

size %>%
  filter(classcode == "PANINT") %>%
  group_by(reserve, classcode) %>%
  # mutate(median = median(size, na.rm = T), 
  #           mean = mean(size, na.rm = T)) %>%
  ggplot(aes(x = size))+
  geom_density(alpha = 0.75, adjust = 1)+
  coord_flip()+
  # geom_vline(aes(xintercept=median, col = reserve))+
  scale_fill_manual(values = c("#7fbf7b", "#af8dc3"))+
  scale_x_log10()+
  # scale_color_manual(values = c("#7fbf7b", "#af8dc3"))+
  labs(x = "Size (cm)", y = "Density", fill = "")+
  facet_wrap(~year+reserve, nrow = 1)+
  geom_vline(xintercept = 8.25, linetype = "dashed")+
  theme_pubclean()





ggsave(here("figures", "size-plot.png"), size.plot, device = "png", width = 8.5, height = 8.5*0.6 )

size %>% group_by(reserve, classcode) %>%
  summarize(median = median(size, na.rm = T), 
            mean = mean(size, na.rm = T))


# Notes: 
#the lter data, supplemented by the data that we collected, is the only data that could be useful for modeling size shift in urchin populations due to lobsters. The pisco data is great for getting abundance, and size distributions for more sites... but really only has 4 years of overlapping data where both lobsters and urchins were sized. Furthermore, the pisco sampling takes sizes of urchins only at the anchor site but not along the transect.

# model shifts in median lobster size inside and outside of reserves

df <- size %>% filter(classcode == "PANINT") #%>%
#   group_by_at(vars(-size)) %>%
#   summarize(size.median = median(size, na.rm = T))

lme <- lmer(size ~ reserve + (1|year) + (1|site), data = df[df$year >= 2012, ])
summary(lme)
boxplot(size ~ reserve, data = df)
coefs <- data.frame(reserve = c("IN", "OUT"), y = c(8.1772, 8.1772-1.5775), high = c(8.1772 + 0.2946, (8.1772-1.5775 + 0.3954)), low = c(8.1772 - 0.2946, (8.1772-1.5775 - 0.3954)))

lme2 <- lmer(size ~ reserve *year + (1|site), data = df[df$year >= 2011, ])
summary(lme2)

newdat <- expand.grid(reserve = c("IN", "OUT"), year = seq(2012, 2019, length.out = 1000))
newdat$predicted <- predict(lme2, newdata = newdat, type = "response", re.form = NA)


ggplot(df[df$year >= 2012, ], aes(x = year, y = size))+
  geom_jitter(aes(color = reserve))+
  scale_color_manual(values = c("#7fbf7b", "#af8dc3"))+
  geom_line(data = newdat, aes(x = year, y = predicted, group = reserve))+
  theme_pubclean()


lob.size <- ggplot(df[df$year >= 2012, ], aes(x = reserve, y = size))+
  geom_jitter(aes(color = reserve), show.legend = FALSE)+
  scale_color_manual(values = c("#7fbf7b", "#af8dc3"))+
  geom_point(data = coefs, aes(x = reserve, y = y))+
  geom_linerange(data = coefs, aes(x = reserve, y = y, ymin = low, ymax = high))+
  labs(x = NA, y = "Lobster carapace length (cm)") +
  annotate("text", label = "***", x = 1, y = 12.5, cex = 8)+
  theme_pubclean()

ggsave(here("figures", "lob-size.png"), lob.size, device = "png", width = 4, height = 4*1.3, bg = "transparent" )


df <- size %>% filter(classcode == "STRPURAD")# %>%
  # group_by_at(vars(-size)) %>%
  # summarize(size.median = median(size, na.rm = T))

lme <- lmer(size ~ reserve + (1|year) + (1|site), data = df)
summary(lme)
boxplot(size.median ~ reserve, data = df)

#----------------------------------------------------------------------------------
## Analyze density data for patterns inside and outside of MPA's, and for correlations between lobster abundance and urchin abundance

# some random plots to explore the data
den %>%
  ggplot(aes(x = log10(PANINT+1), y = log10(STRPURAD+1)))+
  geom_jitter(aes(color = reserve))

den %>%
  ggplot(aes(x = PANINT, y = STRPURAD))+
  geom_jitter(aes(color = reserve), pch = 21)


den %>% group_by(reserve) %>%
  summarize(lobster = mean(PANINT, na.rm = T), 
            urchin = mean(STRPURAD, na.rm = T))

ggplot(den, aes(x = PANINT))+
  geom_density(aes(fill = reserve))


ggplot(den, aes(x = PANINT, y = STRPURAD))+
  geom_jitter(aes(color = reserve), pch = 21)+
  geom_line(data = data.frame(x = exp(lme4@frame$`log(PANINT + 1)`)-1,
                              y = exp(predict(lme4, re.form = NA, type = "response"))-1), 
            aes(x = x, y = y))


# urchin ~ lobster densities plot
lme4 <- lmer(log(STRPURAD+1) ~ log(PANINT+1) + (1|year) + (1|site), den)
summary(lme4)
newdat <- data.frame(PANINT = seq(0, log(max(den$PANINT)+1), length.out = 1000))
newdat$pred <- predict(lme4, newdata = newdat, re.form = NA)

lob.urc <- ggplot(den, aes(x = log(PANINT+1), y = log(STRPURAD+1)))+
  geom_jitter(aes(color = reserve))+
  scale_color_manual(values = c("#7fbf7b", "#af8dc3"))+
  geom_line(data = newdat, aes(y = pred, x = PANINT))+
  labs(y = "Purple urchin density (ind. m2)", x = "Lobster density (ind. m2)")

ggsave(here("figures", "lobster-urchin.png"), lob.urc, device = "png", width = 7, height = 5.63)

# lobster density ~ reserve status plot 

lme3 <- lmer(PANINT ~ reserve + (1|year) + (1|site), den[den$year >= 2012, ])
summary(lme3)
inn <- 0.022053
out <- 0.013674

coefs <- data.frame(reserve = c("IN", "OUT"), y = c(inn, inn-out), high = c(inn + 0.004132, (inn-out + 0.003874)), low = c(inn - 0.004132, (inn-out - 0.003874)))

lob.reserve <- ggplot(den[den$year >= 2012, ], aes(x = reserve, y = PANINT))+
  geom_jitter(aes(color = reserve), show.legend = FALSE)+
  scale_color_manual(values = c("#7fbf7b", "#af8dc3"))+
  geom_point(data = coefs, aes(x = reserve, y = y))+
  geom_linerange(data = coefs, aes(x = reserve, y = y, ymin = low, ymax = high))+
  labs(x = NA, y = "Lobster density (ind. m2)") +
  annotate("text", label = "~61% increase in \ndensity inside reserves", x = 1.5, y = 0.5, cex = 4)+
  annotate("text", label = "***", x = 1, y = 0.3, cex = 8) 

ggsave(here("figures", "lobster-reserve.png"), lob.reserve, device = "png", width = 3, height = 5.63)

lme5 <- lmer(STRPURAD ~ reserve + (1|year) + (1|site), den)
summary(lme5)
boxplot(STRPURAD ~ reserve, den)


#------------------------------------------------------------------------------------
## Bodysize ratio histograms
#------------------------------------------------------------------------------------


size %>%
  filter(year >=2010) %>%
  group_by(reserve, classcode) %>%
  ggplot(aes(x = size, group = classcode))+
  geom_density(aes(fill = classcode), alpha = 0.75, adjust = 1.75)+
  # geom_vline(aes(xintercept=median, col = reserve))+
  scale_fill_manual(values = c("#7fbf7b", "#af8dc3"))+
  # scale_color_manual(values = c("#7fbf7b", "#af8dc3"))+
  coord_flip()+
  labs(x = "Size (cm)", y = "Density", fill = "")+
  facet_wrap(~year, nrow = 1)+
  theme_pubclean()

size %>%
  filter(year >=2010) %>%
  group_by(reserve, classcode, site) %>%
  ggplot(aes(x = size))+
  geom_density(aes(fill = site), alpha = 0.75, adjust = 1.75, show.legend = F)+
  # geom_vline(aes(xintercept=median, col = reserve))+
  #scale_fill_manual(values = c("#7fbf7b", "#af8dc3"))+
  # scale_color_manual(values = c("#7fbf7b", "#af8dc3"))+
  coord_flip()+
  labs(x = "Size (cm)", y = "Density", fill = "")+
  facet_wrap(year ~ classcode, nrow = 1)+
  theme_pubclean()

size %>%
  filter(year >=2012, site %in% c("NAPL", "IVEE", "CARP", "MOHK", "AQUE")) %>%
  group_by(reserve, classcode, site) %>%
  ggplot(aes(x = size))+
  geom_density(aes(fill = classcode, group = site), alpha = 0.5, adjust = 1.75, show.legend = F)+
  scale_fill_manual(values = c("#7fbf7b", "#af8dc3"))+
  coord_flip()+
  labs(x = "Size (cm)", y = "Density", fill = "")+
  facet_wrap(year ~ classcode, nrow = 1)+
  theme_pubclean()


size %>%
  filter(year >=2012, site %in% c("NAPL", "IVEE", "CARP", "MOHK", "AQUE")) %>%
  group_by(reserve, classcode, site) %>%
  ggplot(aes(x = size))+
  geom_density(aes(fill = site), alpha = 0.5, adjust = 1.75, show.legend = F)+
  scale_fill_manual(values = c('#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3'))+
  labs(x = "Size (cm)", y = "Density", fill = "")+
  facet_wrap( ~ classcode, scales = "free", ncol = 1)+
  theme_pubclean()





size %>%
  filter(year >=2012, site %in% c("NAPL", "IVEE", "CARP", "MOHK", "AQUE")) %>%
  group_by(reserve, classcode) %>%
  ggplot(aes(x = size, group = classcode))+
  geom_density(alpha = 0.5, adjust = 1.75, show.legend = F)+
  # geom_vline(aes(xintercept=median, col = reserve))+
  scale_fill_manual(values = c("#7fbf7b", "#af8dc3"))+
  scale_color_manual(values = terrain.colors(5))+
  coord_flip()+
  labs(x = "Size (cm)", y = "Density", fill = "")+
  facet_wrap(site + year ~ classcode, nrow = 5)+
  theme_pubclean()


size %>%
  group_by(year, classcode, site) %>%
  filter(year >=2010) %>%
  summarize(median = median(size, na.rm = T)) %>%
  spread(classcode, median) %>%
  mutate(ratio = PANINT/STRPURAD) %>%
  ggplot(aes(x = ratio))+
    geom_histogram()


#----------------------------------------------------------

size.temp <- size %>%
  mutate(size = as.numeric(size) 
         ,size = ifelse(classcode == "PANINT", 0.001352821*(size*10)^2.913963113,
                        ifelse(classcode == "STRPURAD", 0.000592598*(size*10)^2.872636198*1.01, NA)))


lob.list <- size.temp %>% filter(year >= 2012, classcode == "PANINT", site %in% c("NAPL", "IVEE", "CARP", "MOHK", "AQUE") ) %>% 
  dplyr::select(size, site, year) %>%
  mutate(id = paste(site, year, sep = "")) %>%
  mutate(id = as.factor(id)) %>%
  as.data.frame()

lob.list <- split(lob.list$size, lob.list$id)

urc.list <- size.temp %>% filter(year >= 2012, classcode == "STRPURAD", site %in% c("NAPL", "IVEE", "CARP", "MOHK", "AQUE") ) %>% 
  dplyr::select(size, site, year) %>%
  mutate(id = paste(site, year, sep = "")) %>%
  mutate(id = as.factor(id)) %>%
  as.data.frame()

urc.list <- split(urc.list$size, urc.list$id)


numid <- length(lob.list)
out <- list()
for(i in 1:numid){
  temp.lob <- lob.list[[i]]
  temp.urc <- urc.list[[i]]
  out[[i]] <- matrix(nrow = length(temp.urc), ncol = length(temp.lob))
  for(k in 1:length(temp.lob)){
    out[[i]][,k] <- temp.lob[k] / temp.urc
  }
}

names(out) <- names(lob.list)

temp.list <- list()
for(i in 1:numid){
  temp.list[[i]] <- data.frame(ratio = as.numeric(na.omit(unlist(out[i]))), id = names(out[i]))
}

ratio <- do.call(rbind, temp.list) %>%
  separate(id, into = c("site", "year"), sep = "(?<=[A-Za-z])(?=[0-9])")

ggplot(ratio, aes(x = ratio))+
  geom_density(aes(fill = site), alpha = 0.5, adjust = 1.75)+
  scale_x_log10()


fig1 <- ggplot(ratio, aes(x = ratio))+
  geom_density(aes(fill = site), alpha = 0.5, adjust = 1.75, show.legend = T)+
  scale_fill_manual(values = c('#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3'))+
  scale_x_log10()+
  coord_flip()+
  labs(x = "Predator:prey ratio", y = "Density", fill = "")+
  facet_wrap(~year, nrow = 1)+
  theme_pubclean()+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggsave(here("figures", "ratio-histogram.png"), fig1, device = "png", width = 8.5, height = 6)

ggplot(ratio, aes(x = ratio))+
  geom_density(aes(group = site), alpha = 0.5, adjust = 1.75, show.legend = T)+
  scale_fill_manual(values = c('#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3'))+
  scale_x_log10()+
  labs(x = "Predator:prey ratio", y = "Density", fill = "")+
  facet_wrap(year~site, nrow = 8, ncol = 5)+
  theme_pubclean()+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


ggplot(ratio, aes(x = ratio, y = year))+
  geom_density_ridges(aes(fill = site), rel_min_height = 0.01, alpha = 0.8, scale = 4)+
  scale_fill_manual(values = c('#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3'))+
  scale_x_log10()+
  facet_wrap(~site)+
  theme_pubclean()

ggplot(lincoln_weather, aes(x = `Mean Temperature [F]`, y = `Month`, fill = ..x..)) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
  scale_fill_viridis(name = "Temp. [F]", option = "C") +
  labs(title = 'Temperatures in Lincoln NE in 2016')



library(dplyr)
ratio %>% group_by(site, year) %>%
  do(ggplot2:::compute_density(log10(.$ratio), NULL, adjust = 1.75)) %>%
  rename(ratio = x) %>%
  mutate(ratio = 10^ratio) -> ratio_densities
head(ratio_densities)


fig1 <- ggplot(ratio_densities, aes(x = ratio, y = year, height = density))+
  geom_density_ridges(aes(fill = site), stat = "identity", scale = 2, alpha = 0.9)+
  scale_fill_manual(values = c('#8dd3c7','#ffffb3','#bebada','#fb8072','#80b1d3'))+
  scale_x_log10(breaks = c(0, 1, 10, 100), labels = c(0, 1,10, 100))+
  coord_cartesian(xlim = c(0.4, 500))+
  labs(y = "", x = "Potential predator:prey body size ratios")+
  facet_wrap(~site, nrow = 1)+
  theme_pubclean()+
  theme(strip.background = element_blank(),
        strip.text.x = element_blank())

ggsave(here("figures", "ratio-ridgeplot.png"), fig1, device = "png", width = 8.5, height = 8.5/2.5)

# ggplot(ratio_densities, aes(x = ratio, y = year, height = density))+
#   geom_density_ridges_gradient(aes(fill = log(ratio)), stat = "identity", scale = 2, alpha = 0.8, rel_min_height = 0.001, show.legend = F)+
#   scale_fill_gradient(low = "#40004b", high = "#00441b")+
#   scale_x_log10(breaks = c(0, 1, 10, 100), labels = c(0, 1,10, 100), lim = c(0.001, 500))+
#   labs(y = "", x = "Potential predator:prey body size ratios")+
#   facet_wrap(~site, nrow = 1)+
#   theme_pubclean()
#   # theme(strip.background = element_blank(),
#   #       strip.text.x = element_blank())




ggplot(lincoln_weather, aes(x = `Mean Temperature [F]`, y = `Month`, fill = ..x..)) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
  scale_fill_viridis(name = "Temp. [F]", option = "C") +
  labs(title = 'Temperatures in Lincoln NE in 2016')


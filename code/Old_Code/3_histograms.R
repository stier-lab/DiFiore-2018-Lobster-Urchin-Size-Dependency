####################################################################################
## Size Frequency Distributions
####################################################################################
library(here)
source(here("code", "1_setup.R"))

# Get data

lob <- read.csv("data/LTER/Lobster_Abundance_All_Years.csv", header = T) %>% 
  na_if(-99999) %>% select(YEAR, MONTH, DATE, SITE, SBC_LTER_TRANSECT, SIZE, COUNT) %>%
  add_column(COMMON_NAME = "Spiny Lobster") %>%
  rename(TRANSECT = SBC_LTER_TRANSECT) %>%
  mutate(PROTECTION = ifelse(SITE == "IVEE" | SITE == "NAPL", "MPA", "FISHED"), SIZE = SIZE*0.1) %>%
  mutate(bins = cut_width(SIZE, width = .5))


urc <- read.csv("data/LTE_Urchin_All_Years_20190611.csv", header = T) %>% 
  filter(TREATMENT == "CONTROL", COMMON_NAME == "Purple Urchin") %>% select(YEAR, MONTH, DATE, SITE, TRANSECT, SIZE, COUNT, COMMON_NAME) %>% mutate(PROTECTION = ifelse(SITE == "IVEE" | SITE == "NAPL", "MPA", "FISHED"))%>% mutate(bins = cut_width(SIZE, width = .5))

names(urc) <- tolower(names(urc))

urc <- urc %>% group_by(year, month, date, site, transect, common_name, size) %>%
  summarize(count = sum(count)) %>%
  spread(size, count) %>%
  gather(size, count, -c(year, month, date, site, transect, common_name)) %>%
  group_by(year, month, date, site, transect, common_name, size) %>%
  summarize(count = sum(count)) %>%
  group_by(year, month, date, site, transect, common_name, size) %>%
  drop_na(count)%>%
  complete(count = full_seq(1:count, 1))%>%
  select(-count) %>%
  ungroup() %>%
  mutate(size = as.numeric(size))
  

quants <- quantile(urc$size, prob = c(0.25, 0.975))

bysite <- ggplot(urc, aes(x = size))+
  geom_density(aes(fill = site), alpha = 0.5, adjust = 2)

ggsave("figures/urchinsize_bysite.png", bysite, device = "png")

hist <- ggplot(urc, aes(x = size))+
  geom_histogram(stat = "count")+
  geom_vline(xintercept = quantile(urc$size, prob = c(0.25, 0.975)), color = "red", linetype = "dashed")

ggsave("figures/urchinsizehist.png", hist, device = "png")


p1 <- lob %>%
  group_by(bins, PROTECTION) %>% 
  summarize(freq = sum(COUNT)) %>%
  ggplot(aes(x=bins, y=freq))+
  geom_bar(aes(fill = PROTECTION), stat = "identity", position = "dodge") +
  scale_x_discrete(labels = seq(min(lob$SIZE, na.rm = T), max(lob$SIZE, na.rm = T), by = .5))+
  labs(x = "Carapace length (cm)", y = "Frequency", title = "Lobster size distribution") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


p2 <- urc %>%
  group_by(bins, PROTECTION) %>% 
  summarize(freq = sum(COUNT)) %>%
  ggplot(aes(x=bins, y=freq))+
  geom_bar(aes(fill = PROTECTION), stat = "identity", position = "dodge") +
  scale_x_discrete(labels = seq(min(urc$SIZE, na.rm = T), max(urc$SIZE, na.rm = T), by = .5))+
  labs(x = "Test diameter (cm)", y = "Frequency", title = "Urchin size distribution")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

png("figures/fig_hist.png", width = 1000, height = 350)
cowplot::plot_grid(p1+guides(fill=FALSE), p2+labs(y = "")+theme(legend.position = c(0.7, 0.5)))
dev.off()


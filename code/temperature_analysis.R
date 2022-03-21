library(lubridate)

df <- read.csv("data/raw/temperature_log.csv") %>%
  mutate(dmy = lubridate::dmy(dmy)) %>%
  separate(dmy, into = c("year", "month", "day"), sep = "[-]") %>%
  mutate(year = as.numeric(year), 
         month = as.numeric(month), 
         day = as.numeric(day)) %>%
  filter(year == 2018, 
         month >=6 & month <= 9)

df %>%
  summarize(mean_temp = mean(temp), 
            sd_temp = sd(temp))

summary(df$temp)

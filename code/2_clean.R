##############################################################################
## Clean raw csv file to prepare for fitting functional responses
##############################################################################
library(here)
source(here("code","1_setup.R"))
source(here("code","Base_functions/functions.R"))

df <- read.table(here("data/raw","loburc_raw.csv"), header = T, sep = ",")

df <- df %>% 
  filter(finished == "finished", discard != "yes") %>%
  drop_na(num_offered) %>%
  select(-last_molt, -sex, -finished, -starve_hours, -food_remaining, -discard, -notes)


###############################################################################
## Table to visualize replication
###############################################################################

table(df$id, df$num_offered) # So there was obviously an issue and some lobsters did not recieve all of the densities. Not sure of the implications of this...


table(df$class, df$num_offered, df$treatment) # We have a problem with replication in the jumbo size class which we knew was going to happen.

################################################################################
## Clean up the data and add in sizes
################################################################################

df <- df %>%
  filter(trial != "special")

df <- df[,c("id", "size", "treatment", "num_offered", "num_left", "num_dead")]

df[is.na(df)] <- 0

df$killed <- ifelse((df$num_offered - df$num_dead) - df$num_left < 0, 0, (df$num_offered - df$num_dead) - df$num_left)
df$initial <- df$num_offered - df$num_dead
df <- df[,c("id", "size", "treatment", "initial", "killed")]

df$udiam <- ifelse(df$treatment == "urc_medium", 40,
                   ifelse(df$treatment == "urc_large", 60, 20))


################################################################

# Get weights (g) for urchins

df$mr <- 0.000592598*(df$udiam^2.872636198)*1.01

# get lobster weights

lw <- read.csv(here("data/raw", "lobster_weights.csv"))

names(lw) <- c("id", "w1", "w2", "w3")

lw <- lw %>% group_by(id) %>%
  gather(testnum, weight, -id) %>%
  summarize(mc = mean(weight))

df <- left_join(df, lw) %>% mutate(id = as.factor(id))




###############################################################################
## Write csv
###############################################################################

write.table(df, here("data/cleaned","loburc_cleaned.csv"), sep = ",", quote = F, row.names = F)

#################################################################
## Lobster metadata file
#################################################################

md <- distinct(df, id, size, treatment, udiam, mr, mc)

write.csv(md, here("data", "lob-metadata.csv"), row.names = F)


summary(md$mc)









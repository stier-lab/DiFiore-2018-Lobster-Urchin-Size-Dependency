##############################################################################
## Clean raw csv file to prepare for fitting functional responses
##############################################################################

source(here("code","setup.R"))
source(here("code","functions.R"))

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


###############################################################################
## Write csv
###############################################################################

write.table(df, here("data/cleaned","loburc_cleaned.csv"), sep = ",", quote = F, row.names = F)

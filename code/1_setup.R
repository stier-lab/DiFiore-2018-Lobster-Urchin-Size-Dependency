##########################################
## Libraries
##########################################

library(tidybayes)
library(ggplot2)
library(tidyr)
library(dplyr)
library(cowplot)
library(here)
library(lme4)
library(lmerTest)
library(purrr)
library(rethinking)
library(rstanarm)

colors <- c('#AF8DC3','#C3AF8D','#8DC3AF')

# size of experimental tanks (sides and bottom -- i.e. areal NOT volumetric)
tsize <- 137/2 * 76 /10000



set.seed(112618)

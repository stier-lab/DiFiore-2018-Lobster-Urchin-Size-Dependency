#----------------------------------------------
## Figures for committee meeting
#----------------------------------------------

library(here)
source(here("code", "setup.R"))
source(here("code", "functions.R"))
library(frair)

#####################################
## Get data
#####################################


df <- read.table(here("data/cleaned","loburc_cleaned.csv"), header = T, sep = ",") %>% arrange(treatment, id)

# plot global fit

global <- mle2(hollingsII_nll,
             start=list(a=0.02,h=1),
             data=list(Y = df$killed, X = df$initial, T=48),
             method = "Nelder-Mead")



png(here("figures", "globalfit.png"), width = 1000*2, height = 563*2, res = 300)
par(mar = c(4,5,1,1))
plot(I(killed/48) ~ jitter(initial), df, ylab = "Consumption rate \n(ind. per predator per hour)", xlab = "Initial number of prey", ylim = c(0,1))
curve(holling2(x,coef(global)[1],coef(global)[2],P=1,T=1),add=TRUE,lty=1, lwd = 2, col = "darkblue") 
dev.off()

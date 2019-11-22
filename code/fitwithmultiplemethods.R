
########################################################
## read in the lobster data 
#########################################################

df <- read.table(here("data/cleaned","loburc_cleaned.csv"), header = T, sep = ",")

x <- data.frame(killed = df$killed, initial = df$initial, ind = as.numeric(df$id))

df <- fit.individual(x, T = 48,  plot = F) 
fit.individual(x , plot = T, T = 48)



# ------------------------------------------------------------
# examine variation in estimates between different packages


# bring in JAGS data
jg <- read.csv(here("data", "JAGSparam-estimates.csv"))

forbind <- jg[,c("ind", "package", "mean.a", "mean.h")]
names(forbind) <- c("ind", "package", "a", "h")

df <- df %>% bind_rows(forbind) %>% group_by(ind) %>% arrange(ind, package)
write.csv(df, here("data", "paramest-all.csv"), row.names = F)

p1 <- ggplot(df)+
  coord_flip() + 
  geom_hline(yintercept = 0, colour = gray(1/2), lty = 2)+
  # geom_linerange(aes(x = vnames, ymin = Coefficient - SE, ymax = Coefficient + SE, color = Model), lwd = 1, position = position_dodge(width = 1/2))+
  geom_point(aes(x = as.factor(ind), y = a, color = package), lwd = 1, position = position_dodge(width = 0.5), shape = 21)+
  scale_y_continuous(limits = c(0,0.25) )+
  theme_bw()+
  theme(axis.title.y = element_blank(), axis.text = element_text(size = 12), axis.title = element_text(size = 14))

p2 <- ggplot(df)+
  coord_flip() + 
  geom_hline(yintercept = 0, colour = gray(1/2), lty = 2)+
  # geom_linerange(aes(x = vnames, ymin = Coefficient - SE, ymax = Coefficient + SE, color = Model), lwd = 1, position = position_dodge(width = 1/2))+
  geom_point(aes(x = as.factor(ind), y = h, color = package), lwd = 1, position = position_dodge(width = 0.5), shape = 21)+
  theme_bw()+
  theme(axis.title.y = element_blank(), axis.text = element_text(size = 12), axis.title = element_text(size = 14))

cowplot::plot_grid(p1, p2, nrow = 1)


fp <- df %>% group_by(ind) %>%
  gather(parameter, estimate, -c(ind, package)) %>%
  mutate(paste = paste(package, parameter, sep = ""), 
         parameter = NULL, 
         package = NULL) %>%
  spread(paste, estimate)

png(here("figures", "paramter-compare.png"), width = 400*3, height = 1200*3, res =150)
d <- par(mfcol = c(6,2), mar = c(4,4,1,0.5))
plot(mle2a ~ rogersa, fp, ylim = c(0,0.1), xlim = c(0, 0.1))
abline(a=0,b=1)
plot(mle2a ~ nlsa, fp, ylim = c(0,0.1), xlim = c(0, 0.1))
abline(a=0,b=1)
plot(rogersa ~ nlsa, fp)
abline(a=0,b=1)
plot(JAGSa ~ mle2a, fp, xlim = c(0, 0.1), ylim = c(0,1))
abline(a=0,b=1)
plot(JAGSa ~ rogersa, fp, xlim = c(0, 0.1), ylim = c(0,1))
abline(a=0,b=1)
plot(JAGSa ~ nlsa, fp, xlim = c(0, 0.1), ylim = c(0,1))
abline(a=0,b=1)

plot(mle2h ~ rogersh, fp)
abline(a=0,b=1)
plot(mle2h ~ nlsh, fp)
abline(a=0,b=1)
plot(rogersh ~ nlsh, fp)
abline(a=0,b=1)
plot(JAGSh ~ mle2h, fp)
abline(a=0,b=1)
plot(JAGSh ~ rogersh, fp)
abline(a=0,b=1)
plot(JAGSh ~ nlsh, fp)
abline(a=0,b=1)

par(d)

dev.off()

psych::pairs.panels(fp[,c("JAGSa", "mle2a", "rogersa", "nlsa")], ellipses = F)
psych::pairs.panels(fp[,c("JAGSh", "mle2h", "rogersh", "nlsh")], ellipses = F)

# ------------------------------------------------------------










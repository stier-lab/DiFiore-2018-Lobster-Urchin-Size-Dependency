library(here)
source(here("code", "setup.R"))
source(here("code", "functions.R"))


#####################################
## Get data
#####################################

df <- read.table(here("data/cleaned","loburc_cleaned.csv"), header = T, sep = ",") %>% 
  mutate(lobcat = ifelse(size <=70, "small", 
                         ifelse(size >70 & size <=90, "medium", 
                                ifelse(size > 90, "large", NA))), 
         newtreat = as.factor(paste(treatment, lobcat, sep = "-")),
         m.int = initial*mr, 
         m.kill = killed*mr,
         prop.m.kill = m.kill/m.int)  %>%
  arrange(newtreat, id)



##################################################################
## the model
##################################################################

model.loc=here("code","heirarchical_jags.txt")
jags.params=c("a", "h", "t.a", "t.h", "mu.a", "mu.h", 
              "sigma_int.a", "sigma_int.h", "sigma_t.a", "sigma_t.h", 
              "loga", "logh", "t.logit.a", "t.log.h"
)


tind <- distinct(df, newtreat, id) %>%
  arrange(id)

tind <- as.factor(as.vector(tind$newtreat))

jags.data = list("initial"= round(df$m.int, 0),
                 "killed" = round(df$m.kill, 0),
                 "id" = df$id,
                 "num.ind" = length(unique(df$id)),
                 "n" = length(df$initial), 
                 "tind" = tind, 
                 "T" = 48, 
                 "Ntreats" = length(unique(df$newtreat))
) # named list

n.chains = 3
n.burnin = 10000
n.thin = 2
n.iter = 25000
model = R2jags::jags(jags.data,parameters.to.save=jags.params,inits=NULL,
                     model.file=model.loc, n.chains = n.chains, n.burnin=n.burnin,n.thin=n.thin, n.iter=n.iter, DIC=TRUE)

a = MCMCsummary(model,params='a', round = 4)
h = MCMCsummary(model,params='h', round = 4)

df.ind <- data.frame(ind = 1:46, package = "JAGS", median.a = a[,4],
                     a.low = a[,3],
                     a.high = a[,5],
                     median.h = h[,4],
                     h.low = h[,3],
                     h.high = h[,5])
# write.csv(df.ind, here("data", "JAGSparam-estimates.csv"), row.names = F)

t.a = MCMCsummary(model,params='t.a', round = 4)
t.h = MCMCsummary(model,params='t.h', round = 4)

mu.a = MCMCsummary(model,params='mu.a', round = 4)
mu.h = MCMCsummary(model,params='mu.h', round = 4)



#----------------------------------------------------------------------------------
# for poster
#----------------------------------------------------------------------------------


newdat <- expand.grid(initial = seq(min(df$m.int, na.rm = T), max(df$m.int, na.rm = T), length.out = 100), newtreat = unique(df$newtreat)) %>% arrange(newtreat) %>%
  separate(newtreat, into = c("treatment", "lobcat"), sep = "[-]")

df %>% group_by(treatment) %>%
  summarize(min = min(m.int), 
            max = max(m.int) 
            )
initial <- c(seq(min(df$m.int[df$treatment == "urc_large"], na.rm =T), max(df$m.int[df$treatment == 'urc_large'], na.rm = T), length.out = 100), seq(min(df$m.int[df$treatment == 'urc_medium'], na.rm =T), max(df$m.int[df$treatment == 'urc_medium'], na.rm = T), length.out = 100), seq(min(df$m.int[df$treatment == 'urc_small'], na.rm =T), max(df$m.int[df$treatment == 'urc_small'], na.rm = T), length.out = 100))

newdat <- data.frame(treatment = rep(unique(df$treatment), each = 100), initial = initial)

newdat <- rbind(newdat, newdat, newdat)
newdat$lobcat <- rep(unique(df$lobcat), each = 300)
newdat <- arrange(newdat, treatment, lobcat)
newdat$a <- rep(t.a$`50%`, each = 100)
newdat$h <- rep(t.h$`50%`, each = 100)
newdat$killed <- holling2(N = newdat$initial, a = newdat$a, h = newdat$h, P = 1, T = 48)



jags.samples(model, c(t.a, t.h), thin = 10)
temp <- MCMCpstr(model, params = c("t.a", "t.h"), type = "chains")
temp <- MCMCchains(model, params = c("t.a"))

ggplot(df, aes(x = m.int, y = m.kill))+
  geom_jitter(aes(size = treatment, fill = lobcat), pch = 21, show.legend = F, stroke = 1)+
  scale_size_manual(values = c(4,2.5,1.5))+
  scale_fill_manual(values = c('#d53e4f','#fc8d59','#fee08b'))+
  geom_line(data = newdat, aes(x = initial, y = killed, color = lobcat), size = 1.5, show.legend = F)+
  scale_color_manual(values = c('#d53e4f','#fc8d59','#fee08b'))+
  facet_wrap(~treatment, scales = "free_x")

#ggsave(here("figures", "forposter1.svg"), poster1, device = "svg", width = 10, height = 10/3)

# --------------------------------------------------------------------------

df.ind %>% ggplot()



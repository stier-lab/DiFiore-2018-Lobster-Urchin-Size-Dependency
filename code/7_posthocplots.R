#-------------------------------------------
## Setup
#-------------------------------------------

library(here)
source(here("code", "Base_functions/Functions.R"))
source(here("code", "1_setup.R"))

#-------------------------------------------
## Get data
#-------------------------------------------

post.a <- read.csv(here::here("data/cleaned/posteriors", "posteriors_posthoc_a.csv")) %>%
  spread(parameter,estimate) %>% filter(model == "model_wpuncert")
post.h <- read.csv(here::here("data/cleaned/posteriors", "posteriors_posthoc_h.csv")) %>%
  spread(parameter,estimate) %>% filter(model == "model_wpuncert")

post.null <- read.csv(here::here("data/cleaned/posteriors", "posteriors_null.csv")) %>%
  median_qi()
post.mu <- read.csv(here::here("data/cleaned/posteriors", "posteriors_population.csv")) %>%
  median_qi()

meta <- read.csv(here::here("data/", "lob-metadata.csv"))

df <- read.csv(here::here("data/cleaned/posteriors", "posteriors_individuals.csv")) %>% 
  left_join(meta) %>%
  filter(id != "N07") %>%
  group_by(id) %>%
  sample_draws(500)


forplot <- rbind(post_hoc_CI(mr = unique(df$mr)[1], data = post.h), post_hoc_CI(mr = unique(df$mr)[2], data = post.h), post_hoc_CI(mr = unique(df$mr)[3], data = post.h))
forplot$mr <- rep(unique(df$mr), each = 100)
forplot$treatment <- rep(unique(df$treatment), each = 100)
names(forplot)[1:2] <- c("mc", "h")

sum.forplot <- df %>% ungroup() %>% group_by(mc, treatment) %>%
median_qi(h, .width = c(0.75))


p1 <- ggplot(df, aes(x = mc, y = h))+
  geom_jitter(aes(color = treatment), shape = 1, alpha = 0.5, show.legend = F)+
  geom_pointinterval(data = sum.forplot, aes(x = mc, y = h), color = "gray")+
  geom_point(data = sum.forplot, aes(x = mc, y = h), color = "black")+
  geom_line(data = forplot, aes(x = mc, y = h, color = treatment), show.legend = F)+
  scale_color_manual(values = c('#AF8DC3','#C3AF8D','#8DC3AF'))+
  geom_ribbon(data = forplot, aes(ymin = mu.lower, ymax = mu.upper), fill = "gray", alpha = 0.5)+
  scale_y_log10()+
  facet_wrap(~treatment)+
  labs(x = "Predator body mass (g)", y = "Handling time (h)", color = "")

forplot.a <- rbind(post_hoc_CI(mr = unique(df$mr)[1], data = post.a), post_hoc_CI(mr = unique(df$mr)[2], data = post.a), post_hoc_CI(mr = unique(df$mr)[3], data = post.a))
forplot.a$mr <- rep(unique(df$mr), each = 100)
forplot.a$treatment <- rep(unique(df$treatment), each = 100)
names(forplot.a)[1:2] <- c("mc", "a")

sum.forplot.a <- df %>% ungroup() %>% group_by(mc, treatment) %>%
  median_qi(a, .width = c(0.75))


p2 <- ggplot(df, aes(x = mc, y = a))+
  geom_jitter(aes(color = treatment), shape = 1, alpha = 0.5, show.legend = F)+
  geom_pointinterval(data = sum.forplot.a, aes(x = mc, y = a), color = "gray")+
  geom_point(data = sum.forplot.a, aes(x = mc, y = a), color = "black")+
  geom_line(data = forplot.a, aes(x = mc, y = a, color = treatment), show.legend = F)+
  scale_color_manual(values = c('#AF8DC3','#C3AF8D','#8DC3AF'))+
  geom_ribbon(data = forplot.a, aes(ymin = mu.lower, ymax = mu.upper), fill = "gray", alpha = 0.5)+
  scale_y_log10()+
  facet_wrap(~treatment)+
  labs(x = "Predator body mass (g)", y = expression(paste("Attack rate (ind. m"^{-2},"h"^{-1},")")), color = "")

plot2 <- cowplot::plot_grid(p2, p1)

ggsave(here::here("figures/", "posthoc-plot2.png"), plot2, width = 10, height = 5)

#---------------------------------------------------------------------
## Scrap
#---------------------------------------------------------------------


# Jensen's inequality exploration
forplot <- rbind(post_hoc_CI(mr = unique(df$mr)[1], data = post.h), post_hoc_CI(mr = unique(df$mr)[2], data = post.h), post_hoc_CI(mr = unique(df$mr)[3], data = post.h))
forplot$mr <- rep(unique(df$mr), each = 100)
forplot$treatment <- rep(unique(df$treatment), each = 100)
names(forplot)[1:2] <- c("mc", "h")

sum.forplot <- df %>% ungroup() %>% group_by(mc, treatment) %>%
  median_qi(h, .width = c(0.75))

# What is the expectation for an average sized lobster foraging on an average sized prey?

mean.lob = df %>% group_by(treatment) %>% summarize(mean = mean(mc))
mean.lob  = mean(df$mc)
mean.urc = df %>% group_by(treatment) %>% summarize(mean = mean(mr))
mean.urc = mean(df$mr)


f_x.bar = df %>% group_by(treatment) %>%
  summarize(mean.lob = mean(mc),
            mean.urc = mean(mr)) %>%
  mutate(f_x.bar = exp(mean(post.h$alpha))*mean.lob^mean(post.h$beta1)*mean.urc^mean(post.h$beta2))


f_x.bar <- mean(exp(post.h$alpha)*mean.lob^post.h$beta1*mean.urc^post.h$beta2)
mean.f_x <- mean(forplot$h)

mean.f_x <- forplot %>% group_by(treatment) %>% summarize(mean.f_x = mean(h))


ggplot(df, aes(x = mc, y = h))+
  geom_jitter(aes(color = treatment), shape = 1, alpha = 0.5)+
  geom_pointinterval(data = sum.forplot, aes(x = mc, y = h), color = "gray")+
  geom_point(data = sum.forplot, aes(x = mc, y = h), color = "black")+
  geom_line(data = forplot, aes(x = mc, y = h, color = treatment))+
  scale_color_manual(values = c('#AF8DC3','#C3AF8D','#8DC3AF'))+
  geom_ribbon(data = forplot, aes(ymin = mu.lower, ymax = mu.upper), fill = "gray", alpha = 0.5)+
  # Add in the unstructured and mu expectations!
  geom_hline(data = post.null, aes(yintercept = h), color = "gray", linetype = "dashed")+
  geom_hline(data = post.null, aes(yintercept = h.lower), color = "gray", linetype = "dashed") +
  geom_hline(data = post.null, aes(yintercept = h.upper), color = "gray", linetype = "dashed") +
  geom_hline(data = post.mu, aes(yintercept = mu.h), color = "gray")+
  geom_hline(data = post.mu, aes(yintercept = mu.h.lower), color = "gray")+
  geom_hline(data = post.mu, aes(yintercept = mu.h.upper), color = "gray")+
  scale_y_log10()+
  facet_wrap(~treatment)

forplot <- rbind(post_hoc_CI(mr = unique(df$mr)[1], data = post.a), post_hoc_CI(mr = unique(df$mr)[2], data = post.a), post_hoc_CI(mr = unique(df$mr)[3], data = post.a))
forplot$mr <- rep(unique(df$mr), each = 100)
forplot$treatment <- rep(unique(df$treatment), each = 100)

sum.forplot <- df %>% ungroup() %>% group_by(mc, treatment) %>%
  median_qi(a, .width = c(0.75))


ggplot(df, aes(x = mc, y = a))+
  geom_jitter(aes(color = treatment), shape = 1)+
  geom_pointinterval(data = sum.forplot, aes(x = mc, y = a), color = "gray", fill = "black")+
  geom_point(data = sum.forplot, aes(x = mc, y = a), color = "black")+
  geom_line(data = forplot, aes(x = mc.seq, y = mu))+
  scale_y_log10()+
  facet_wrap(~treatment)


#---------------------------------------------------------



calculate.CI <- function(mr, mc, prob = 0.95){
  require(rethinking)
  a <- exp(post.a$alpha)*mc^post.a$beta1*mr^post.a$beta2
  h <- exp(post.h$alpha)*mc^post.h$beta1*mr^post.h$beta2
  
  predicted.k <- function(N){a*N/(1 + a*h*N)}
  N.seq <- seq(0, 26)
  mu <- sapply(N.seq, predicted.k)
  
  mu.median <- apply( mu , 2 , median ) # calculate the median predicted value for each N
  mu.PI <- t(apply( mu , 2 , PI , prob=prob )) # calculate the credible interval for each value of N
  
  return(data.frame(N.seq = N.seq, mc = mc, mr = mr, mu = mu.median, mu.lower = mu.PI[,1], mu.upper = mu.PI[,2])) # return a data.frame to organize output!
}

ob <- read.table(here("data/cleaned","loburc_cleaned.csv"), header = T, sep = ",") %>%
  arrange(id, treatment)

forplot <- calculate.CI(mr = min(df$mr), mc = max(df$mc))

out <- list()
mt <- meta %>% drop_na(mc)
for(i in 1:length(mt$id)){
  out[[i]] <- calculate.CI(mr = mt$mr[i], mc = mt$mc[i])
}
# This applies names to each slot in the list.

forplot <- bind_rows(out) %>% left_join(mt)


ggplot(ob, aes(x = initial, y = killed/48))+
  geom_jitter()+
  geom_line(data = forplot, aes(x = N.seq, y = mu, color = as.factor(mc)), show.legend = F)+
  facet_wrap(~treatment)










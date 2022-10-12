source("code/1_setup.R")
source("code/theme.R")

#---------------------------------------------------------------------
## Get data
#--------------------------------------------------------------------

# Posteriors
df.ind <- read.csv(here::here("data/cleaned/posteriors", "allometric_individualSTAN.csv")) %>% as_tibble()
df.pop <- read.csv(here::here("data/cleaned/posteriors", "allometric_populationSTAN.csv")) %>% as_tibble() %>% pivot_wider(names_from = .variable, values_from = .value)

# Summary data 

read.csv(here::here("data/cleaned/posteriors", "allometric_populationSTAN.csv")) %>% as_tibble() %>%
  group_by(.variable) %>%
  median_qi(.value)

# Experimental data
df <- read.table(here("data/cleaned","loburc_cleaned.csv"), header = T, sep = ",") %>%
  drop_na(mc)


# Function to estimate mean and CI's of posteriors

# assumes no variance between individuals...
p_link <- function(N, mc, mr, data = df.pop){
  loga <- with(data, alphaa + beta1a*log(mc) + beta2a*log(mr))
  logh <- with(data, alphah + beta1h*log(mc) + beta2h*log(mr))
  a <- exp(loga)
  h <- exp(logh)
  
  a*N / (1 + a*h*N)
}

allometric_CI <- function(mc, mr, prob = 0.95, ...){
  temp.mr <- as.numeric(mr)
  temp.mc <- as.numeric(mc)
  
  N.vec <- seq(0, 26, length.out = 100)
  p_sim_output <- sapply(N.vec, function(i) p_link(i, mc = temp.mc, mr = temp.mr))
  p_mu <- apply(p_sim_output, 2 ,median, na.rm = T)
  p_ci <- t(apply( p_sim_output , 2 , PI, prob = prob))
  
  return(data.frame(N = N.vec, mu = p_mu, mu.lower = p_ci[,1], mu.upper = p_ci[,2], mc = temp.mc, mr = temp.mr))
}

predict <- expand.grid(mc = c(median(df$mc, na.rm = T), quantile(df$mc, probs = c(0.1, 0.9), na.rm = T)), mr = unique(df$mr))

forplot <- predict %>%
  purrr::pmap(allometric_CI) %>%
  set_names(1:9) %>%
  do.call(rbind, .) %>%
  rename(initial = N)

forplot$treatment.cleaned <- forcats::fct_rev(as.factor(rep(c("Medium", "Large", "Small"), each = 300)))
forplot <- arrange(forplot, mc)
forplot$lob.sizeclass <- rep(rev(c("Large", "Medium", "Small")), each = 300)

df <- df %>% 
  mutate(treatment.cleaned = case_when(treatment == "urc_small" ~ "Small", 
                                       treatment == "urc_medium" ~ "Medium", 
                                       treatment == "urc_large" ~ "Large"))


plot4 <- ggplot(df, aes(x = initial, y = killed))+
  geom_jitter(show.legend = F, height = 0, alpha = 0.5, size = 2.5)+
  geom_ribbon(data = forplot, aes(ymin = mu.lower, ymax = mu.upper, y = mu, group = lob.sizeclass), alpha = 0.2)+
  geom_line(data = forplot, aes(x = initial, y = mu, color = lob.sizeclass, linetype = lob.sizeclass), size = 1.5)+
  scale_color_manual(values = rev(c('#d53e4f','#fc8d59','#fee08b')))+
  facet_wrap(~forcats::fct_rev(treatment.cleaned))+
  labs(x = "Number of urchins offered", y = "Number of urchins consumed", color = "Lobster size", title = "Urchin size class", linetype = "Lobster size")+
  theme_bd()+
  theme(legend.position = c(0.8, 0.8))

ggsave(here::here("figures/", "allometric_fr.png"), plot4, width = 10, height = 6)
ggsave(here::here("figures/", "allometric_fr.svg"), plot4, width = 10, height = 6)

ggplot(forplot, aes(x = initial, y = mu))+
  geom_line(aes(color = lob.sizeclass))+
  facet_wrap(~treatment)

#-------------------------------------------------------------------------
## Posterior vs. prior plots
#-------------------------------------------------------------------------


# Handling times
prior.beta1.h <- rnorm(length(df.pop$beta1h), mean = -0.75, sd = rgamma(length(df.pop$beta1h), shape = 1, rate = 1))
plot(density(prior.beta1.h))


prior.beta2.h <- rnorm(length(df.pop$beta2h), mean = 0.5, sd = rgamma(length(df.pop$beta2h), shape = 1, rate = 1))
plot(density(prior.beta2.h))

prior.alpha.h <- rnorm(length(df.pop$alphah), mean = 0, sd = rgamma(length(df.pop$alphah), shape = 1, rate = 1))
plot(density(prior.alpha.h))

# Attack rates
prior.beta1.a <- rnorm(length(df.pop$beta1a), mean = 0.75, sd = rgamma(length(df.pop$beta1a), shape = 1, rate = 1))
plot(density(prior.beta1.a))


prior.beta2.a <- rnorm(length(df.pop$beta2a), mean = 0.5, sd = rgamma(length(df.pop$beta2a), shape = 1, rate = 1))
plot(density(prior.beta2.a))

prior.alpha.a <- rnorm(length(df.pop$alphaa), mean = 0, sd = rgamma(length(df.pop$alphaa), shape = 1, rate = 1))
plot(density(prior.alpha.a))

df.plot <- data.frame(prior.beta1.h, 
                      posterior.beta1.h = df.pop$beta1h, 
                      prior.beta2.h, 
                      posterior.beta2.h = df.pop$beta2h,
                      posterior.alpha.h = df.pop$alphah,
                      prior.alpha.h,
                      prior.beta1.a, 
                      prior.beta2.a, 
                      prior.alpha.a,
                      posterior.beta1.a = df.pop$beta1a, 
                      posterior.beta2.a = df.pop$beta2a, 
                      posterior.alpha.a = df.pop$alphaa) %>%
  gather(dist, value) %>%
  separate(dist, into = c("model", "scaling_exponent", "parameter"), sep = "[.]")%>%
  mutate(parameter = paste(scaling_exponent, parameter, sep = " "))

# meta <- df.plot %>% group_by(model, parameter) %>% 
#   summarize(mean = mean(value))
# meta$text <- c(NA, NA, "MTE\n-0.75", "FP\n0-1")



p5 <- ggplot(df.plot, aes(x = value, fill = model))+
  geom_density(alpha = 0.5)+
  facet_wrap( ~ parameter, scales = "free", nrow = 3)+
  labs(x = "Parameter estimate", y = "Density", fill = "")+
  theme_classic()
  
  ggsave(here::here("figures/", "posterior-prior.png"), p5, width = 8.5, height = 8.5)


#-------------------------------------------------------------------------
## Individual effects
#-------------------------------------------------------------------------  

df.ind <- df.ind %>%
  pivot_wider(names_from = .variable, values_from = .value)

CI_ind <- function(a, h, prob = 0.95, ...){
    N.vec <- seq(0, 26, length.out = 100)
    p_sim_output <- sapply(N.vec, function(i) a*i / (1 + a*h*i))
    p_mu <- apply(p_sim_output, 2 ,median, na.rm = T)
    p_ci <- t(apply( p_sim_output , 2 , PI, prob = prob))
    
    return(data.frame(N = N.vec, mu = p_mu, mu.lower = p_ci[,1], mu.upper = p_ci[,2], id = ids[i]))
  }
  
ids <- unique(df.ind$id)

out <- list()
for(i in 1:length(ids)){
  temp.a <- df.ind$a[df.ind$id == ids[i]]
  temp.h <- df.ind$h[df.ind$id == ids[i]]
  out[[i]] <- CI_ind(a = temp.a, h = temp.h)
}

formerge <- df %>% select(id, treatment, mc)
forplot <- data.frame(do.call(rbind, out)) %>% left_join(formerge) %>%
  mutate(id.ord = forcats::fct_reorder(id, mc, .desc = T))

df %>% 
  mutate(id.ord = forcats::fct_reorder(id, mc, .desc = T)) %>%
ggplot(aes(x = initial, y = killed))+
  geom_jitter(aes(size = treatment), pch = 21, show.legend = F, stroke = 1)+
  scale_size_manual(values = c(4,2.5,1.5))+
  geom_line(data = forplot, aes(x = N, y = mu))+
  geom_ribbon(data = forplot, aes(x = N, y = mu, ymin = mu.lower, ymax = mu.upper), alpha = 0.5)+
  facet_wrap(~id)+
  labs(x = "Initial prey abundance", y = "Number of prey consumed")+
  theme_classic()


p1 <- df %>% 
  filter(treatment == "urc_large") %>%
  mutate(id.ord = forcats::fct_reorder(id, mc, .desc = T)) %>%
  ggplot(aes(x = initial, y = killed))+
  geom_jitter(pch = 21, size = 4)+
  geom_line(data = forplot %>% filter(treatment == "urc_large"), aes(x = N, y = mu))+
  geom_ribbon(data = forplot %>% filter(treatment == "urc_large"), aes(x = N, y = mu, ymin = mu.lower, ymax = mu.upper), alpha = 0.5)+
  facet_wrap(~id.ord, ncol = 4)+
  coord_cartesian(ylim = c(0, 26))+
  labs(x = "Initial prey abundance", y = "Number of prey consumed", title = "Large urchins")+
  theme_classic()

p2 <- df %>% 
  filter(treatment == "urc_medium") %>%
  mutate(id.ord = forcats::fct_reorder(id, mc, .desc = T)) %>%
  ggplot(aes(x = initial, y = killed))+
  geom_jitter(pch = 21, size = 2.5)+
  geom_line(data = forplot %>% filter(treatment == "urc_medium"), aes(x = N, y = mu))+
  geom_ribbon(data = forplot %>% filter(treatment == "urc_medium"), aes(x = N, y = mu, ymin = mu.lower, ymax = mu.upper), alpha = 0.5)+
  facet_wrap(~id.ord, ncol = 4)+
  coord_cartesian(ylim = c(0, 26))+
  labs(x = "Initial prey abundance", y = "Number of prey consumed", title = "Medium urchins")+
  theme_classic()

p3 <- df %>% 
  filter(treatment == "urc_small") %>%
  mutate(id.ord = forcats::fct_reorder(id, mc, .desc = T)) %>%
  ggplot(aes(x = initial, y = killed))+
  geom_jitter(pch = 21, size = 1.5)+
  geom_line(data = forplot %>% filter(treatment == "urc_small"), aes(x = N, y = mu))+
  geom_ribbon(data = forplot %>% filter(treatment == "urc_small"), aes(x = N, y = mu, ymin = mu.lower, ymax = mu.upper), alpha = 0.5)+
  facet_wrap(~id.ord, ncol = 4)+
  coord_cartesian(ylim = c(0, 26))+
  labs(x = "Initial prey abundance", y = "Number of prey consumed", title = "Small urchins")+
  theme_classic()

cowplot::plot_grid(p1+ labs(x = "", y = ""), p2 + labs(x = ""), p3 + labs(y = ""), nrow = 3)
ggsave("figures/ind_FR.png", width = 8.5*1.5, height = 11*1.5)

#-------------------------------------------------------------------------
## Allometric scaling of parameters
#-------------------------------------------------------------------------  
background_points <- read.csv(here::here("data/cleaned/posteriors", "allometric_individualSTAN.csv")) %>% as_tibble() %>%
  pivot_wider(names_from = .variable, values_from = .value) %>%
  mutate(a = a/48, 
         h = h*48) %>%
  left_join(formerge) %>%
  group_by(id) %>%
  tidybayes::sample_draws(n = 100)



urchin_cols <- rev(c("#65aecd","#436fe2","#7b47c1"))

plot(1:3, col = urchin_cols, pch = 19, cex = 5)

post_hoc_CI <- function(mr, data, prob = 0.95){
  temp.mr <- as.numeric(mr)
  mu.link <- function(mc, mr){
    with(data, exp(alpha)*mc^beta1*mr^beta2)
  } # defines a function to predict the prey killed at combination of a and h in the posteriors
  
  mc.seq <- seq( from=min(df[, "mc"]) , to=max(df[, "mc"]) , length.out = 100) #define a sequence of consumer masses
  mu <- sapply( mc.seq, function(i) mu.link(i, mr = temp.mr)) # apply the mu.link funciton to each N in the sequence
  
  mu.median <- apply( mu , 2 , median ) # calculate the median predicted value for each N
  mu.PI <- t(apply( mu , 2 , PI , prob=0.95 )) # calculate the credible interval for each value of N
  
  return(data.frame(mc.seq = mc.seq, mu = mu.median, mu.lower = mu.PI[,1], mu.upper = mu.PI[,2]))
}


df.pop.a <- df.pop %>%
  select(alphaa, beta1a, beta2a) %>%
  rename(alpha = alphaa, beta1 = beta1a, beta2 = beta2a) %>%
  mutate(parameter = "a")

forplot.a <- rbind(post_hoc_CI(mr = unique(df$mr)[1], data = df.pop.a), post_hoc_CI(mr = unique(df$mr)[2], data = df.pop.a), post_hoc_CI(mr = unique(df$mr)[3], data = df.pop.a))
forplot.a$mr <- rep(unique(df$mr), each = 100)
forplot.a$treatment <- rep(unique(df$treatment), each = 100)
forplot.a$parameter <- "a"
names(forplot.a)[1] <- c("mc")

df.pop.h <- df.pop %>%
  select(alphah, beta1h, beta2h) %>%
  rename(alpha = alphah, beta1 = beta1h, beta2 = beta2h) %>%
  mutate(parameter = "h")

forplot.h <- rbind(post_hoc_CI(mr = unique(df$mr)[1], data = df.pop.h), post_hoc_CI(mr = unique(df$mr)[2], data = df.pop.h), post_hoc_CI(mr = unique(df$mr)[3], data = df.pop.h))
forplot.h$mr <- rep(unique(df$mr), each = 100)
forplot.h$treatment <- rep(unique(df$treatment), each = 100)
forplot.h$parameter <- "h"
names(forplot.h)[1] <- c("mc")


forplot <- bind_rows(forplot.a, forplot.h) %>% 
  as_tibble() %>%
  mutate(mu = case_when(parameter == "a" ~ mu/48, 
                        parameter == "h" ~ mu*48), 
         mu.lower = case_when(parameter == "a" ~ mu.lower/48, 
                              parameter == "h" ~ mu.lower*48), 
         mu.upper = case_when(parameter == "a" ~ mu.upper/48, 
                              parameter == "h" ~ mu.upper*48))

p1 <- background_points %>% 
  ggplot(aes(x = mc, y = a))+
  geom_point(aes(color = treatment), shape = 1, alpha = 0.1, show.legend = T)+
  geom_ribbon(data = forplot[forplot$parameter == "a", ], aes(ymin = mu.lower, ymax = mu.upper, y = mu, x = mc, group = treatment), fill = "gray", alpha = 0.25)+
  geom_line(data = forplot[forplot$parameter == "a",], aes( x= mc, y = mu, color = treatment, linetype = treatment), show.legend = T, lwd = 2)+
  scale_color_manual(values = urchin_cols)+
  scale_y_log10()+
  scale_x_log10()+
  labs(x = "Predator body mass (g)", y = expression(paste("Attack rate (",m^-2,h^-1,")")), color = "", linetype = "")+
  cowplot::theme_cowplot()+
  theme(legend.position = "none")

p2 <- background_points %>% 
  ggplot(aes(x = mc, y = h))+
  geom_point(aes(color = treatment), shape = 1, alpha = 0.1, show.legend = T)+
  geom_ribbon(data = forplot[forplot$parameter == "h", ], aes(ymin = mu.lower, ymax = mu.upper, y = mu, x = mc, group = treatment), fill = "gray", alpha = 0.25)+
  geom_line(data = forplot[forplot$parameter == "h",], aes( x= mc, y = mu, color = treatment, linetype = treatment), show.legend = T, lwd = 2)+
  scale_color_manual(values = urchin_cols)+
  scale_y_log10()+
  scale_x_log10()+
  labs(x = "Predator body mass (g)", y = "Handling time (h)", color = "", linetype = "")+
  cowplot::theme_cowplot()+
  theme(legend.position = c(0.6, 0.8))

fig_s2 <- cowplot::plot_grid(p1, p2, labels = "AUTO")

ggsave(here::here("figures/", "fig4_posthocandh.png"), fig_s2, width = 10, height = 5)

fig_s2 <- cowplot::plot_grid(p1, p2)

ggsave("figures/fig4_posthocandh.svg", fig_s2, device = "svg", width = 10, height = 5)






#-------------------------------------------
## Setup
#-------------------------------------------

library(here)
source(here("code", "Base_functions/Functions.R"))
source(here("code", "1_setup.R"))

#-------------------------------------------
## Get data
#-------------------------------------------

post.a <- read.csv(here::here("data/cleaned/posteriors", "posteriors_posthoc_a.csv"))
post.h <- read.csv(here::here("data/cleaned/posteriors", "posteriors_posthoc_h.csv"))

post.allometric <- read.csv(here::here("data/cleaned/posteriors", "allometric_population.csv"))

    # Summary stats for paper
      
      post.a %>% pivot_longer(cols = c(alpha:beta2)) %>%
        group_by(name) %>%
        median_qi(value)
      
      post.h %>% pivot_longer(cols = c(alpha:beta2)) %>%
        group_by(name) %>%
        median_qi(value)
      
      post.allometric %>% pivot_longer(cols = c(mu.alpha.h:beta2.h)) %>%
        group_by(name) %>%
        median_qi(value)



meta <- read.csv(here::here("data/", "lob-metadata.csv"))

df <- read.csv(here::here("data/cleaned/posteriors", "posteriors_individuals.csv")) %>% 
  left_join(meta) %>%
  filter(id != "N07") %>%
  group_by(id) %>%
  sample_draws(100)


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
  scale_x_log10()+
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
  scale_x_log10()+
  facet_wrap(~treatment)+
  labs(x = "Predator body mass (g)", y = expression(paste("Attack rate (ind. m"^{-2},"h"^{-1},")")), color = "")

plot2 <- cowplot::plot_grid(p2, p1)

ggsave(here::here("figures/", "posthoc-plot2.png"), plot2, width = 10, height = 5)



# Alternate version of figure

forplot <- forplot %>% mutate(parameter = "h", 
                   mu = h, 
                   h = NULL) %>% 
  bind_rows(forplot.a %>% mutate(parameter = "a", mu = a, a = NULL)) %>% 
  as_tibble()

sum.forplot <- df %>% 
  pivot_longer(cols = c(a,h), names_to = "parameter", values_to = "posterior_estimate") %>%
  group_by(parameter, mc, treatment) %>%
  median_qi(posterior_estimate, .width = c(0.75))

df %>% 
  pivot_longer(cols = c(a,h), names_to = "parameter", values_to = "posterior_estimate") %>%
  ggplot(aes(x = mc, y = posterior_estimate))+
  geom_point(aes(color = treatment), shape = 1, alpha = 0.5, show.legend = T)+
  geom_point(data = sum.forplot, aes(x = mc, y = posterior_estimate), pch = 21, fill = "gray30")+
  geom_line(data = forplot, aes( x= mc, y = mu, color = treatment), show.legend = T)+
  scale_color_manual(values = c('#AF8DC3','#C3AF8D','#8DC3AF'))+
  geom_ribbon(data = forplot, aes(ymin = mu.lower, ymax = mu.upper, y = mu, x = mc, group = treatment), fill = "gray", alpha = 0.5)+
  scale_y_log10()+
  scale_x_log10()+
  facet_wrap(~ parameter, scales = "free")+
  labs(x = "Predator body mass (g)", y = "Parameter", color = "")+
  cowplot::theme_cowplot()

# Version w/out facet wrap and correct axis labels

p1 <- df %>% 
  ggplot(aes(x = mc, y = a))+
  geom_point(aes(color = treatment), shape = 1, alpha = 0.5, show.legend = T)+
  geom_point(data = sum.forplot[sum.forplot$parameter == "a",], aes(x = mc, y = posterior_estimate), pch = 21, fill = "gray30")+
  geom_line(data = forplot[forplot$parameter == "a",], aes( x= mc, y = mu, color = treatment), show.legend = T)+
  scale_color_manual(values = c('#AF8DC3','#C3AF8D','#8DC3AF'))+
  geom_ribbon(data = forplot[forplot$parameter == "a", ], aes(ymin = mu.lower, ymax = mu.upper, y = mu, x = mc, group = treatment), fill = "gray", alpha = 0.5)+
  scale_y_log10()+
  scale_x_log10()+
  labs(x = "Predator body mass (g)", y = expression(paste("Attack rate (",m^-2,h^-1,")")), color = "")+
  cowplot::theme_cowplot()+
  theme(legend.position = c(0.6, 0.2))

p2 <- df %>% 
  ggplot(aes(x = mc, y = h))+
  geom_point(aes(color = treatment), shape = 1, alpha = 0.5, show.legend = T)+
  geom_point(data = sum.forplot[sum.forplot$parameter == "h",], aes(x = mc, y = posterior_estimate), pch = 21, fill = "gray30")+
  geom_line(data = forplot[forplot$parameter == "h",], aes( x= mc, y = mu, color = treatment), show.legend = T)+
  scale_color_manual(values = c('#AF8DC3','#C3AF8D','#8DC3AF'))+
  geom_ribbon(data = forplot[forplot$parameter == "h", ], aes(ymin = mu.lower, ymax = mu.upper, y = mu, x = mc, group = treatment), fill = "gray", alpha = 0.5)+
  scale_y_log10()+
  scale_x_log10()+
  labs(x = "Predator body mass (g)", y = "Handling time (h)", color = "")+
  cowplot::theme_cowplot()+
  theme(legend.position = "none")

fig_s2 <- cowplot::plot_grid(p1, p2)

ggsave(here::here("figures/", "fig4_posthocandh.png"), fig_s2, width = 10, height = 5)


##-----------------------------
# Posterior-prior comparisons
##-----------------------------


# On handling time
nsim = 5000

sigma.beta1 <- rgamma(nsim, 2,2)
sigma.beta2 <- rgamma(nsim, 2,2)
beta1 <- rnorm(nsim, mean = -0.75, sd = sigma.beta1)
hist(beta1) # prior predictive distribution
beta2 <- rnorm(nsim, mean = 0.50, sd = sigma.beta2 )
alpha <- rnorm(nsim, 0, tau2sd(0.01))

sigmay <- runif(nsim, 0, 10)
tauy <- 1/(sigmay*sigmay)

# mu <- alpha + beta1 + beta2
# handling_time <- rnorm(nsim, mu, tauy)

priors <- data.frame(beta1 = beta1,
                     beta2 = beta2, 
                     alpha = alpha) %>% 
  pivot_longer(cols = c(beta1, beta2, alpha), names_to = "parameter", values_to = "prior_estimate")

post.h %>% 
  sample_draws(nsim) %>%
  pivot_longer(cols = c(beta1, beta2, alpha), names_to = "parameter", values_to = "posterior_estimate") %>%
  bind_cols(select(priors, prior_estimate)) %>%
  pivot_longer(cols = c(posterior_estimate, prior_estimate)) %>%
  mutate(name = case_when(
    name == "posterior_estimate" ~ "posterior", 
    name == "prior_estimate" ~ "prior"
  ), 
  name = forcats::fct_rev(name)) %>%
  ggplot(aes(x = value))+
  geom_density(aes(fill = name, y = ..scaled..), alpha = 1/2, color = NA)+
  facet_wrap(~parameter, scales = "free")+
  theme_classic()


# On attack rate
nsim = 5000

sigma.beta1 <- rgamma(nsim, 2,2)
sigma.beta2 <- rgamma(nsim, 2,2)
beta1 <- rnorm(nsim, mean = beta1.theory, sd = sigma.beta1)
hist(beta1) # prior predictive distribution
beta2 <- rnorm(nsim, mean = beta2.theory, sd = sigma.beta2 )
alpha <- rnorm(nsim, 0, tau2sd(0.01))

sigmay <- runif(nsim, 0, 10)
tauy <- 1/(sigmay*sigmay)

priors <- data.frame(beta1 = beta1,
                     beta2 = beta2, 
                     alpha = alpha) %>% 
  pivot_longer(cols = c(beta1, beta2, alpha), names_to = "parameter", values_to = "prior_estimate")

post.a %>% 
  sample_draws(nsim) %>%
  pivot_longer(cols = c(beta1, beta2, alpha), names_to = "parameter", values_to = "posterior_estimate") %>%
  bind_cols(select(priors, prior_estimate)) %>%
  pivot_longer(cols = c(posterior_estimate, prior_estimate)) %>%
  mutate(name = case_when(
    name == "posterior_estimate" ~ "posterior", 
    name == "prior_estimate" ~ "prior"
  ), 
  name = forcats::fct_rev(name)) %>%
  ggplot(aes(x = value))+
  geom_density(aes(fill = name, y = ..scaled..), alpha = 1/2, color = NA)+
  facet_wrap(~parameter, scales = "free")+
  theme_classic()


















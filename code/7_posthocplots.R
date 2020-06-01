

post.a <- read.csv(here::here("data/cleaned", "posteriors_posthoc_a.csv")) %>%
  spread(parameter,estimate) %>% filter(model == "model_wpuncert")
post.h <- read.csv(here::here("data/cleaned", "posteriors_posthoc_h.csv")) %>%
  spread(parameter,estimate) %>% filter(model == "model_wpuncert")

meta <- read.csv(here::here("data/", "lob-metadata.csv"))
df <- read.csv(here::here("data/cleaned", "posteriors_individuals.csv")) %>% 
  left_join(meta) %>%
  filter(id != "N07") %>%
  group_by(id) %>%
  sample_draws(500)


temp <- calculate.CI(mr = mean(df$mr), mc = mean(df$mc), alpha.h = post.h$alpha, alpha.a = post.a$alpha, beta1.h = post.h$beta1, beta1.a = post.a$beta1, beta2.h = post.h$beta2, beta2.a = post.h$beta2)

ggplot(temp, aes(x = N.seq, y = mu))+
  geom_point()

ggplot(df, aes( x = initial, y = killed))+
  geom_jitter()+
  geom_line(data = temp, aes(x = N.seq, y = mu*48))





post_hoc_CI <- function(mr, data){
  temp.mr <- as.numeric(mr)
mu.link <- function(mc, mr = temp.mr, alpha = data[, "alpha"], beta1 = data[, "beta1"], beta2 = data[, "beta2"]){
  exp(alpha)*mc^beta1*mr^beta2
} # defines a function to predict the prey killed at combination of a and h in the posteriors

mc.seq <- seq( from=min(df$mc) , to=max(df$mc) , length.out = 100) # define a sequence of initial densities
mu <- sapply( mc.seq, mu.link) # apply the mu.link funciton to each N in the sequence

mu.median <- apply( mu , 2 , median ) # calculate the median predicted value for each N
mu.PI <- t(apply( mu , 2 , PI , prob=prob )) # calculate the credible interval for each value of N

return(data.frame(mc.seq = mc.seq, mu = mu.median, mu.lower = mu.PI[,1], mu.upper = mu.PI[,2]))
}


forplot <- rbind(post_hoc_CI(mr = unique(df$mr)[1], data = post.h), post_hoc_CI(mr = unique(df$mr)[2], data = post.h), post_hoc_CI(mr = unique(df$mr)[3], data = post.h))
forplot$mr <- rep(unique(df$mr), each = 100)
forplot$treatment <- rep(unique(df$treatment), each = 100)

sum.forplot <- df %>% ungroup() %>% group_by(mc, treatment) %>%
median_qi(h, .width = c(0.75))


ggplot(df, aes(x = mc, y = h))+
  geom_jitter(aes(color = treatment), shape = 1)+
  geom_pointinterval(data = sum.forplot, aes(x = mc, y = h), color = "gray", fill = "black")+
  geom_point(data = sum.forplot, aes(x = mc, y = h), color = "black")+
  geom_line(data = forplot, aes(x = mc.seq, y = mu))+
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


























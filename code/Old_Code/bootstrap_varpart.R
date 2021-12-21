
variance.partition <- function(ndraws = 1000, .urcdensity, .lobdensity, .urcmass, .lobmass){
  full_meandensity <- r.s %>%
    ungroup() %>%
    mutate(urc_density = .urcdensity,
           lob_density = .lobdensity) %>%
    purrr::pmap(allometricFR, 
                a0. = post.a$alpha, 
                h0. = post.h$alpha, 
                beta1a. = post.a$beta1, 
                beta2a. = post.a$beta2, 
                beta1h. = post.h$beta1, 
                beta2h. = post.h$beta2) %>% 
    purrr::flatten() %>%
    set_names(names) %>%
    as_tibble() %>%
    gather(id, prediction) %>%
    separate(id, into = c("year", "site"), sep = "[-]") %>%
    mutate(prediction = prediction*24/tsize,
           estimate = "full_meandensity")
  
  full_meanbodysize <- r.s %>%
    ungroup() %>%
    select(-urc_mass, -lob_mass) %>% 
    mutate(urc_mass = purrr::map(.urcmass, rep, ndraws), 
           lob_mass = purrr::map(.lobmass, rep, ndraws))%>%
    purrr::pmap(allometricFR, 
                a0. = post.a$alpha, 
                h0. = post.h$alpha, 
                beta1a. = post.a$beta1, 
                beta2a. = post.a$beta2, 
                beta1h. = post.h$beta1, 
                beta2h. = post.h$beta2) %>% 
    set_names(names) %>%
    as_tibble() %>%
    gather(id, prediction) %>%
    separate(id, into = c("year", "site"), sep = "[-]") %>%
    mutate(prediction = prediction*24/tsize,
           estimate = "full_meanbodysize")
  
  # Due to bodysize
  duetobodysize <- var(full_meandensity$prediction) / (var(full_meandensity$prediction) + var(full_meanbodysize$prediction))
  
  # Due to density
  duetodensity <- var(full_meanbodysize$prediction) / (var(full_meandensity$prediction) + var(full_meanbodysize$prediction))
  
  out <- list(duetobodysize, duetodensity)
  names(out) <- c("duetobodysize", "duetodensity")
  
  out
}


# test the function
variance.partition(ndraws = 1000, .urcdensity = mean_urc_density, .lobdensity = mean_lob_density, .urcmass = mean_urc_mass, .lobmass = mean_lob_mass) # works

system.time(t <- variance.partition(ndraws = 1000, .urcdensity = mean_urc_density, .lobdensity = mean_lob_density, .urcmass = mean_urc_mass, .lobmass = mean_lob_mass)) 


summary(unnest(r.s, cols = c(urc_mass))$mass)

pred <- expand.grid(urcdensity = seq(min(r.s$urc_density), max(r.s$urc_density), length.out = 10), 
                    lobdensity = seq(min(r.s$lob_density), max(r.s$lob_density), length.out = 10), 
                    urcmass = seq(min(unnest(r.s, cols = c(urc_mass))$mass), max(unnest(r.s, cols = c(urc_mass))$mass), length.out = 10), 
                    lobmass = seq(min(unnest(r.s, cols = c(lob_mass))$mass), max(unnest(r.s, cols = c(lob_mass))$mass), length.out = 10)) %>% 
  mutate(.id = 1:n()) %>%
  tidybayes::sample_draws(n = 1000, draw = ".id")

out <- vector()
system.time(
for(i in 1:dim(pred)[1]){
  out[i] <- variance.partition(ndraws = 100, 
                     .urcdensity = pred$urcdensity[i], 
                     .lobdensity = pred$lobdensity[i], 
                     .urcmass =  pred$urcmass[i], 
                     .lobmass = pred$lobmass[i])$duetobodysize
})

hist(out) # This shows that the technique is entirely dependent on the fixed values chosen and has no real bearing on how much of the total variation is due to bodysize vs. density. 


mean(out)




#------------------------------------------------------------------------

variance.partition <- function(ndraws = 1000, .urcdensity, .lobdensity, .urcmass, .lobmass){
  
  full <- r.s %>%
    purrr::pmap(allometricFR,
                a0. = post.a$alpha,
                h0. = post.h$alpha,
                beta1a. = post.a$beta1,
                beta2a. = post.a$beta2,
                beta1h. = post.h$beta1,
                beta2h. = post.h$beta2) %>%
    purrr::flatten() %>%
    set_names(names) %>%
    as_tibble() %>%
    gather(id, prediction) %>%
    separate(id, into = c("year", "site"), sep = "[-]") %>%
    mutate(prediction = prediction*24/tsize,
           estimate = "full")
  
  full_meandensity <- r.s %>%
    ungroup() %>%
    mutate(urc_density = .urcdensity,
           lob_density = .lobdensity) %>%
    purrr::pmap(allometricFR, 
                a0. = post.a$alpha, 
                h0. = post.h$alpha, 
                beta1a. = post.a$beta1, 
                beta2a. = post.a$beta2, 
                beta1h. = post.h$beta1, 
                beta2h. = post.h$beta2) %>% 
    purrr::flatten() %>%
    set_names(names) %>%
    as_tibble() %>%
    gather(id, prediction) %>%
    separate(id, into = c("year", "site"), sep = "[-]") %>%
    mutate(prediction = prediction*24/tsize,
           estimate = "full_meandensity")
  
  full_meanbodysize <- r.s %>%
    ungroup() %>%
    select(-urc_mass, -lob_mass) %>% 
    mutate(urc_mass = purrr::map(.urcmass, rep, ndraws), 
           lob_mass = purrr::map(.lobmass, rep, ndraws))%>%
    purrr::pmap(allometricFR, 
                a0. = post.a$alpha, 
                h0. = post.h$alpha, 
                beta1a. = post.a$beta1, 
                beta2a. = post.a$beta2, 
                beta1h. = post.h$beta1, 
                beta2h. = post.h$beta2) %>% 
    set_names(names) %>%
    as_tibble() %>%
    gather(id, prediction) %>%
    separate(id, into = c("year", "site"), sep = "[-]") %>%
    mutate(prediction = prediction*24/tsize,
           estimate = "full_meanbodysize")
  
  df3 <- rbind(full, full_meandensity, full_meanbodysize) %>% 
    group_by(year, site, estimate) %>%
    mutate(.id = 1:n()) %>%
    pivot_wider(id_cols = c(year, site, .id), names_from = estimate, values_from = prediction) %>%
    mutate(.id = NULL, 
           ID = paste(year, site, sep = "-"))
  
  lm.full <- lm(full ~ full_meandensity + full_meanbodysize, df3)
  summary(lm.full)
  
  lm.1 <- lm(full ~ full_meandensity, df3)
  summary(lm.1)
  
  lm.2 <- lm(full ~ full_meanbodysize, df3)
  summary(lm.2)
  
  out <- vector()
  
  #proportion of variation due to body size 
  out[1] <- summary(lm.full)$r.squared - summary(lm.2)$r.squared
  
  #proportion of variation due to density
  out[2] <- summary(lm.full)$r.squared - summary(lm.1)$r.squared
  
  #proportion of variation shared between body size and density
  out[3] <- summary(lm.1)$r.squared + summary(lm.2)$r.squared - summary(lm.full)$r.squared
  
  # total variation explained by the model
  out[4] <- summary(lm.full)$r.squared
  
  # proportion of variation due to body size
  out[5] <- summary(lm.1)$r.squared
  
  # proportion of variation due to density
  out[6] <- summary(lm.2)$r.squared
  
  out
}

pred <- expand.grid(urcdensity = seq(min(r.s$urc_density), max(r.s$urc_density), length.out = 10), 
                    lobdensity = seq(min(r.s$lob_density), max(r.s$lob_density), length.out = 10), 
                    urcmass = seq(min(unnest(r.s, cols = c(urc_mass))$mass), max(unnest(r.s, cols = c(urc_mass))$mass), length.out = 10), 
                    lobmass = seq(min(unnest(r.s, cols = c(lob_mass))$mass), max(unnest(r.s, cols = c(lob_mass))$mass), length.out = 10)) %>% 
  mutate(.id = 1:n()) %>%
  tidybayes::sample_draws(n = 10, draw = ".id")

out <- matrix(nrow = dim(pred)[1], ncol = 6)
system.time(
  for(i in 1:dim(pred)[1]){
    out[i, ] <- variance.partition(ndraws = 100, 
                                 .urcdensity = pred$urcdensity[i], 
                                 .lobdensity = pred$lobdensity[i], 
                                 .urcmass =  pred$urcmass[i], 
                                 .lobmass = pred$lobmass[i])
  })


output <- data.frame(out)
names(output) <- c("duetobodysize", "duetodensity", "shared", "model_total", "new_duetobodysize", "new_duetodensity")
output$prop_duetobodyszie <- output$duetobodysize/output$model_total
output$prop_duetodensity <- output$duetodensity/output$model_total

mean(output$prop_duetobodyszie)
mean(output$prop_duetodensity)


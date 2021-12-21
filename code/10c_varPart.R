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
  
  df3 <- rbind(full, full_meanbodysize) %>% 
    group_by(year, site, estimate) %>%
    mutate(.id = 1:n()) %>%
    pivot_wider(id_cols = c(year, site, .id), names_from = estimate, values_from = prediction) %>%
    mutate(.id = NULL, 
           ID = paste(year, site, sep = "-"))

  
  lm <- lm(full ~ full_meanbodysize, df3)
  
  # proportion of variation due to density
  out <- summary(lm)$r.squared
  
  out
}

# From extensive simulation, I know that I can't accurately assess the variance due to bodysize directly because I can't seperate from density. But I can accurately assess the proportion of variance in the full simulation due to density, by fixing body size. Here, I fix body size to different values of predator and prey size across the range, then estimate the R2 of the regression between the full simulation ~ simulations with fixed body size across each different combination of predator and prey size. Hopefully there will be low variance in the R2, suggesting the amount of variance due to body size. 


pred <- expand.grid(urcmass = seq(min(unnest(r.s, cols = c(urc_mass))$mass), max(unnest(r.s, cols = c(urc_mass))$mass), length.out = 25), 
                    lobmass = seq(min(unnest(r.s, cols = c(lob_mass))$mass), max(unnest(r.s, cols = c(lob_mass))$mass), length.out = 25))

out <- vector()
system.time(
  for(i in 1:dim(pred)[1]){
    out[i] <- variance.partition(ndraws = 100, 
                                .urcmass =  pred$urcmass[i], 
                                .lobmass = pred$lobmass[i])
  })


out
summary(out)
hist(out)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.1071  0.1452  0.1860  0.1710  0.1979  0.2019 

# Ok, so depending on which value of predator and prey size the proportion of variance in the full simulation explained by density ranges from ~10-20%. This suggests that variation in body size within and between sites explains ~80-90% of the variance in interaction strength. The conclusion that I draw from this is that failing to account for variation in body size will result in a loss of 80-90% of the actual variation in how strongly species interact. 









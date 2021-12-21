

temp <- r.s %>% left_join(full %>% group_by(year, site) %>% 
                    nest(prediction = prediction) %>% 
                    mutate(year = as.integer(year), 
                           estimate = NULL)) %>% 
  unnest() %>%
  rename(urc_mass = mass, lob_mass = mass1)



lm.full <- lm(prediction ~ lob_mass + urc_mass + lob_density + urc_density, temp)
summary(lm.full)

var(full$prediction) / 
  (var(full$prediction) + var(full$prediction - full_meandensity$prediction))



var(full_meandensity$prediction - full$prediction) / var(full$prediction)
var(full$prediction - rep(full_meanbodysize$prediction, each = 10000)) / var(full$prediction)

summary(full_meandensity$prediction - full$prediction)


var(scale(full$prediction))

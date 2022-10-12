
source(here::here("code", "10a_clean-obsdata.R"))

f=expression((0.6312*mc^0.0498*mr^0.0925 * N) / 
               (1 + 0.6312*mc^0.0498*mr^0.0925 * 357.8*mc^-1.61*mr^1.30 * N ))

D(f, "N")

derivative <- function(mc, mr, N){
  0.6312 * mc^0.0498 * mr^0.0925/(1 + 0.6312 * mc^0.0498 * mr^0.0925 * 
  357.8 * mc^-1.61 * mr^1.3 * N) - (0.6312 * mc^0.0498 * mr^0.0925 * 
  N) * (0.6312 * mc^0.0498 * mr^0.0925 * 357.8 * mc^-1.61 * 
  mr^1.3)/(1 + 0.6312 * mc^0.0498 * mr^0.0925 * 357.8 * mc^-1.61 * mr^1.3 * N)^2
}

mc <- 100:3000
mr <- c(3, 30, 90)
N <- 1:26

temp <- derivative(mc = 3000, mr = 3, N = 1:26)
plot(temp ~ N)


r.s %>% ungroup() %>% summarize(max_urc_mass = map_dbl(urc_mass, ~min(.x$mass)), 
                                min_lob_mass = map_dbl(lob_mass, ~ max(.x$mass)), 
                                min_urc_density = max(urc_density)) %>%
  summary()

plot(derivative(mc = 6184.0, mr = 8.161, N = 1:200) ~ c(1:200))



output <- r.s %>%
  ungroup() %>%
  group_by(year, site, urc_density) %>% 
  summarize(mean_urc_mass = map_dbl(urc_mass, ~mean(.x$mass)), 
            mean_lob_mass = map_dbl(lob_mass, ~mean(.x$mass))) %>%
  mutate(derivative = derivative(mc = mean_lob_mass, mr = mean_urc_mass, N = urc_density),
         max_derivative = derivative(mc = mean_lob_mass, mr = mean_urc_mass, N = 0.001),
         per_of_max = (1-(max_derivative-derivative) / max_derivative)*100)

hist(output$derivative)
hist(output$per_of_max, xlab = "% of maximum rate of change")


# Here I estimate the maximum rate of change in consumption with prey density as the derivative evaluated at ~zero prey density. We know that the derivative of a M-M funciton is always decreasing. Thus the maximum is at zero. I then estimated, given the average urchin and lobster body size at a site, the derivative at the observed urchin density. And compared the rate of change at the observated prey density with the maximum possible rate of change.

# From this I conclude that all site/years had urchin densities where the rate of change in interaction strength with prey density was less than 12 % of the maximum rate of change. Furthermore, I predict that at 60% of the site/years the rate of change in interaction strength with prey density was less than 1% of the maximum at low prey desities. This means that at most sites, its unlikely that losters where limited by prey availability. Rather they were limited by their capacity to handle prey -- a component of predation which we know to be strongly body size dependent.


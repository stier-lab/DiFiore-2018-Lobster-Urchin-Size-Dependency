#-------------------------------------------------------------------------------------------
## Build a map of interaction strengths through time and space
#-------------------------------------------------------------------------------------------



df <- read.csv(here("data/cleaned", "all_densities.csv"))


sz <- read.csv(here("data/cleaned", "size_frequencies.csv")) %>%
  select(-month, -day, -mpa_status, -reserve, -year_mpa, -count, -area) %>%
  mutate(transect = paste(`side.LTERtransect`, `zone.LTERreplicate`, transect, sep = "-")) %>%
  select(-c(`side.LTERtransect`, `zone.LTERreplicate`, transect)) %>%
  group_by(year, site, classcode) %>%
  summarize(median.size = median(size, na.rm = T)) %>%
  mutate(mass = ifelse(classcode == "PANINT", 0.001352821*(median.size*10)^2.913963113,
                       ifelse(classcode == "STRPURAD", 0.000592598*(median.size*10)^2.872636198*1.01, NA)), 
         median.size = NULL
  )%>%
  spread(classcode, mass, fill = 0)

names(sz) <- c("year", "site", "lob.size", "urc.size")

d <- df %>%
  group_by(year, site, mpa_status, reserve, year_mpa) %>%
  summarize(lob.den = mean(PANINT),
            urc.den = mean(STRPURAD)) %>% 
  left_join(sz)

b.hol <- function(n, a, p, t, mc, mr){
  h <- 6.94936*exp(-1.07875*log10((mc/mr)+1))
  a*n*p*t/(1+a*h*n)
}

d$est <- b.hol(n = d$urc.den, a = mean(a[,1]), p = d$lob.den, mc = d$lob.size, mr = d$urc.size, t = 1)


d %>% group_by(reserve) %>%
  summarize(mean(est, na.rm = T))
boxplot(sqrt(est) ~ reserve, d)


lme6 <- lmer(sqrt(est) ~ reserve + (1|year) + (1|site), d[d$year >= 2010, ])
summary(lme6)

summary(lm(sqrt(est) ~ reserve, d))

hist(d$est)

ggplot(d, aes(x = year, y = est)) +
  geom_line(aes(color = site))



#-----------------------------------------------------------------

## Make the map

library(MARSS)
library(R2jags)
library(coda)
library(rgdal)
library(rgeos)
library(sp)
library(maptools)
library(raster)
library(sf)
library(marmap)
library(grid)
library(gridBase)
library(rworldmap)
library(RColorBrewer)
library(maps)  
library(geosphere)
library(ggplot2)
library(zoo)
library(fields)
library(akima)
library(colorRamps)
library(dplyr)
library(gstat) # For semivariograms and kriging
library(spdep) # Useful package for Moran's I and similar methods
library(lattice) # Used occasionally for plotting spatial
library(ape) # What we will use for Moran's I today
library(dismo) # for Voronoi polygons




pis.sites <- read.csv(here("data/PISCO", "Swath_Transects_2003-2018.csv")) %>%
  distinct(site, LAT_WGS84, LON_WGS84, side) %>%
  group_by(site) %>% 
  filter(row_number() == 1L) %>% #THIS IS WRONG!!!!!!
  select(-side)

names(pis.sites) <- c("site", "lat", "long")
pis.sites <- pis.sites[, c("lat", "long", "site")]


sites <- read.csv(here("data/spatial", "lter_waypoints.csv")) %>%
  group_by(site) %>% 
  filter(row_number() == 1L) %>% #THIS IS WRONG!!!!!! 
  select(-c(transect, survey, lte_survey, treatment)) %>%
  bind_rows(pis.sites)


mp <- d %>% left_join(sites) %>%
  drop_na(lat, long)

mp <- as.data.frame(mp) 

coordinates(mp) <- ~long + lat

proj4string(mp) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"



cal <- readOGR(here("data/spatial", "caloutline.shp"))
all_mpas <- readOGR(here("data/spatial", "state_mpas.shp"))

d <- par(las = 1, mgp = c(3, 0.75, 0))
# blues <- colorRampPalette(c("red","purple","blue", # color ramp
#                             "cadetblue1","white"))
# plot(bath, image = TRUE, land = TRUE, lwd = 0.5, xlim = c(-120.55,-117.3), ylim = c(32.5,35), xlab = "", ylab = "", cex.axis = 1.5) # bathy layer
plot(cal, col = "#FFEB9B", xlim = c(-120.55,-118.3), ylim = c(33.8,34.8), xlab = "", ylab = "", cex.axis = 1.5, axes = TRUE) # cal state
plot(all_mpas, pch = 4, col = "#99ebff", add = T) 
plot(mp, pch = 21, cex = mp$est*500, add = T, lwd = 1.5,  bg = alpha("#e31a1c", .5))
#plot(kelp, col = "#00802b")
#plot(sites, bg = "red", add = T, pch = 21, cex = 2) 
par(d)

for(i in 2010:2018){
  
myfile <- file.path("figures/", paste("year", "_", i, ".png"))
png(myfile, width = 1000*3, height = 561*3, res = 300)
d <- par(las = 1, mgp = c(3, 0.75, 0), mar = c(4,5,3,1))
plot(cal, col = "#FFEB9B", xlim = c(-120.55,-118.3), ylim = c(33.8,34.8), xlab = "", ylab = "", cex.axis = 1.5, axes = TRUE, main = paste(i)) # cal state
plot(all_mpas, pch = 4, col = "#99ebff", add = T) 
plot(mp[mp$year == i,], pch = 21, cex = mp$est[mp$year == i]*500, add = T, lwd = 1.5,  bg = alpha("#e31a1c", .5))
par(d)
dev.off()


}
















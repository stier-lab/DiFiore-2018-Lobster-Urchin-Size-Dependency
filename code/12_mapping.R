library(here)
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
library(tidyverse)

#shoreline and MPAS
us <- readOGR(here("data/spatial", "nos80k (1)"))
all_mpas <- readOGR(here("data/spatial", "state_mpas.shp"))
all_mpas <- all_mpas[all_mpas$NAME == "Naples SMCA" | all_mpas$NAME == "Campus Point SMCA (No-Take)", ]

proj <- proj4string(all_mpas)
us <- spTransform(us, proj)
clipper <- extent(-120.197, -119.4563, 34.33406, 34.5)
shoreline <- crop(us, clipper)

# Kelp patches
patches <- readOGR(here("data/spatial", "patches.shp"))
proj4string(patches) <- "+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs"
patches <- spTransform(patches, proj)
patches <- crop(patches, clipper)

# Site data
sites <- read.csv(here("data/spatial", "lter_waypoints.csv")) %>% filter(lte_survey == "lte" & treatment == "control")
coordinates(sites) <- ~long + lat

# Depth contours
contours <- readOGR(here("data/spatial", "contours.shp"))
contours <- spTransform(contours, proj)
contours <- crop(contours, clipper)

contours <- contours[contours$ELEV %in% c(seq(-10, -100, by = -10), -200, -300, -400, -500), ]

#Predicted consumption data

dat <- read.csv(here::here("data/cleaned/posteriors/", "observational_predictions.csv")) %>% group_by(estimate) %>% filter(estimate == "mu") %>%
  group_by(site) %>%
  summarize(mean = mean(prediction)*10000)

sites <- merge(sites, dat, by = "site")
sites <- sites[order(sites$site), ]
sites$col <- c('#AF8DC3','#C3AF8D', '#c3958d', '#8DC3AF', '#c38db5')


#c(bottom, left, top, right)

# Build the map
png(here::here("figures/", "fig6_pA.png"), width = 1200, height = 375)
d <- par(las = 1, mgp = c(3, 0.75, 0), mar = c(3, 6, 1, 1), tcl = -0.5)
plot(shoreline, col = "gray80",  xlab = "", ylab = "", cex.axis = 1.5, axes = TRUE, lwd = 0.01, 
     xlim = c(-120.197, -119.4563), ylim = c(34.33406, 34.5)) # cal state
plot(all_mpas, pch = 4, col = "#99ebff", add = T, border = "#99ebff") 
plot(contours, add = T, lwd = 0.25, ext = clipper, col = "gray80")
plot(patches, col = alpha("#00802b", 0.9), add = T, border = alpha("#00802b", 0))
plot(sites, pch = 21, add = T, lwd = 1.5, cex = sites$mean*1.75, bg = alpha(sites$col, 0.85))
axis(tcl = -0.5, side = 1, labels = F)
axis(side = 2, tcl = -0.5, labels = F)
scalebar(xy = xy, type = "bar", d = 8, divs = 4, lonlat = T, below = "kilometers")
par(d)
dev.off()






















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
clipper <- extent(-120.7, -119, 34, 35)
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

source(here::here("code", "10b_scenarios.R"))

dat <- s.2a_full_wparameter %>%
  group_by(site) %>%
  mean_qi(prediction)

sites <- merge(sites, dat, by = "site")
sites <- sites[order(sites$site), ]
sites$col <- c('#AF8DC3','#C3AF8D', '#c3958d', '#8DC3AF', '#c38db5')


#c(bottom, left, top, right)
xy <- c(-119.5449, 34.4705)

# Build the map
png(here::here("figures/", "fig6_pA.png"), width = 1200, height = 375)
d <- par(las = 1, mgp = c(3, 0.75, 0), mar = c(3, 6, 1, 1), tcl = -0.5)
plot(shoreline, col = "gray80",  xlab = "", ylab = "", cex.axis = 1.5, axes = TRUE, lwd = 0.01)
     #, xlim = c(-120.197, -119.4563), ylim = c(34.33406, 34.5)) # cal state
plot(contours, add = T, lwd = 0.25, ext = clipper, col = "gray80")
#plot(patches, col = alpha("#00802b", 0.9), add = T, border = alpha("#00802b", 0))
plot(sites, pch = 21, add = T, lwd = 0.1, cex = sites$.upper*500, bg = alpha("black", 0.1))
plot(sites, pch = 21, add = T, lwd = 1.5, cex = sites$prediction*500, bg = alpha(sites$col, 0.85))
axis(tcl = -0.5, side = 1, labels = F)
axis(side = 2, tcl = -0.5, labels = F)
scalebar(xy = xy, type = "bar", d = 8, divs = 4, lonlat = T, below = "kilometers")
par(d)
dev.off()




us <- raster::getData("GADM", country = c("United States"), level = 1)
mex <- raster::getData("GADM", country = c("Mexico"), level = 1)
can <- raster::getData("GADM", country = c("Canada"), level = 1)

clipper_large <- extent(-127.5, -112, 27.5, 52.5)
shoreline_inset <- crop(us, clipper_large)
shoreline_mex <- crop(mex, clipper_large)
shoreline_can <- crop(can, clipper_large)

shore <- bind(shoreline_inset, shoreline_mex, shoreline_can)

#shoreline_can <- crop(can, clipper_large)
png(here::here("figures/", "tempmap.png"))
d <- par(bg = NA)
plot(shore, lwd = 0.01, col = "gray")
par(d)
dev.off()




# Build the map
png(here::here("figures/", "fig6_pA.png"), width = 1200, height = 375)
d <- par(las = 1, mgp = c(3, 0.75, 0), mar = c(3, 6, 1, 1), tcl = -0.5)
plot(shoreline, col = "gray80",  xlab = "", ylab = "", cex.axis = 1.5, axes = TRUE, lwd = 0.01)
# plot(contours, add = T, lwd = 0.25, ext = clipper, col = "gray80")
#plot(patches, col = alpha("#00802b", 0.9), add = T, border = alpha("#00802b", 0))
plot(sites, pch = 21, add = T, lwd = 0.1, cex = sites$.upper*250, bg = alpha("black", 0.1))
plot(sites, pch = 21, add = T, lwd = 1.5, cex = sites$prediction*250, bg = alpha(sites$col, 0.85))
axis(tcl = -0.5, side = 1, labels = F)
axis(side = 2, tcl = -0.5, labels = F)
scalebar(xy = xy, type = "bar", d = 8, divs = 4, lonlat = T, below = "kilometers")
par(d)
dev.off()


temp <- st_as_sf(sites, crs = 4326) %>% st_set_crs(4326)
temp2 <- st_as_sf(shoreline, crs = 4326) %>% st_set_crs(4326)
temp2 <- st_set_crs(temp2, crs = 4326)
ggplot()+
  geom_sf(data = temp2, fill = "white")+
  geom_sf(data = temp, aes(size = prediction))


ggplot()+
  geom_sf(data = temp2, fill = "white")



library(ggspatial)


source(here::here("code/", "10b_scenarios.R"))

#shoreline and MPAS
library(raster)
library(sf)
library(rgdal)
library(here)
library(ggplot2)
library(ggspatial)
library(tidyverse)

clipper <- st_polygon(list(rbind(c(-120.7, 34.65),
                                 c(-119.25, 34.65), 
                                 c(-119.25, 33.8), 
                                 c(-120.7, 33.8), 
                                 c(-120.7, 34.65)))) %>% st_sfc() %>% st_set_crs(4326)

clipper <- st_polygon(list(rbind(c(-121.75, 34.65),
                                 c(-119.25, 34.65), 
                                 c(-119.25, 33.8), 
                                 c(-121.75, 33.8), 
                                 c(-121.75, 34.65)))) %>% st_sfc() %>% st_set_crs(4326)

  

us <- readOGR(here("data/spatial", "nos80k (1)")) %>% st_as_sf() %>% st_transform(4326) %>% st_intersection(clipper)

patches <- readOGR(here("data/spatial", "patches.shp")) %>% st_as_sf() %>% st_set_crs("+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs") %>% st_transform(4326) %>% st_buffer(dist = 0) %>% st_intersection(clipper)

# contours <- readOGR(here("data/spatial", "contours.shp")) %>% st_as_sf() %>% st_transform(4326) %>% filter(ELEV %in% c(seq(-10, -100, by = -10), -200, -300, -400, -500))

dat <- s.2a_full_wparameter %>%
  group_by(site) %>%
  mean_qi(prediction)

sites <- read.csv(here("data/spatial", "lter_waypoints.csv")) %>% filter(lte_survey == "lte" & treatment == "control") %>% st_as_sf(coords = c("long", "lat"), crs = 4326, agr = "constant") %>% left_join(dat)

depth <- marmap::read.bathy("~/Downloads/sbc.xyz")
coords <- st_coordinates(clipper)
depth <- marmap::subsetBathy(depth, x = coords[,1], y = coords[,2], locator = F)
depth.df <- marmap::fortify.bathy(depth) #%>% filter(z <= 0)


blues <- colorRampPalette(colors = c("#7FA8CB", "#DAF0FD" ))
browns <- colorRampPalette(colors = c("#ACD0A5", "#C3A76B"))

breaks <- round(c(seq(min(depth.df$z), 0, length.out = 200), seq(1, max(depth.df$z), length.out = 200)))


zoom_map <- ggplot()+
  geom_tile(data = filter(depth.df, z < 0), aes(x = x, y = y, fill = z))+
  scale_fill_gradientn(colours = blues(10))+
  geom_contour(data = filter(depth.df, z < 10), aes(x = x, y = y, z = z), color = "black", binwidth = 100, alpha = 0.25)+
  geom_sf(data = us, fill = "gray", lwd = 0.01)+
  geom_sf(data = patches, fill = alpha("#00802b", 0.9), col = alpha("#00802b", 0))+
  geom_sf(data = sites, aes(size = .upper), alpha = 0.1)+
  geom_sf(data = sites, aes(size = prediction, color = site))+
  scale_size(range = c(2,22))+
  scale_color_manual(values = c('#AF8DC3','#C3AF8D', '#c3958d', '#8DC3AF', '#c38db5'))+
  labs(x = "", y = "", size = "")+
  guides(color = FALSE, fill = FALSE, size = guide_legend(direction = "horizontal", label.position = "top", label.vjust = 0, override.aes = list(shape = 21)))+
  theme(legend.position = c(0.3, 0.85))
  

ggsave(filename = here::here("figures/ggmap.png"),plot = zoom_map, device = "png" )


clipper_large <- st_polygon(list(rbind(c(-165, 62),
                      c(-165, 21), 
                      c(-105, 21), 
                      c(-105, 62), 
                      c(-165, 62)))) %>% st_sfc() %>% st_set_crs(4326)

world1 <- sf::st_as_sf(maps::map('world', plot = FALSE, fill = TRUE)) %>% st_transform(4326) %>% st_intersection(clipper_large) %>% st_union()

map <- ggplot() + 
  geom_sf(data = world1) + 
  geom_sf(data = clipper, color = "red", fill = NA, lwd = 1)+
  theme_classic()

ggsave(filename = here::here("figures/ggmap_big.png"),plot = map, device = "png" )


panel_map <- cowplot::plot_grid(map, zoom_map, align = "hv", nrow = 1)

ggsave(here::here("figures/panel_map.pdf"), plot = panel_map, device = "pdf", width = 15, height = 7.5)



map <- ggplot() + 
  geom_sf(data = world1, fill = "white") + 
  geom_sf(data = clipper, color = "red", fill = NA, lwd = 1)+
  theme_void()


gg_inset_map2 = ggdraw() +
  draw_plot(zoom_map) +
  draw_plot(map, x = 0.25, y = 0.65, width = 0.25, height = 0.25)

ggsave(here::here("figures/panel_map.pdf"), plot = gg_inset_map2, device = "pdf")


















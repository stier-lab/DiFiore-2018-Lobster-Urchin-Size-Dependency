source(here::here("code/", "10b_scenarios.R"))

#shoreline and MPAS
library(raster)
library(sf)
library(rgdal)
library(here)
library(ggplot2)
library(ggspatial)
library(tidyverse)
library(ggsn)



#us <- readOGR(here("data/spatial", "nos80k (1)")) %>% st_as_sf() %>% st_transform(4326)

us <- raster::getData("GADM", country = c("United States"), level = 1)


clipper_small <- st_polygon(list(rbind(c(-120.7, 34.65),
                                 c(-119.25, 34.65), 
                                 c(-119.25, 33.8), 
                                 c(-120.7, 33.8), 
                                 c(-120.7, 34.65)))) %>% st_sfc() %>% st_set_crs(4326)

coord_sf(xlim = c(-120.7, -119.25), ylim = c(33.8, 34.65) )+

shore_small <- us %>% sf::st_as_sf() %>% sf::st_intersection(clipper_small) %>% sf::st_transform(4326) %>% sf::st_union()


patches <- readOGR(here("data/spatial", "patches.shp")) %>% st_as_sf() %>% st_set_crs("+proj=utm +zone=10 +datum=WGS84 +units=m +no_defs") %>% st_transform(4326) %>% st_buffer(dist = 0) %>% sf::st_intersection(clipper_small)

dat <- s.2a_full_wparameter %>%
  group_by(site) %>%
  mean_qi(prediction)

sites <- read.csv(here("data/spatial", "lter_waypoints.csv")) %>% filter(lte_survey == "lte" & treatment == "control") %>% st_as_sf(coords = c("long", "lat"), crs = 4326, agr = "constant") %>% left_join(dat)

depth <- marmap::read.bathy("~/Downloads/sbc.xyz")
coords <- st_coordinates(clipper_small)
depth <- marmap::subsetBathy(depth, x = coords[,1], y = coords[,2], locator = F)
depth.df <- marmap::fortify.bathy(depth)


blues <- colorRampPalette(colors = c("#94AAC3", "#F9FAFB")) #Low #94AAC3, high #F9FAFB
browns <- colorRampPalette(colors = c("#ACD0A5", "#C3A76B"))

# breaks <- round(c(seq(min(depth.df$z), 0, length.out = 200), seq(1, max(depth.df$z), length.out = 200)))


zoom_map <- ggplot()+
  geom_tile(data = filter(depth.df, z < 0), aes(x = x, y = y, fill = z))+
  scale_fill_gradientn(colours = blues(10))+
  geom_contour(data = filter(depth.df, z < 10), aes(x = x, y = y, z = z), color = "black", binwidth = 100, alpha = 0.25)+
  geom_sf(data = shore_small, fill = "#596778", lwd = 0.01)+
  geom_sf(data = patches, fill = alpha("#00802b", 0.9), col = alpha("#00802b", 0))+
  geom_sf(data = sites, aes(size = .upper), alpha = 0.1)+
  geom_sf(data = sites, aes(size = prediction, color = site, alpha = 0.95), show.legend = F)+
  scale_size(range = c(2,22))+
  scale_color_manual(values = c('#AF8DC3','#C3AF8D', '#c3958d', '#8DC3AF', '#c38db5'))+
  labs(x = "", y = "", size = expression(paste("Interaction\nstrength (ind. m"^-2,"d"^-1,")")))+
  coord_sf(xlim = c(-120.7, -119.25), ylim = c(33.8, 34.65), expand = F)+
  scale_x_continuous(breaks = -1*seq(119.5, 120.5, by = 0.5))+
  annotation_scale(location = "bl", style = "ticks",  pad_x = unit(0.25, "cm"), pad_y = unit(0.25, "cm"))+
  annotation_north_arrow(location = "tr", style = north_arrow_nautical, height = unit(0.75, "cm"), width = unit(0.75, "cm"), pad_x = unit(0.25, "cm"), pad_y = unit(0.25, "cm"))+
  guides(color = FALSE, fill = FALSE, size = guide_legend(override.aes = list(shape = 21, color = "black", alpha = 1)))+
  theme(legend.title=element_text(size=14))
  #guides(color = FALSE, fill = FALSE, size = guide_legend(direction = "horizontal", label.position = "top", label.vjust = 0, override.aes = list(shape = 21)))+
  #theme(legend.position = c(0.4, 0.9), panel.border = element_rect(colour = "black", fill=NA, size=1))

ggsave(filename = here::here("figures/ggmap.png"),plot = zoom_map, device = "png", width = 12*0.85, height = 8*0.85)




clipper_large <- st_polygon(list(rbind(c(-135, 21),
                                       c(-135, 62), 
                                       c(-105, 62), 
                                       c(-105, 21), 
                                       c(-135, 21)))) %>% st_sfc() %>% st_set_crs(4326)

world1 <- sf::st_as_sf(maps::map('world', region = c("USA", "Canada", "Mexico"), plot = FALSE, fill = TRUE)) %>% st_transform(4326)  %>% st_union() %>% st_intersection(clipper_large)


map <- ggplot() + 
  geom_sf(data = world1, fill = "#C2BDBD", col = "#C2BDBD") + 
  geom_sf(data = clipper_small, color = "red", fill = NA, lwd = 1)+
  coord_sf(xlim = c(-135, -105), ylim = c(22,61), expand = F)+
  scale_x_continuous(breaks = c(-110, -120, -130))+
  scale_y_continuous(breaks = c(seq(25, 55, by = 15)))+
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
  theme(axis.text = element_text(size = 10))


ggsave(filename = here::here("figures/ggmap_big.png"),plot = map, device = "png", width = 12*0.25, height = 8*0.85 )



panel_map <- cowplot::plot_grid(map, zoom_map, nrow = 1, rel_widths = c(0.325, 1), align = "v")

ggsave(here::here("figures/panel_map.png"), plot = panel_map, device = "png")
                                


                                
                                
                                
#-------------------
## Scrap
#-------------------


panel_map <- cowplot::plot_grid(map, zoom_map, nrow = 1, rel_widths = c(1, 2), align = "h")

save_plot(here::here("figures/panel_map.png"), plot = panel_map, device = "png")





gg_inset_map2 = ggdraw() +
  draw_plot(zoom_map) +
  draw_plot(map, x = .75, y = 0.75, width = 0.25, height = 0.25)

ggsave(here::here("figures/panel_map.pdf"), plot = gg_inset_map2, device = "pdf")
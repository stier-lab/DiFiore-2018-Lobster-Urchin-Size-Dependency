



forage <- read.csv("~/Downloads/resource_map_doi_10_5063_F17H1GTQ/data/FoRAGE_db_12_19_18_data_set.csv") %>% janitor::clean_names() %>%
  mutate(maxC = 1/fittted_h_day) %>%
  drop_na(predator_mass_mg, prey_mass_mg, temp_c)

lm1 <- lm(log(maxC) ~ log(predator_mass_mg) + log(prey_mass_mg) + temp_c + I(temp_c^2), forage)

car::crPlot(model = lm1, variable = "log(predator_mass_mg)")

lm1.1 <- lm(log(maxC) ~ log(predator_mass_mg) + log(prey_mass_mg) + temp_c + I(temp_c^2), forage)
lm1.2 <- lm(log(maxC) ~ log(prey_mass_mg) + temp_c + I(temp_c^2), forage)

plot(residuals(lm1.2) ~ log(predator_mass_mg) , forage)
summary(lm1)

forage$reslm1.1 <- residuals(lm1.1)
forage$reslm1.2 <- residuals(lm1.2)


col <- colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))

p1 <- forage %>% filter(exp(reslm1.2) > 10e-6, major_grouping_1 != "Rotifer") %>%
  ggplot(aes(y = reslm1.2, x = predator_mass_mg))+
  geom_point(color = "lightgray", alpha = 0.75, size = 1.5)+
  scale_x_log10(breaks = c(0.0001, 1, 100, 1000, 10000), labels = c(0.0001, 1, 100, 1000, 10000))+
  coord_cartesian(ylim = c(-10, 25))+
  labs(x = "Predator size (mg)", y = "Consumption rate")+
  theme_classic()+
  theme(text = element_text(size = 24))+
  theme(legend.position = c(0.25, 0.75), text = element_text(size = 18, color = "lightgray"), axis.line = element_line(color = "lightgray"), axis.ticks = element_line(color = "lightgray"))+
  theme(
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill='transparent'),
    legend.box.background = element_rect(fill='transparent')
  )

ggsave(filename = "figures/pres-p1.png", plot = p1, device = "png", bg = "transparent")


p2 <- forage %>% filter(exp(reslm1.2) > 10e-6, major_grouping_1 != "Rotifer") %>%
  ggplot(aes(y = reslm1.2, x = predator_mass_mg))+
  geom_point(color = "lightgray", alpha = 0.75, size = 1.5)+
  geom_smooth(method = "lm", col = "#ee9a00", lwd = 1.5)+
  scale_x_log10(breaks = c(0.0001, 1, 100, 10000), labels = c(0.0001, 1, 100, 10000))+
  coord_cartesian(ylim = c(-10, 25))+
  labs(x = "Predator size (mg)", y = "Consumption rate")+
  theme_bd()+
  theme(legend.position = c(0.25, 0.75), text = element_text(size = 18, color = "lightgray"), axis.text = element_text(size = 14, color = "lightgray"), axis.line = element_line(color = "lightgray"), axis.ticks = element_line(color = "lightgray"))+
  theme(
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill='transparent'),
    legend.box.background = element_rect(fill='transparent')
  )

ggsave(filename = "figures/pres-p2.png", plot = p2, device = "png", bg = "transparent", width = 6, height = 10/3)

p3 <- forage %>% filter(exp(reslm1.2) > 10e-6, major_grouping_1 != "Rotifer") %>%
  ggplot(aes(y = reslm1.2, x = predator_mass_mg))+
  geom_point(aes(color = major_grouping_1), show.legend = F, alpha = 0.75, size = 1.5)+
  scale_color_manual(values = col(19))+
  geom_smooth(method = "lm", aes(col = major_grouping_1), alpha = 0.25, show.legend = F)+
  geom_smooth(method = "lm", col = "darkred", se = F, lwd = 1.5)+
  scale_x_log10(breaks = c(0.0001, 1, 100, 1000, 10000), labels = c(0.0001, 1, 100, 1000, 10000))+
  coord_cartesian(ylim = c(-10, 25))+
  labs(x = "Predator size (mg)", y = "Consumption rate")+
  theme_classic()+
  theme(text = element_text(size = 30), 
        axis.text = element_text(size = 24))+
  theme(legend.position = c(0.25, 0.75), text = element_text(size = 18, color = "lightgray"), axis.line = element_line(color = "lightgray"), axis.ticks = element_line(color = "lightgray"))+
  theme(
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill='transparent'),
    legend.box.background = element_rect(fill='transparent')
  )

ggsave(filename = "figures/pres-p3.png", plot = p3, device = "png")


forgg <- forage %>% filter(exp(reslm1.2) > 10e-6, major_grouping_1 != "Rotifer")

bo <- data.frame(mc = forgg$predator_mass_mg, mr = forgg$prey_mass_mg)
bo$MaxC <- -5.04 + 1.44*log(bo$mc) + 0.27*log(bo$mr)
bo$MaxC <- exp(bo$MaxC)

bo1.2 <- lm(log(MaxC) ~ log(mr), bo)
bo$bo1.2 <- residuals(bo1.2)

rall <- data.frame(mc = forgg$predator_mass_mg, mr = forgg$prey_mass_mg)
rall$log.h <- 10.38 + -0.76*log(rall$mc) + 0.76*log(rall$mr)
rall$h <- exp(rall$log.h)*60*60*24
rall$MaxC <- 1/rall$h

rall1.2 <- lm(log(MaxC) ~ log(mr), rall)
rall$rall1.2 <- residuals(rall1.2)
plot(log(MaxC) ~ log(mc), rall)




#9f5e85
p4 <- forage %>% filter(exp(reslm1.2) > 10e-6, major_grouping_1 != "Rotifer") %>%
  ggplot(aes(y = reslm1.2, x = predator_mass_mg))+
  geom_point(color = "lightgray", alpha = 0.75, size = 1.5)+
  scale_color_manual(values = col(19))+
  geom_smooth(method = "lm", alpha = 0.25, show.legend = F, col = "#ee9a00", lwd = 1.5)+
  geom_smooth(data = bo, aes(y = bo1.2, x = mc), method  = "lm", col = "#859f5e", lwd = 1.5)+
  geom_smooth(data = rall, aes(y = rall1.2, x = mc), method = "lm", col = "#5e859f", lwd = 1.5)+
  scale_x_log10(breaks = c(0.0001, 1, 100, 10000), labels = c(0.0001, 1, 100, 10000))+
  coord_cartesian(ylim = c(-10, 25))+
  labs(x = "Predator size (mg)", y = "Consumption rate")+
  theme_bd()+
  theme(legend.position = c(0.25, 0.75), text = element_text(size = 18, color = "lightgray"), axis.text = element_text(size = 14, color = "lightgray"), axis.line = element_line(color = "lightgray"), axis.ticks = element_line(color = "lightgray"))+
  theme(
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill='transparent'),
    legend.box.background = element_rect(fill='transparent')
  )

ggsave(filename = "figures/pres-p4.png", plot = p4, device = "png", bg = "transparent", width = 6, height = 10/3)



























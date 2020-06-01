################################################
## GLMM 
################################################

library("here")
source(here("code","setup.R"))
source(here("code","functions.R"))
source(here("code","clean.R"))


df$urc.density <- df$num_offered/2
res <- cbind(df$num_consumed, df$num_offered)

glm1 <- glm(res ~ urc.density + size + treatment, data = df, family = "binomial")
summary(glm1)


glmer1 <- glmer(res ~ scale(urc.density) + scale(size) + treatment + (1|id), data = df, family = "binomial")
summary(glmer1)
hist(residuals(glmer1), breaks = 30)

glmer2 <- glmer(res ~ scale(urc.density) + scale(size) * treatment + (1|id), data = df, family = "binomial", control=glmerControl(optCtrl=list(maxfun=50000)))
summary(glmer2)



newdata <- with(df, expand.grid(urc.density = seq(1,13.5, by = 0.1), size = unique(size), treatment=unique(treatment)))
newdata$id <- rep(seq(1,43), each = 378)
newdata$predicted <- predict(glmer1, newdata, type = "response", allow.new.levels = T, re.form = NA)

#c(bottom, left, top, right)
pdf("figures/fig_glmer.pdf", width = 13.3*0.8, height = 7.5*0.8)
d <- par(mfrow = c(1,2), mar = c(4, 4, 1, 0.5), las = 1)
plot(I(num_consumed/num_offered) ~ jitter(urc.density), data = df, xlab = "Urchin density (ind./m2)", ylab = "Proportional mortality")
lines(predicted ~ urc.density, data = newdata[newdata$treatment == "urc_small" & newdata$size == 89.3, ], col = "#bdc9e1", lwd = 1.2)
lines(predicted ~ urc.density, data = newdata[newdata$treatment == "urc_medium" & newdata$size == 89.3, ], col = "#74a9cf", lwd = 2)
lines(predicted ~ urc.density, data = newdata[newdata$treatment == "urc_large" & newdata$size == 89.3, ], col = "#0570b0", lwd = 3.5)
LEGEND(x = 8, y = 1, legend = c("Small urchins", "Medium urchins", "Large urchins"), lty = c(1,1,1), lwd = c(1.2,2,3.5), col = c("#bdc9e1", "#74a9cf", "#0570b0"), xjust = 0, yjust = 1, y.intersp = 1, x.intersp = 0.75, seg.len = 2, bg = alpha("white", .8))



l1 <- newdata[newdata$treatment == "urc_small" & newdata$urc.density == 13.5, ]
l1 <- l1[order(l1$size), ]

l2 <- newdata[newdata$treatment == "urc_medium" & newdata$urc.density == 13.5, ]
l2 <- l2[order(l2$size), ]

l3 <- newdata[newdata$treatment == "urc_large" & newdata$urc.density == 13.5, ]
l3 <- l3[order(l3$size), ]

plot(I(num_consumed/num_offered) ~ jitter(size), data = df, xlab = "Lobster size (mm)", ylab = "")
lines(predicted ~ size, data = l1, col = "#bdc9e1", lwd = 1.2)
lines(predicted ~ size, data = l2, col = "#74a9cf", lwd = 2)
lines(predicted ~ size, data = l3, col = "#0570b0", lwd = 3.5)

par(d)
dev.off()




# library(sjPlot)
# library(sjmisc)
# library(sjlabelled)
# library(sjstats)
# 
# labels <- c(`(Intercept` = "Large sized urchins", `scale(urc.density)` = "Prey density", `scale(size)` = "Predator size", `treatmenturc_medium` = "Medium sized urchins",`tratmenturc_small` = "Small sized urchins")
# 
# tab_model(glmer1, show.std = TRUE, show.se = TRUE, transform = NULL, show.intercept = T, show.ci = F,  pred.labels = labels, show.icc = F, show.r2 = F, file = "figures/table_glmer.html")

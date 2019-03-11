##########################################################
## Plot Jags output
##########################################################

# Write function to plot fits
plot.jags <- function(model, pred_size, urc_size){
  temp <- df %>% filter(class == pred_size, treatment == urc_size)
  killed <- temp$num_consumed
  initial <- temp$num_offered
  plot(killed ~ jitter(initial),data = temp, xlab="Number of prey",ylab="Number consumed", ylim = c(0,26), main = paste(pred_size, "lobster", "--", urc_size, "urchins"))
  a = MCMCsummary(model,params='a')[1] # the alpha estimate here is often bounding up against zero
  h = MCMCsummary(model,params='h')[1]
  curve(holling2(x,a,h,P=1,T=1),add=TRUE,col=1,lty=1) #true curve
  
}

#############################################################
## Plot 12 panel of fits
#############################################################

lob.size <- c("small", "medium", "large", "jumbo")
urc.size <- c("urc_small", "urc_medium", "urc_large")
treat2 <- data.frame(model = as.character(mods), lob.size = rep(lob.size, times = 3), urc.size = rep(urc.size, each = 4))

d <- par(mfrow = c(3,4))
for(i in 1:nrow(treat2)){
  temp <- treat2$model[i]
  plot.jags(get(substitute(temp)), treat2$lob.size[i], treat2$urc.size[i])
}
par(d)

plot.jags(mod1, small, small)

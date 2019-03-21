dunn <- read.csv("data/raw/Dunnetal2019_bootstrap.csv") # get the data

#make a histogram 
hist(dunn[, "a"])
hist(dunn[, "h"])

# calculate central tendencies
mean.a <- mean(dunn[,"a"])
median.a <- median(dunn[,"a"])
var.a <- var(dunn[,"a"])

mean.h <- mean(dunn[, "h"])
median.h <- median(dunn[,"h"])
var.h <- var(dunn[,"h"])

sum <- list(mean.a, median.a, var.a, mean.h, median.h, var.h)
names(sum) <- c("mean.a", "median.a", "var.a", "mean.h", "median.h", "var.h")


#from Bolker et al. 2008  scale = v/m,  shape = m2/v; where m is the mean, and v is the variance. 

scale.a <- var.a / mean.a
scale5.a <- (var.a*5)/mean.a # increase the variance to see if we can "flatten" the prior...
shape.a <- mean.a^2 / var.a
shape5.a <- mean.a^2/ (var.a*5) # increase the variance to see if we can "flatten" the prior...

scale.h <- var.h / mean.h
shape.h <- mean.h^2 / var.h

x <- seq(0, 1, by = 0.01)
x2 <- seq(0,2.5, by = 0.01)
# plot to inspect distributions

plot(dgamma(x = x, shape = shape.a, scale = scale.a ) ~ x)
points(dgamma(x = x, shape = shape5.a, scale = scale5.a ) ~ x, col = 2) # changing the variance also shift the mean of the distribution... so likely not a good strategy
plot(dgamma(x = x2, shape = shape.h, scale = scale.h ) ~ x2)



















##############################################################################
## Misc from http://doingbayesiandataanalysis.blogspot.com/2012/08/gamma-likelihood-parameterized-by-mode.html
##############################################################################

# 
# #To parameterize the gamma prior by the mean
# model {
#   for ( i in 1:N ) {
#     y[i] ~ dgamma( sh , ra )
#   }
#   # parameterized by mean (m) and standard deviation (sd)
#   sh <- pow(m,2) / pow(sd,2)
#   ra <-     m    / pow(sd,2)
#   m ~ dunif(0,100)
#   sd ~ dunif(0,100)
# }
# 
# #To parameterize the gamma prior by the mode (supposedly better because of skewedness)
# model {
#   for ( i in 1:N ) {
#     y[i] ~ dgamma( sh , ra )
#   }
#   # parameterized by mode (m) and standard deviation (sd):
#   sh <- 1 + m * ra
#   ra <- ( m + sqrt( m^2 + 4*sd^2 ) ) / ( 2 * sd^2 )
#   m ~ dunif(0,100)
#   sd ~ dunif(0,100)
# }



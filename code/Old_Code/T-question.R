library(here)
source(here("code","functions.R"))
source(here("code","setup.R"))

sam <- read.csv(here("data/samdata", "fr_data.csv"))
sam$id <- as.numeric(sam$lobster_id)

sam$temp2 <- as.factor(paste("t", sam$temp, sep = ""))

s <- sam[,c( "temp2", "lobster_id","Initial", "Killed")]
names(s) <- c("temp", "id", "initial", "killed")

s <- arrange(s, temp, id, initial)

n <- filter(s, id == "N14")


plot(killed ~ jitter(initial), n)


library(frair)

m1 <- nls(
  killed ~ (a*initial)/(1+a*initial*h),
  start = c(a = 0.01, h = 0.001), data = n)
summary(m1)


m2 <- nls(
  killed ~ (a*48*initial)/(1+a*initial*h),
  start = c(a = 0.01, h = 0.001), data = n)
summary(m2)

d <- par(mfrow = c(1,2), mar = c(4,4,1,1))
plot(killed ~ jitter(initial),data = n, xlab="Number of prey",ylab="Number consumed")
curve(holling2(x,a = coef(m1)[1], coef(m1)[2],P=1,T=1),add=TRUE,col = "red",lty=1) #true curve
plot(I(killed/48) ~ jitter(I(initial)),data = n, xlab="Number of prey",ylab="Number consumed per hour")
curve(holling2(x,a = coef(m2)[1], coef(m2)[2],P=1,T=1),add=TRUE,col="green",lty=2) #true curve
par(d)

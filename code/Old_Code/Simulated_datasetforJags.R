mortrisk <- function(N0,h,a){
  risk <- a/(1+a*N0*h)                 ## Holling risk (per capita)
  pmin(1.0,risk)                         ## bound risk <= 1
}






set.seed(1001) ## set random-number seed for reproducibility
## The data will simulate a case where initial attack rate "a" varies
## across 6 replicate blocks
simdata <- function(nind){
  test.vals <- expand.grid(N0=c(2,3,5,10,16,26),
                           ind=1:nind)
  ## attack rate varies randomly by individual with a median of 0.5
  ##  and proportional variation of approx 10%
  a <- rlnorm(nind,meanlog=log(0.5),sdlog=0.1)
  ## handling time varies randomly by individual with a median of 0.1
  ##  and proportional variation of approx 10%
  h <- rlnorm(nind,meanlog=log(0.1),sdlog=0.1)
  p <- with(test.vals,mortrisk(N0=N0,
                               a=a[ind],
                               h=h[ind]))
  z <- rbinom(nrow(test.vals),prob=p,size=test.vals$N0)
  data.frame(test.vals,killed=z)
}
x <- simdata(46)
## Plot results ...
with(x,plot(jitter(N0),killed,col=ind))
## Now fit the data using maximum likelihood without block effects
eps <- 1e-4 ## used to bound probabilities between 0 and 1
##Using the classical unstructured Type II functional response
classic <- mle2(killed~dbinom(prob=pmax(eps,
                                        pmin(1-eps,1/(1/a+h*N0))),
                              size=N0),start=list(a=.01,h=.001),data=x)


names(x) <- c("initial", "ind", "killed")

## linear shift between two samples
ssize <- 100
beta <- 1.5
set.seed(123)
sample1 <- runif(ssize,1,ssize)
sample2 <- beta * sample1 + rnorm(ssize,sd=ssize/10)
plot(sample1,sample2)
abline(0,1,col="red")

## quadratic shift between two samples
set.seed(456)
sample3 <- beta/6 * sample1 + beta/100 * sample1^2 + rnorm(ssize,sd=ssize/10)
plot(sample1,sample3)
abline(0,1,col="red")

## check the fit
sample12 <- sample1^2
LM21 <- lm(sample2 ~ sample1)
LM31 <- lm(sample3 ~ sample1 + sample12)

## create the dataset composed of the three samples
dat <- data.frame(sample1,sample2,sample3,row.names=paste("gene",seq(1,ssize),sep=""))
saveRDS(dat,file=file.path(Sys.getenv("OMPATH"),"data/toy3.RDS"))

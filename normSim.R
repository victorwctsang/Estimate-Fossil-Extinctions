## Testing regression-based inversion code (regInversion.R)
## Test case is estimating the confidence level `q` for the mean
## from a normal random sample with known sd (equal to one)
## using sample mean as the test stat

simFn = function(mu,n,sigma=1){ rnorm(n,mu,sigma) }

q=0.025
n=10
dat=simFn(0,10)

thetas=ft$theta
errors=ft$error


ft_rq = regInversion(data=dat,getT=mean,simulateData=simFn, thetaInits=c(-1,1),n=n,iterMax=400,q=q)
ft_prob = regInversion(data=dat,getT=mean,simulateData=simFn, thetaInits=c(-1,1),n=n,iterMax=400,method="prob",q=q)

thetas = matrix(c(ft_rq$theta,ft_prob$theta),ncol=2)
errors = matrix(c(ft_rq$error,ft_prob$error),ncol=2)
colnames(thetas) = colnames(errors) = c("quantile","glm")

nTimes = 50
for(iTime in 1:nTimes)
{
  ft_rq = regInversion(data=dat,getT=mean,simulateData=simFn, thetaInits=c(-1,1),n=n,iterMax=400,q=q)
  ft_prob = regInversion(data=dat,getT=mean,simulateData=simFn, thetaInits=c(-1,1),n=n,iterMax=400,method="prob",q=q)

  thetas = rbind(thetas,c(ft_rq$theta,ft_prob$theta))
  errors = rbind(errors,c(ft_rq$error,ft_prob$error))
  print(paste0("Completed ",iTime, " of ", nTimes, " runs..."))
}

n = length(dat)
thetaHat = mean(dat)+qnorm(q)*1/sqrt(n)

boxplot(thetas,ylim=thetaHat+c(-0.5,0.5)) # values should be in this range...
abline(h=thetaHat,col="red")
MSEs = apply((thetas-thetaHat)^2,2,mean)*1000
BIASes = apply(thetas-thetaHat,2,mean)
cat("\n MSE \n")
print(MSEs)
cat("\n Biases \n")
print(BIASes)
cat("\n empirical SE of theta \n")
apply(thetas,2,sd)
cat("\n Estimated SE from regression model \n")
print(apply(errors,2,mean))

cat("\n Reporting results for 'clean' cases \n")
oops = apply(abs(thetas-thetaHat),1,sum)>0.5
cleanThetas = thetas[oops==FALSE,]
cleanMSEs = apply((cleanThetas-thetaHat)^2,2,mean)*1000
cleanBIASes = apply(cleanThetas-thetaHat,2,mean)
cat("\n MSE \n")
print(cleanMSEs)
cat("\n Biases \n")
print(cleanBIASes)
cat("\n empirical SE of theta \n")
apply(cleanThetas,2,sd)
cat("\n Estimated SE from regression model \n")
print(apply(errors[oops==FALSE,],2,mean))

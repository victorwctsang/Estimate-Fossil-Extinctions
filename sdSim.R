## Testing regression-based inversion code (regInversion.R) for non-linear case
## Test case is estimating the confidence level `q` for the log(sd)
## from a normal random sample with mean zero, using sample sd as the test stat

simFn = function(lsigma,n,mu=0){ rnorm(n,mu,exp(lsigma)) }

q=0.025
n=10
dat=simFn(0,10)


sumSq = function(dat,n=NULL) {sum(dat^2)}
thetaHat=log(sqrt(sumSq(dat)/qchisq(1-q,n)))
tObs=sumSq(dat)


source("regInversion.R")
ft_rq = regInversion(data=dat, getT=sumSq, simulateData=simFn, thetaInits=seq(-1,1,length=20),n=n,iterMax=1000,q=q)
ft_wrq = regInversion(data=dat, getT=sumSq, simulateData=simFn, thetaInits=seq(-1,1,length=20),n=n,iterMax=1000,method="wrq",q=q)
ft_rq2 = regInversion(data=dat, getT=sumSq, simulateData=simFn, thetaInits=seq(-1,1,length=20),n=n,iterMax=1000,method="rq2",q=q)
str(ft_wrq)

#ft_qgam = regInversion(data=dat,getT=sumSq,simulateData=simFn, thetaInits=seq(-1,1,length=20),n=n,iterMax=400,method="qgam",q=q)

par(mfrow=c(2,1))
plot(T~theta,data=ft_rq$stats)
plot(T~theta,data=ft_wrq$stats)
ftq=quantreg::rq(T~theta,tau=1-q,data=ft_rq$stats)

ft_rq2 = regInversion(data=dat, getT=sumSq, simulateData=simFn, thetaInits=seq(-1,1,length=20),n=n,iterMax=1000,method="rq2",q=q)
plot(T~theta,data=ft_rq2$stats,xlim=c(-3,1),ylim=c(0,20))
#OKstats = ft_rq2$stats$T<5*tObs & ft_rq2$stats$T>tObs/5
#ftq=quantreg::rq(T~poly(theta,2),tau=1-q,data=ft_rq2$stats[OKstats,])
ftq=quantreg::rq(T~poly(theta,2),tau=1-q,data=ft_rq2$stats[OKstats,])
ts=seq(-4,2,length=1000)
lines(ts,predict(ftq,newdata=list(theta=ts)),col="blue")
abline(h=tObs,col="red")
abline(v=thetaHat,col="red")


plot(T~theta,data=ft_qgam$stats)

ftq=quantreg::rq(T~poly(theta,2),tau=q,data=ft_rq2$stats)
ts=seq(-4,0,length=100)
lines(ts,predict(ftq,newdata=list(theta=ts)),col="blue")
abline(h=tObs,col="red")
abline(v=thetaHat,col="red")

str(ft_qgam)

th=ft_rq2$stats$theta
thPoly = poly(ft_rq2$stats$theta,2)
thetaMat = cbind(1,thPoly)
thVec = seq(min(ft_rq2$stats$theta),max(ft_rq2$stats$theta),length=100)
thPr = predict(thPoly,newdata=thVec)
thetaPr = cbind(1,thPr)
ses=diag(thetaPr%*% solve( t(thetaMat)%*% thetaMat  ) %*%t(thetaPr) )
plot(thVec,ses)
abline(v=mean(ft_rq2$stats$theta))
abline(v=median(ft_rq2$stats$theta,col="blue"))

hist(ft_rq2$stats$theta,20)


ft2=lm(th^2~th)
plot(residuals(ft2),thPoly[,2])

thetas = matrix(c(ft_rq$theta,ft_wrq$theta,ft_rq2$theta),ncol=3)
errors = matrix(c(ft_rq$error,ft_wrq$error,ft_rq2$error),ncol=3)
colnames(thetas) = colnames(errors) = c("rq","wrq","rq2")

nTimes = 50
for(iTime in 1:nTimes)
{
  ft_rqi = regInversion(data=dat,getT=sumSq,simulateData=simFn, thetaInits=seq(-1,1,length=20),n=n,iterMax=400,method="rq",q=q)
  ft_wrqi = regInversion(data=dat,getT=sumSq,simulateData=simFn, thetaInits=seq(-1,1,length=20),n=n,iterMax=400,method="wrq",q=q)
  ft_rq2i = regInversion(data=dat,getT=sumSq,simulateData=simFn, thetaInits=seq(-1,1,length=20),n=n,iterMax=400,method="rq2",q=q)
  
  thetas = rbind(thetas,c(ft_rqi$theta,ft_wrqi$theta,ft_rq2i$theta))
  errors = rbind(errors,c(ft_rqi$error,ft_wrqi$error,ft_rq2i$error))
  print(paste0("Completed ",iTime, " of ", nTimes, " runs..."))
}

n = length(dat)

boxplot(thetas,ylim=thetaHat+c(-0.25,0.25)) # values should be in this range...
abline(h=thetaHat,col="red")
MSEs = apply((thetas-thetaHat)^2,2,mean)*100
BIASes = apply(thetas-thetaHat,2,mean)*10
cat("\n MSE \n")
print(MSEs)
cat("\n Biases \n")
print(BIASes)
cat("\n empirical SE of theta \n")
apply(thetas,2,sd)*10
cat("\n Estimated SE from regression model \n")
print(apply(errors,2,mean)*10)
title("a=0 noooooutliers")

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

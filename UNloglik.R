# UNloglik finds the log-likelihood of our uniform-normal model
# given a set of ages and their sds

UNloglik = function(theta,ages,sds,K,B=100,u=matrix(runif(B*length(ages)),ncol=B))
{
  # get logLik
  nFossils = length(ages)
  dUN = rep(NA,nFossils )
  for(iObs in 1:nFossils)
  {
    if(sds[iObs]==0)
      e = 0
    else
      e = extraDistr::qtnorm(p = u[iObs,], mean = 0, sd = sds[iObs], a = -Inf, b = ages[iObs]-theta)
    F.eps.m = pnorm(ages[iObs] - theta, mean = 0, sd = sds[iObs])
    F.eps.K = pnorm(K - theta, mean = 0, sd = sds[iObs])
    dUN[iObs] = mean( F.eps.m / F.eps.K / (K-e-theta) )
  }
  return(ll=sum(log(dUN)))
}

# find MLE
getTheta = function(ages, theta=min(ages), eps.sigma, K, B=100, u=matrix(runif(B*length(ages)),ncol=B) )
{
  if(all(eps.sigma==0))
    thetaMLE = list( par=min(ages), value=length(ages)*log(1/(K-min(ages))), hessian=-Inf )
  else
    thetaMLE = optim(theta,UNloglik,ages=ages,sds=eps.sigma,K=K,u=u,method="Brent",lower=min(theta*c(0,2)),upper=max(theta*c(0,2)),control=list(trace=TRUE,fnscale=-1),hessian=TRUE)
  return(thetaMLE)
}

getLRT = function(theta0,thetaMLE,ages, sds, K, alpha=0.05, B=100, u=matrix(runif(B*length(ages)),ncol=B) )
{
  ll0=UNloglik(theta0,ages,sds,K,u=u)
  return(-2*(ll0-thetaMLE$value)-qchisq(1-alpha,1))
}
  
getUNci = function(theta, ages, sds, K, alpha=0.05, B=100, u=matrix(runif(B*length(ages)),ncol=B), wald=FALSE )
{
  thetaMLE = getTheta(ages, theta, sds, K, u=u)
  SE=NULL
  if(wald==TRUE)
  {
    SE = 1/sqrt(-thetaMLE$hessian)
    lo = list( root=thetaMLE$par - qnorm(1-alpha/2) * SE )
    hi = list( root=thetaMLE$par + qnorm(1-alpha/2) * SE )
  }
  else
  {
    lo = try( uniroot(getLRT,thetaMLE$par*c(0.25,1),thetaMLE,alpha=alpha, ages=ages,sd=sds,K=K,u=u,extendInt="downX") )
    if(inherits(lo,"try-error")) lo=list(root=thetaMLE$par)
    hi = try( uniroot(getLRT,thetaMLE$par*c(1,1.25),thetaMLE,alpha=alpha, ages=ages,sd=sds,K=K,u=u,extendInt="upX") )
    if(inherits(hi,"try-error")) hi=list(root=thetaMLE$par)
  }
  return( list( theta=c(lower=lo$root,point=thetaMLE$par,upper=hi$root), B=c(lower=B,point=B,upper=B), se=SE) )
}

do.test = FALSE
if (do.test)
{
K=25000
theta=10000
ages = runif(20, theta, K)
sds = runif(20, 50, 100)
#sd=rep(0,20)
UNloglik(theta,ages,sds,K)

# get MLE
B=100  
u=matrix(runif(B*length(ages)),ncol=B)
thetaMLE=getTheta(ages,theta,sds,K,u=u)

print(getUNci(theta,ages,sds,K,u=u))

# get LRTs
nTheta = 100
thetas=seq(min(ages)-10*sds[1],min(ages)+0.1*sds[1],length=nTheta)
LRs = rep(NA,nTheta)
for (iTheta in 1:nTheta)
  LRs[iTheta] = getLRT(thetas[iTheta],thetaMLE,ages=ages,sds=sds,K=K,u=u)
plot(LRs~thetas,type="l")
print(min(LRs))

}
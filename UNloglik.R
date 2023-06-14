# UNloglik finds the log-likelihood of our uniform-normal model
# given a set of ages and their sds

UNloglik = function(theta,ages,sd,K,B=100,u=matrix(runif(B*length(ages)),ncol=B))
{
  # get logLik
  nFossils = length(ages)
  dUN = rep(NA,nFossils )
  for(iObs in 1:nFossils)
  {
    if(sd[iObs]==0)
      e = 0
    else
      e = extraDistr::qtnorm(p = u[iObs,], mean = 0, sd = sd[iObs], a = -Inf, b = ages[iObs]-theta)
    F.eps.m = pnorm(ages[iObs] - theta, mean = 0, sd = sd[iObs])
    F.eps.K = pnorm(K - theta, mean = 0, sd = sd[iObs])
    dUN[iObs] = mean( F.eps.m / F.eps.K / (K-e-theta) )
  }
  return(ll=sum(log(dUN)))
}

# find MLE
getTheta = function(theta, ages, sd, K, B=100, u=matrix(runif(B*length(ages)),ncol=B) )
{
  if(all(sd==0))
    thetaMLE = list( par=min(ages), value=sum(log(1/(K-min(ages)))) )
  else
    thetaMLE = optim(theta,UNloglik,ages=ages,sd=sd,K=K,u=u,method="Brent",lower=0*theta,upper=2*theta,control=list(trace=TRUE,fnscale=-1))
#   thetaMLE = optim(theta,UNloglik,ages=ages,sd=sd,K=K,u=u,method="BFGS",control=list(trace=TRUE,fnscale=-1))
}

getLRT = function(theta0,thetaMLE,ages, sd, K, alpha=0.05, B=100, u=matrix(runif(B*length(ages)),ncol=B) )
{
  ll0=UNloglik(theta0,ages,sd,K,u=u)
  return(-2*(ll0-thetaMLE$value)-qchisq(1-alpha,1))
}
  
getUNci = function(theta, ages, sd, K, alpha=0.05, B=100, u=matrix(runif(B*length(ages)),ncol=B) )
{
  thetaMLE = getTheta(theta, ages, sd, K, u=u)
  lo = try( uniroot(getLRT,thetaMLE$par*c(0.5,1),thetaMLE,alpha=alpha, ages=ages,sd=sd,K=K,u=u,extendInt="downX") )
  if(inherits(lo,"try-error")) lo=list(root=thetaMLE$par)
  hi = try( uniroot(getLRT,thetaMLE$par*c(1,2),thetaMLE,alpha=alpha, ages=ages,sd=sd,K=K,u=u,extendInt="upX") )
  if(inherits(hi,"try-error")) hi=list(root=thetaMLE$par)
  return( list(theta=c(lower=lo$root,point=thetaMLE$par,upper=hi$root),B=c(lower=B,point=B,upper=B)) )
}

do.test = TRUE
if (do.test)
{
K=25000
theta=10000
ages = runif(20, theta, K)
sd = runif(20, 50, 100)
#sd=rep(0,20)
UNloglik(theta,ages,sd,K)

minmi(
  ages,
  sd = sd,
  alpha = alpha,
  K = K
)

# get MLE
B=100  
u=matrix(runif(B*length(ages)),ncol=B)
thetaMLE=getTheta(theta,ages,sd,K,u=u)

print(getUNci(theta,ages,sd,K,u=u))

# get LRTs
nTheta = 100
thetas=seq(6000,14000,length=nTheta)
LRs = rep(NA,nTheta)
for (iTheta in 1:nTheta)
  LRs[iTheta] = getLRT(thetas[iTheta],thetaMLE,ages=ages,sd=sd,K=K,u=u)
plot(LRs~thetas,type="l")
print(min(LRs))

}
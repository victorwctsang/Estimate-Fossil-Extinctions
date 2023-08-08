# UNloglik finds the log-likelihood of our uniform-normal model
# given a set of ages and their sds

UNaltloglik = function(theta,ages,sds,K)
{
  # get F(K-theta), f(K-theta), f(w-theta)
  lF.eps.w = pnorm(ages - theta, mean = 0, sd = sds, log.p=TRUE)
  F.eps.K = pnorm(K - theta, mean = 0, sd = sds)
  f.eps.K = dnorm(K - theta, mean = 0, sd = sds)
  C = (K-theta) * F.eps.K + sds^2 * f.eps.K
  dlUN = lF.eps.w - log(C)
  return(ll=sum(dlUN))
}

UNaltscore = function(theta,ages,sds,K,sumScores=TRUE)
{
  # get logLik
  F.eps.w = pnorm(ages - theta, mean = 0, sd = sds)
  F.eps.K = pnorm(K - theta, mean = 0, sd = sds)
  f.eps.w = dnorm(ages - theta, mean = 0, sd = sds)
  f.eps.K = dnorm(K - theta, mean = 0, sd = sds)

  C = (K-theta) * F.eps.K + sds^2 * f.eps.K
  score = -f.eps.w / F.eps.w + F.eps.K / C
  if(sumScores==TRUE) score=sum(score)
  return(score)
}

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

UNscore = function(theta,ages,sds,K,B=100,u=matrix(runif(B*length(ages)),ncol=B),sumScores=TRUE)
{
  # get logLik
  nFossils = length(ages)
  dlUN = rep(NA,nFossils )
  for(iObs in 1:nFossils)
  {
    if(sds[iObs]==0)
      e = 0
    else
      e = extraDistr::qtnorm(p = u[iObs,], mean = 0, sd = sds[iObs], a = -Inf, b = ages[iObs]-theta)
    F.eps.m = pnorm(ages[iObs] - theta, mean = 0, sd = sds[iObs])
    F.eps.K = pnorm(K - theta, mean = 0, sd = sds[iObs])
    f.eps.m = dnorm(ages[iObs] - theta, mean = 0, sd = sds[iObs])
    f.eps.K = dnorm(K - theta, mean = 0, sd = sds[iObs])
    dlUNworking = mean( F.eps.m / F.eps.K / (K-e-theta) * ( f.eps.K/F.eps.K + 1/(K-e-theta) ) ) - f.eps.m/F.eps.K/(K-ages[iObs])
    dlUN[iObs] = dlUNworking / UNloglik(theta,ages[iObs],sds[iObs],K,u=as.matrix(u[iObs,]))
  }
  if(sumScores==TRUE) dlUN = sum(dlUN)
  return(dlUN)
}

pUNeps = function(eps,u,theta,K,eps.sigma,n=length(eps.sigma),tol=sqrt(.Machine$double.eps))
# function to compute marginal cdf of epsilon, minus u, to solve for eps
{
  # get CDF denominator
  F.K    = pnorm(K - theta, mean = 0, sd = eps.sigma)
  f.K    = dnorm(K - theta, mean = 0, sd = eps.sigma)
  C      = (K-theta) * F.K + eps.sigma^2 * f.K

  # get CDF-u
  F.eps  = pnorm(eps, mean = 0, sd = eps.sigma)
  f.eps  = dnorm(eps, mean = 0, sd = eps.sigma)
  cdfEps = ( (K-theta)*F.eps + eps.sigma^2*f.eps ) / C - u
  return(cdfEps)
}
  
rUN = function(theta, K, eps.sigma, n=length(eps.sigma), tol=sqrt(.Machine$double.eps),nIter=50)
{
  #ensure sds is the right length
  nSD = length(eps.sigma)
  if(nSD!=n)
  {
    eps.sigma = rep(eps.sigma,length=n)
    warning("length of 'eps.sigma' is not equal to 'n' so will extend 'eps.sigma' as needed.")
  }
  if(nSD==1)  eps.sigma = rep(eps.sigma,length=n)
  
  # compute C and u
  F.K = pnorm(K - theta, mean = 0, sd = eps.sigma)
  f.K = dnorm(K - theta, mean = 0, sd = eps.sigma)
  C = (K-theta) * F.K + eps.sigma^2 * f.K
  u = runif(n)
  uOnC = u * C
  
  # now get cracking finding eps
  epsOld = qnorm(u,mean=0,sd=eps.sigma) # starting estimate
  eps = quantU = rep(0,n) 
  quantUMax = pnorm(K-theta,sd=eps.sigma) # quantU can't get any larger than this or eps will be larger than K-theta
  iter=0
  isDiff = eps.sigma!=0 #only do the below when eps.sigma is non-zero
  while(any(isDiff) & iter<nIter)
  {
    quantU[isDiff] = ( uOnC[isDiff] - eps.sigma[isDiff]^2*dnorm(epsOld[isDiff], mean = 0, sd = eps.sigma[isDiff]) ) / (K-theta)
    quantU[isDiff] = pmax(quantU[isDiff], sqrt(tol)) # correction for wild estimates sending quantU negative
    quantU[isDiff] = pmin(quantU[isDiff], quantUMax[isDiff]-sqrt(tol)) # correction for wild estimates sending eps over K-theta
    eps[isDiff]    = qnorm(quantU[isDiff], mean=0, sd=eps.sigma[isDiff])
    epsDiff        = abs(eps-epsOld)
    isDiff         = epsDiff>tol
    epsOld[isDiff] = eps[isDiff]
    iter           = iter+1
  }

  if(iter==nIter)
  {
    for(iObs in which(isDiff))
    {
      epsTry = try( uniroot( pUNeps, interval=c(qnorm(sqrt(tol),mean=0,sd=eps.sigma[iObs]),K-theta),u=u[iObs],theta=theta,eps.sigma=eps.sigma[iObs],K=K,extendInt="upX") )
      if(inherits(epsTry,"try-error")==FALSE)
      {
        eps[iObs] = epsTry$root
        isDiff[iObs] = FALSE
      }
    }
    if(any(isDiff)) warning(paste0("non-convergence for ",sum(isDiff)," observations"))
  }
  X = runif(n,min=theta,max=K-eps)
  W = X + eps
  return(W)
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

getThetaAlt = function(ages, theta=min(ages), eps.sigma, K )
{
  if(all(eps.sigma==0))
    thetaMLE = list( par=min(ages), value=length(ages)*log(1/(K-min(ages))), hessian=-Inf )
  else
    thetaMLE = optim(theta,UNaltloglik,ages=ages,sds=eps.sigma,K=K,method="Brent",lower=min(theta*c(0,2)),upper=max(theta*c(0,2)),control=list(trace=TRUE,fnscale=-1),hessian=TRUE)
  return(thetaMLE)
}


getLRTalt = function(theta0,thetaMLE,ages, sds, K, alpha=0.05)
{
  ll0=UNaltloglik(theta0,ages,sds,K)
  return(-2*(ll0-thetaMLE$value)-qchisq(1-alpha,1))
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

getUNciAlt = function(theta, ages, sds, K, alpha=0.05, wald=FALSE )
{
  thetaMLE = getThetaAlt(ages, theta, sds, K)
  SE=NULL
  if(wald==TRUE)
  {
    SE = 1/sqrt(-thetaMLE$hessian)
    lo = list( root=thetaMLE$par - qnorm(1-alpha/2) * SE )
    hi = list( root=thetaMLE$par + qnorm(1-alpha/2) * SE )
  }
  else
  {
    lo = try( uniroot(getLRTalt,thetaMLE$par*c(0.25,1),thetaMLE,alpha=alpha, ages=ages,sd=sds,K=K,extendInt="downX") )
    if(inherits(lo,"try-error")) lo=list(root=thetaMLE$par)
    hi = try( uniroot(getLRTalt,thetaMLE$par*c(1,1.25),thetaMLE,alpha=alpha, ages=ages,sd=sds,K=K,extendInt="upX") )
    if(inherits(hi,"try-error")) hi=list(root=thetaMLE$par)
  }
  return( list( theta=c(lower=lo$root,point=thetaMLE$par,upper=hi$root), B=c(lower=NA,point=NA,upper=NA), se=SE) )
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

UNaltloglik(theta,ages,sds,K)


# get MLE
B=100  
u=matrix(runif(B*length(ages)),ncol=B)
thetaMLE=getTheta(ages,theta,sds,K,u=u)
thetaaltMLE=getThetaAlt(ages,theta,sds,K)

print(getUNci(theta,ages,sds,K,u=u))

print( sum( UNscore(thetaMLE$par,ages,sds,K) ) )
print( sum( UNaltscore(thetaaltMLE$par,ages,sds,K) ) )


# get LRTs
nTheta = 100
thetas=seq(min(ages)-10*sds[1],min(ages)+0.1*sds[1],length=nTheta)
LRs = scores = LRalts = scoreAlts = rep(NA,nTheta)
for (iTheta in 1:nTheta)
{
  LRs[iTheta] = getLRT(thetas[iTheta],thetaMLE,ages=ages,sds=sds,K=K,u=u)
  scores[iTheta] = sum(UNscore(thetas[iTheta],ages=ages,sds=sds,K=K,u=u))
  LRalts[iTheta] = getLRTalt(thetas[iTheta],thetaMLE,ages=ages,sds=sds,K=K)
  scoreAlts[iTheta] = sum(UNaltscore(thetas[iTheta],ages=ages,sds=sds,K=K))
}
par(mfrow=c(2,2),mar=c(3,3,1,1),mgp=c(1.75,0.75,0))
plot(LRs~thetas,type="l")
plot(scores~thetas,type="l")
abline(h=0,col="red")
abline(v=thetaMLE$par,col="red")

plot(LRalts~thetas,type="l")
plot(scoreAlts~thetas,type="l")
abline(h=0,col="red")
abline(v=thetaaltMLE$par,col="red")
print(min(LRs))

rUN(10000,25000,sds,1000)
}
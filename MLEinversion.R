## Testing regression-based inversion code (regInversion.R)

source("regInversion.R")
source("UNloglik.R")

simFn = function (theta, K, eps.mean = 0, eps.sigma, n=length(eps.sigma),B=NULL,u=NULL,trans=FALSE)
{
  # Simulate fossils assuming:
  # - Gaussian measurement error (truncated at K-theta)
  # - Uniform deposition from theta to K-eps
  if(trans) theta=K-exp(-theta)
  if(theta>K) #trying to game it to push estimates away from K 
    W=rep(theta,n)
  else
  {
    eps = extraDistr::rtnorm(n, mean = 0, sd = eps.sigma, a = -Inf, b = K-theta)
    eps[eps.sigma==0] = 0 #dealing with NAs from rtnorm
    X = runif(n, min = theta, max = K-eps)
    W = X + eps
  }
  return(W)
}

getTh = function (ages, theta=NULL, eps.sigma, K, B=100,u=matrix(runif(B*length(ages)),ncol=B), trans=FALSE )
{
  if(is.null(theta))  theta=min(ages)
  else
  {
    if(trans) theta=K-exp(-theta)
  }
  if(all(ages>K))
    theta = min(ages)
  else
    theta = getTheta(ages, theta=theta, eps.sigma, K, B=100, u=u )$par
  return(theta)
}

do.test = FALSE
if (do.test)
{
  K=25000
  theta=10000
  trans=FALSE
  ages = runif(25, theta, K)
  sds = runif(25, 50, 100)
  B=100
  u=matrix(runif(B*length(ages)),ncol=B)
  #sd=rep(0,20)
#  thetaInits = min(ages)+max(sds)*seq(-20,10,length=10)
  thetaInits = min(ages) + mean(sds)*seq(-5,5,length=10)
  if(trans) thetaInits = -log(K-thetaInits)
  time.start = Sys.time()
  ft=regInversion(ages,getT=getTh,simulateData=simFn,thetaInits=thetaInits,
                 q=0.025,iterMax=500,K=K,eps.sigma=sds, u=u, trans=trans,method="rq")
  time.stop = Sys.time()
  if(trans) ft$theta=K-exp(-ft$theta)
  print(list(elapsed=time.stop-time.start,thetaEst=ft$theta,err=ft$err,iterations=ft$iter,converged=ft$converged))
}

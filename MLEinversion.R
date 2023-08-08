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

simFnAlt = function (theta, K, eps.sigma, n=length(eps.sigma))
{
  # Simulate fossils assuming:
  # - Gaussian measurement error (truncated at K-theta)
  # - Uniform deposition from theta to K-eps
  if(theta>K) #trying to game it to push estimates away from K 
    W=rep(theta,n)
  else
    W = rUN(theta, K, eps.sigma, n=length(eps.sigma))
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

getThAlt = function (ages, theta=NULL, K, eps.sigma )
{
  if(is.null(theta))  theta=min(ages)
  if(all(ages>K))
    theta = min(ages)
  else
    theta = getThetaAlt(ages, theta=theta, eps.sigma, K )$par
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
  thetaInits = min(ages) + mean(sds)*seq(-5,5,length=20)
  if(trans) thetaInits = -log(K-thetaInits)
  time.start = Sys.time()
  q=0.025
  ft=regInversion(ages,getT=getTh,simulateData=simFn,thetaInits=thetaInits,
                 q=q,iterMax=500,K=K,eps.sigma=sds, u=u, trans=trans,method="rq")
  fta=regInversion(ages,getT=getThAlt,simulateData=rUN,thetaInits=thetaInits,
                  q=q,iterMax=500,K=K,eps.sigma=sds, method="rq")
  ft2=regInversion(ages,getT=getTh,simulateData=simFn,thetaInits=thetaInits,
                  q=q,iterMax=500,K=K,eps.sigma=sds, u=u, trans=trans,method="wrq")
  fta2=regInversion(ages,getT=getThAlt,simulateData=rUN,thetaInits=thetaInits,
                   q=q,iterMax=500,K=K,eps.sigma=sds, method="wrq")
  time.stop = Sys.time()
  if(trans) ft$theta=K-exp(-ft$theta)
  print(list(elapsed=time.stop-time.start,thetaEst=ft$theta,err=ft$err,iterations=ft$iter,converged=ft$converged))

  print(c(ft$theta,fta$theta,ft2$theta,fta2$theta))
  par(mfrow=c(2,2))
  plot(T~theta,data=ft$stats)
  plot(T~theta,data=fta$stats)
  plot(T~theta,data=ft2$stats)
  plot(T~theta,data=fta2$stats)
}

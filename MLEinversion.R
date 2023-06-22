## Testing regression-based inversion code (regInversion.R)

source("regInversion.R")
source("UNloglik.R")

simFn = function (eta, K, eps.mean = 0, eps.sigma = 0, n=length(eps.sigma))
{
  # Simulate fossils assuming:
  # - Gaussian measurement error (truncated at K-theta)
  # - Uniform deposition from theta to K-eps
  theta=K-exp(-eta)
  eps = extraDistr::rtnorm(n, mean = 0, sd = eps.sigma, a = -Inf, b = K-theta)
    X = runif(n, min = theta, max = K-eps)
    W = X + eps
  return(W)
}

getTh = function (ages, eta=NULL, eps.sigma = 0, K, u=matrix(runif(B*length(ages)),ncol=B) )
{
  if(is.null(eta))
    theta=min(ages)
  else
    theta=K-exp(-eta)
  return( getTheta(ages, theta=theta, eps.sigma, K, B=100, u=matrix(runif(B*length(ages)),ncol=B) )$par )
}
  
do.test = TRUE
if (do.test)
{
  K=25000
  theta=10000
  ages = runif(25, theta, K)
  sds = runif(25, 50, 100)
  B=100
  u=matrix(runif(B*length(ages)),ncol=B)
  #sd=rep(0,20)
#  thetaInits = min(ages)+max(sds)*seq(-20,10,length=10)
  thetaInits = seq(min(ages),max(ages),length=10)
  etaInits = -log(K-thetaInits)
  time.start = Sys.time()
  ft=regInversion(ages,getT=getTh,simulateData=simFn,thetaInits=etaInits,
                 q=0.025,iterMax=500,K=K,eps.sigma=sds)
  time.stop = Sys.time()
  print(list(elapsed=time.stop-time.start,theta=K-exp(-ft$theta),iterations=ft$iter,converged=ft$converged))
}

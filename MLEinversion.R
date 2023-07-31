## Testing regression-based inversion code (regInversion.R)

source("regInversion.R")
source("UNloglik.R")

simFn = function (theta, K, eps.mean = 0, eps.sigma, B=NULL,u=NULL,trans=FALSE, n=NULL)
{
  # Simulate fossils assuming:
  # - Gaussian measurement error (truncated at K-theta)
  # - Uniform deposition from theta to K-eps
  if(trans) theta=K-exp(-theta)
  if(theta>K) #trying to game it to push estimates to more sensible places 
    dat=list(W=theta,eps=0)
  else
  {
    dat = simulate_dataset(theta, K, eps.mean=0, eps.sigma=eps.sigma, n=n)
  }
  return(dat)
}

getTh = function (dat, theta=NULL, eps.sigma, K, B=100,u=matrix(runif(B*length(dat$W)),ncol=B), trans=FALSE, n=length(dat$W) )
{
  if(is.null(theta))  theta=min(dat$W)
  if(trans) theta=K-exp(-theta)
  if(all(dat$W>K))
    theta = min(dat$W)
  else
    theta = getTheta(dat$W, theta=theta, dat$eps, K, B=100, u=u )$par
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
  ftSD = glm(sds~ages,family=Gamma("log"))
  B=100
  u=matrix(runif(B*length(ages)),ncol=B)
  #sd=rep(0,20)
#  thetaInits = min(ages)+max(sds)*seq(-20,10,length=10)
  thetaInits = min(ages) + mean(sds)*seq(-5,5,length=20)
  if(trans) thetaInits = -log(K-thetaInits)
  time.start = Sys.time()
  q=0.025
  dat = list(W=ages,eps=sds)
  ftS=regInversion(dat,getT=getTh,simulateData=simFn,thetaInits=thetaInits,
                  q=q,iterMax=500,K=K,eps.sigma=ftSD, u=u, trans=trans,method="rq",aMean=0, n=length(dat$W))
  ft=regInversion(dat,getT=getTh,simulateData=simFn,thetaInits=thetaInits,
                  q=q,iterMax=500,K=K,eps.sigma=sds, u=u, trans=trans,method="rq",aMean=0)
  ftwS=regInversion(dat,getT=getTh,simulateData=simFn,thetaInits=thetaInits,
                  q=q,iterMax=500,K=K,eps.sigma=ftSD, u=u, trans=trans,method="wrq",aMean=0, n=length(dat$W))
#  ftwS=regInversion(dat,getT=getTh,simulateData=simFn,thetaInits=thetaInits,
#                   q=q,iterMax=500,K=K,eps.sigma=sds, u=u, trans=trans,method="wrq",aMean=0)
  time.stop = Sys.time()
  if(trans) ft$theta=K-exp(-ft$theta)
  print(list(elapsed=time.stop-time.start,thetaEst=ft$theta,err=ft$err,iterations=ft$iter,converged=ft$converged))

  print(c(ft$theta,ftS$theta,ftwS$theta))
  par(mfrow=c(2,2))
  plot(T~theta,data=ft$stats,main="rq")
  plot(T~theta,data=ftS$stats,main="rq,ftSD")
  plot(T~theta,data=ftwS$stats,main="wrq,ftSD")
#  plot(T~theta,data=ftp$stats,main="wrq,a=1")
}

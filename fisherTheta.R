fisherTheta = function(data, getT, simulateData, thetaInits, q=0.5, iterMax=1000, eps=1.e-3, method="rq", ...)
{
  t_obs = getT(data, ...)

  # get initial stats
  stats=data.frame(theta=NULL,T=NULL,thetaEst=NULL)
  for(i in 1:length(thetaInits) )
  {
    newDat = simulateData(thetaInits[i], ...)
    stats = rbind(stats, c(theta=thetaInits[i], T = getT(newDat, ...), thetaEst=thetaInits[i]) )
  }      
  names(stats)=c("theta","T","thetaEst")
  
  # define getPred and getErr functions (to avoid repeated if statements in estimation)
  if(method=="rq")
  {
    q = 1-q # flipping quantile is needed when using this approach
    getPred = function(newtheta,qfit,t_obs,q)
    {
      return(predict(qfit,newdata=list(theta=newtheta))-t_obs)
    }
    getErr = function(newtheta,qfit,t_obs,q) # getting the SE of predictions 
    {
      prEst = try(predict(qfit,newdata=list(theta=newtheta),interval="confidence"),silent=TRUE)
      if(inherits(prEst,"try-error"))
        err = Inf
      else
        err = (prEst[3]-prEst[2])/qnorm(0.975)/2/coef(qfit)[2]
      return(err)
    }
  }
  else
  {
    getPred = function(newtheta,qfit,t_obs,q)
    {
      return(predict(qfit,newdata=list(theta=newtheta),type="response")-q)
    }
    getErr = function(newtheta,qfit,t_obs,q)
    {
      prEst = predict(qfit,newdata=list(theta=newtheta),se.fit=TRUE)
      err = prEst$se.fit/coef(qfit)[2]
      return(err)
    }
  }
  
  iter = length(thetaInits)
  res = updateTheta(stats, t_obs, q, method=method, getPred=getPred)
  isConverged = FALSE
  while(isConverged == FALSE & iter<iterMax)
  {
    iter=iter+1
    thetaNew = iter*res$theta - sum(stats$theta)
    newDat = simulateData(thetaNew, ...)
    stats = rbind(stats, c(theta=thetaNew, T = getT(newDat, ...), thetaEst=res$theta) )
    res = updateTheta(stats, t_obs, q, qfit=res$qfit, method=method, getPred=getPred)
    err = getErr(res$theta,res$qfit,t_obs,q)
    isConverged = err<eps
  }
  return(list(theta=res$theta,error=err,iter=iter,converged=isConverged,stats=stats))
}

updateTheta = function(stats,t_obs,q,qfit=NULL,method="rq",getPred=getPred)
# dat is a dataframe containing theta and T (parameter and statistic)
{
  if(is.null(qfit))
  {
    if(method=="rq")
      ft = quantreg::rq(T~theta,tau=q,data=stats)
    else
      ft = glm(I(T>t_obs)~theta,family=binomial("probit"),data=stats)
  }
  else
     ft = update(qfit,data=stats)
  theta = uniroot(getPred,lower=min(stats$theta), upper=max(stats$theta), extendInt="upX",qfit=ft,t_obs=t_obs,q=q) #find the value matching to observed t
  return(list(theta=theta$root,qfit=ft))
}

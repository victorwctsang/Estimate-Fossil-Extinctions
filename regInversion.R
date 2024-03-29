regInversion = function(data, getT, simulateData, thetaInits, q=0.5, iterMax=1000, eps=1.e-5, method="rq", aMean=0, stats = NULL, ...)
{
  t_obs = getT(data, ...)

  # get initial stats
  if(is.null(stats))
  {
    stats=data.frame(theta=NULL,T=NULL,thetaEst=NULL)
    for(i in 1:length(thetaInits) )
    {
      newDat = simulateData(thetaInits[i], ...)
      stats = rbind(stats, c(theta=thetaInits[i], T = getT(newDat, ...), thetaEst=thetaInits[i], wt=NULL) )
    }      
    names(stats)=c("theta","T","thetaEst")
  }
  
  # define getPred and getErr functions (to avoid repeated if statements in estimation)
  if(method=="rq"||method=="wrq"||method=="rq2")
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
      {
        if(method=="rq2")
          grad = 2*coef(qfit)[2]*newtheta + coef(qfit)[3]
        else
          grad = coef(qfit)[2]
        if(grad>0) # find the value matching to observed t if decent fit
          err = (prEst[3]-prEst[2])/qnorm(0.975)/2/grad
        else
          err = Inf
      }
      if(is.na(err)) err=Inf
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
      if(coef(qfit)[2]>0 & qfit$converged) # find the value matching to observed t if decent fit
      {
        prEst = try(predict(qfit,newdata=list(theta=newtheta),se.fit=TRUE),silent=TRUE)
        if(inherits(prEst,"try-error"))
          err = Inf
        else
          err = prEst$se.fit/coef(qfit)[2]
      }
      else 
        err = Inf
      return(err)
    }
  }
  
  iter = dim(stats)[1]
  res = updateTheta(stats, t_obs, q, method=method, getPred=getPred)
  isConverged = FALSE
  while(isConverged == FALSE & iter<iterMax)
  {
    iter=iter+1
    thetaNew = getThetaSim(iter,thetaEst=res$theta,thetaSims=stats$theta,aMean=aMean)
    newDat = simulateData(thetaNew, ...)
    Titer = getT(newDat, ...)
    if(is.na(Titer)|Titer==Inf) 
    {
      warning("getT not a number at current value of theta, resetting to median")
      thetaNew = median(stats$theta)
      newDat = simulateData(thetaNew, ...)
      Titer = getT(newDat, ...)
    }
    stats = rbind(stats, c(theta=thetaNew, T=Titer, thetaEst=res$theta) )
    res = updateTheta(stats, t_obs, q, qfit=res$qfit, method=method, getPred=getPred)
    err = getErr(res$theta,res$qfit,t_obs,q) / abs(res$theta)
    isConverged = err < eps
  }
  return(list(theta=res$theta,error=err,iter=iter,converged=isConverged,stats=stats,fit=res$qfit))
}

updateTheta = function(stats,t_obs,q,qfit=NULL,method="rq",getPred=getPred,screenStats=NULL)
# dat is a dataframe containing theta and T (parameter and statistic)
{
  if(method=="wrq")
  {
    lm_ft = lm(T~theta,data=stats)
    infl = influence(lm_ft)$hat
    stats$wt[is.na(stats$T)==FALSE] = min(infl,2/length(stats$theta)) / infl
  }

  if(is.null(qfit))
  {
    ft = switch(method,
           "rq" = quantreg::rq(T~theta,tau=q,data=stats), #consider varying rq method, "fn" or "pfn"
           "rq2" = quantreg::rq(T~poly(theta,2,raw=TRUE),tau=q,data=stats), #consider varying rq method, "fn" or "pfn"
#           "qgam" = qgam::qgam(T~s(theta), qu=q, data=stats),
           "wrq" = quantreg::rq(T~theta,tau=q,weights=wt,data=stats),
               glm(I(T>t_obs)~theta,family=binomial("probit"),data=stats)
    )
  }
  else
  {
    #update fit only using data points that are "local" to observed stat (predicted T within factor of screenStats of observed T)
    if(is.null(screenStats)==FALSE) 
    {
      prRatio = abs( log(predict(qfit,newdata=stats)/t_obs) )
      prRatio[is.na(prRatio)] = Inf
      OKstats = prRatio<log(screenStats)
      if(sum(OKstats)<10) #if this does not keep enough, take 5 closest
        OKstats=sort(prRatio,index.return=TRUE)$ix[1:10]
      statsWorking = stats[OKstats,]
    }
    else
      statsWorking=stats
    #    if(method=="qgam")   
#      ft = update(qfit,data=stats,family=qfit$family)
#    else
      ft = update(qfit,data=statsWorking)
  }
  if(method=="rq2")
  {
    # solve quadratic equation for larger solution if it exists and if stationary point is below first quartile
    Delta = coef(ft)[2]^2 - 4 * coef(ft)[3] * ( coef(ft)[1]-t_obs )
#    statPt = -coef(ft)[2]/2/coef(ft)[3]
#    cond = Delta>0 & statPt<quantile(stats$theta,0.25)
    cond = Delta>0
    if(cond)
      theta = ( -coef(ft)[2] + sqrt(Delta) ) / 2 / coef(ft)[3]
  }
  else
  {
    # solve linear equation if increasing slope
    if(method=="prob") t_obs=qnorm(q) # if using probit regression
    cond = coef(ft)[2]>0
    if(cond) # find the value matching to observed t if increasing fit
      theta = (t_obs-coef(ft)[1])/coef(ft)[2] 
  }
  if(cond == FALSE) # if fit is screwy use current best estimate + jitter until something changes!
     theta = stats$thetaEst[nrow(stats)]*(1+runif(1,-0.01,0.01))
  # don't accept any crazy shit... move them back towards the data
  # theta = noOutliers(theta,stats$theta)
  return(list(theta=c(theta,use.names=FALSE),qfit=ft))
}

getThetaSim = function(iter,thetaEst,thetaSims,aMean)
{
  # aMean tells us how much to weight mean of all obs compared to current obs
  thetaNew = (aMean*iter+1-aMean)*thetaEst - aMean*sum(thetaSims)
  # put bounds on thetaNew so no crazy shit
  if(aMean>0) thetaNext = noOutliers(thetaNew,thetaSims)
  thetaNext=thetaNew
  return(thetaNext)
}
noOutliers = function(thetaNew,thetas)
# don't accept any crazy shit from theta - truncate using an asymmetric 3IQR rule
{
  quants = quantile(thetas,c(0.25,0.5,0.75))
  outRange = diff(quants)*3
  if ( thetaNew>quants[3]+outRange[2] )
    thetaNew = quants[3]+outRange[2]
  if ( thetaNew<quants[1]-outRange[1] )
    thetaNew = quants[1]-outRange[1]
  return(thetaNew)
}
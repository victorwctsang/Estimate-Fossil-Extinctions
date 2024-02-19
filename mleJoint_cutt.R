mleJoint_cutt = function(ages, sd, K, dfInit=4, alpha=0.05, q=c(alpha/2,1-alpha/2), wald=FALSE, coord=TRUE, ...)
{  
  nSD = length(sd)
  n = length(ages)
  if(nSD==1) sd = rep(sd,n)
  if(nSD!=n & nSD!=1)
  {
    sd = rep(sd,length=n)
    warning("lengths of 'ages' and sd' are not equal so will extend 'sd' as needed.")
  }
  
  if(all(ages>K)) stop("'ages' need to be no larger than 'K'")
  if(any(ages>K))
  {
    warning("Some ages exceed 'K', these will be ignored")
    sd   = sd[ages<=K]
    ages = ages[ages<=K]
  }
  mles = getMLE(ages=ages, theta=min(ages), sd=sd, K=K, df=dfInit, coord=coord)
  vr = solve(-mles$hessian)
  SE = if(is.nan(vr[1,1])) 0 else sqrt(vr[1,1])

  if(is.null(alpha))
  {
    result = list(theta=mles$par[1], df=1/mles$par[2], se=SE, call=match.call())
  }
  else
  {
    # set up result list.
    ci=rep(NA,length(q))
    if(is.null(names(q)))
      names(ci)=paste0("q=",q)
    else
      names(ci)=names(q)
    
    if(wald==TRUE)
    {
      ci = mles$par[1] + qnorm(q) * SE
    }
    else
    {
      nQ = length(q)
      q2Tail      = 2*pmin(q,1-q)
      # set search limits so that we look above MLE if q>0.5 and below otherwise 
      is_SE_bad   = is.nan(SE) | is.infinite(SE) | SE==0
      searchLim   = ifelse( is_SE_bad, IQR(ages)*0.5, SE*5 )
      qLo = qHi   = rep(mles$par[1],nQ)
      qLo[q<=0.5] = mles$par[1]-searchLim
      qHi[q>=0.5] = min(mles$par[1]+searchLim,K)
      # note LRT function is increasing for q>0.5 
      dir         = rep("downX",nQ)
      dir[q>=0.5] ="upX"
      for(iQ in 1:nQ)
      {
        thLim = try( uniroot(cutt_LRTprofile, c(qLo[iQ],qHi[iQ]), mles, alpha=q2Tail[iQ], ages=ages, sd=sd, K=K, extendInt=dir[iQ], ...) )
        if(inherits(thLim,"try-error"))
          ci[iQ] = mles$par[1]
        else
          ci[iQ] = thLim$root
      }
    }
    result = list( theta=mles$par[1], ci=ci, se=SE, df=1/mles$par[2], q=q, call=match.call())
  }
  class(result)="mle_cutt"
  return( result )
}

getMLE = function(ages, theta=min(ages), sd, K, df=4, coord=TRUE, nIter=10 )
{
  if(all(sd==0))
    MLE = list( par=c(min(ages),Inf), value=length(ages)*log(1/(K-min(ages))), hessian=matrix(-Inf,2,2) )
  else
  {
    if(coord==TRUE)
    {
      iIter = 1
      cond = FALSE
      MLE  = optim(theta, reginv:::cutt_LogLik, 
                ages=ages, sd=sd, K=K, df=df, method="Brent",
                lower=-1/sqrt(.Machine$double.eps), upper=K, control=list(trace=TRUE,fnscale=-1))
      while(cond==FALSE)
      {
        pre = MLE
        df  = getDF( ages, theta=MLE$par, sd=sd, K=K, dfInvInit=1/df )$par
        MLE = optim( pre$par, reginv:::cutt_LogLik, 
                 ages=ages, sd=sd, K=K, df=df, method="Brent",
                 lower=-1/sqrt(.Machine$double.eps), upper=K, control=list(trace=TRUE,fnscale=-1) )
        eps = abs(pre$value-MLE$value)
        cond = eps>1.e-5 & iIter<=nIter
        iIter = iIter+1
      }
      if(cond==FALSE) MLE$convergence = 1
#      MLE = optim(theta, cutt_LRTprofile, mles=list(value=0,pars=c(0,1/df)), ages=ages, sd=sd, K=K, alpha=0, method="Brent",
#                   lower=-1/sqrt(.Machine$double.eps), upper=K, control=list(trace=TRUE,fnscale=1))
#      df  = getDF(ages,theta=MLE$par, sd=sd, K=K, dfInvInit = 1/df)$par
      MLE$par = c( MLE$par, 1/df )
      MLE$hessian = optimHess( MLE$par, cutt_LogLikJoint, ages=ages, sd=sd, K=K, control=list(trace=TRUE,fnscale=-1) )
    }
    
    else
    {    
    MLE = optim(c(theta,1/df), cutt_LogLikJoint, ages=ages, sd=sd, K=K, method="L-BFGS-B",
                lower=c(-1/sqrt(.Machine$double.eps),0.005),
                upper=c(K,0.5-sqrt(.Machine$double.eps)),
                control=list(trace=TRUE,fnscale=-1),hessian=TRUE)
    }
  }
  return(MLE)
}

getDF = function( ages, theta, sd, K, dfInvInit=1/4 )
{
  if(all(sd==0))
    res = list( par=Inf, value=length(ages)*log(1/(K-min(ages))) )
  else
  {
    dfInv = optim( dfInvInit, cutt_LogLikT, ages=ages, sd=sd, theta=theta, K=K, method="Brent", 
               lower=0.005, upper=0.5-sqrt(.Machine$double.eps), control=list(trace=TRUE,fnscale=-1) )
    res = list(par=1/dfInv$par,value=dfInv$value)
  }
  return(res)
}

cutt_LogLikJoint = function(params,ages,sd,K) #parameters are (theta, dfInv) 
{
  if(params[2]>=0.5)
    ll=sum(dcutt(x=ages,theta=params[1],K=K,sd=sd,df=2+sqrt(.Machine$double.eps),log=TRUE))*(params[2]+0.5) # game it away from df=2
  else
    ll=sum(dcutt(x=ages,theta=params[1],K=K,sd=sd,df=1/params[2],log=TRUE))
  return(ll)
}

cutt_LogLikT = function(dfInv,ages,sd,theta,K)
{
  if(dfInv>=0.5)
    ll=sum(dcutt(x=ages,theta=theta,K=K,sd=sd,df=2+sqrt(.Machine$double.eps),log=TRUE))*(dfInv+0.5) # game it away from df=2
  else
    ll=sum(dcutt(x=ages,theta=theta,K=K,sd=sd,df=1/dfInv,log=TRUE))
  return(ll)
}

cutt_LRTprofile = function(theta0,mles,ages, sd, K, alpha=0.05)
{
  ll0  = getDF(ages,theta=theta0, sd=sd, K=K, dfInvInit = mles$par[2])$value
  return(-2*(ll0-mles$value)-qchisq(1-alpha,1))
}

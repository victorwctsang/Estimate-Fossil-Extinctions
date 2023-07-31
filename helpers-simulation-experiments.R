source("MLEinversion.R")

library(readxl)

source("GRIWM.R")

getFossilData = function (path, sheet, range, col_names, col_types) {
  read_excel(path, sheet, range, col_names, col_types)
}

boundFossilData = function (fossil.df, K) {
  return(fossil.df[fossil.df$age < K,])
}

getStdDevFromFossilData = function (path, K, n, ...) {
  df = getFossilData(path, ...)
  df = boundFossilData(df, K)
  return(df$sd[1:n])
}

create_result_df = function (n.sims, method, error_factor) {
  n.methods = length(method)
  n.error_factors = length(error_factor)
  df.length = n.sims * n.methods * n.error_factors
  
  data.frame(
    sim_id = rep(1:n.sims, each = (n.methods * n.error_factors)),
    method = rep(method, each = n.error_factors, length.out = df.length),
    error_factor = rep(error_factor, length.out = df.length),
    lower = rep(NA, df.length),
    upper = rep(NA, df.length),
    mean = rep(NA, df.length),
    width = rep(NA, df.length),
    contains_theta = rep(NA, df.length),
    compute_time = rep(NA, df.length)
  )
}

simulate_sd = function(sd_fit,thetas)
{
  newDat = list(thetas)
  names(newDat) = names(model.frame(sd_fit))[2]
  mus    = predict(sd_fit,newdata=newDat,type="response")
# I think we should use predicted sd not observed (treating random noise as random)
#  shape <- MASS::gamma.shape(sd_fit)$alpha
#  rgamma(length(mus), shape = shape, rate = shape/mus)
}
  
simulate_dataset = function (theta.true,
                             K,
                             eps.mean,
                             eps.sigma,
                             n=NULL,
                             error_factor=1)
{
  
if(inherits(eps.sigma,"glm")) #if eps.sigma is a model then use it to simulate data
{
  mf = model.frame(eps.sigma)
  if(is.null(n)) n = length(mf$sd)
  if(error_factor==0)
  {
           W = runif(n, min = theta.true, max = K)
    sigmaNew = rep(0,n)
  }
  else
  {
    if(error_factor!=1) #if there is an error factor, multiply sds and refit model
    {
      mf[,1]    = mf[,1] * error_factor
      eps.sigma = update(eps.sigma,data=mf)
    }
    # compute 3 predicted sds above K (will simulate to here)
    newDat = list(K)
    names(newDat) = names(mf)[2]
    sdK    = predict(eps.sigma,newdata=newDat,type="response")
    KLim   = K + 3*sdK
    if(theta.true>K) # to discourage impossible values of theta
    {
             W = rep(theta.true,n)
      sigmaNew = rep(0,n)
    }   
    else
    {
      sigmaNew = X = W = eps = rep(NA,n)
         isOut = rep(TRUE,n)
          while(any(isOut))
      {
                 nOut = sum(isOut)
             X[isOut] = runif(nOut, min = theta.true, max = KLim)
      sigmaNew[isOut] = simulate_sd(eps.sigma,X[isOut])
           eps[isOut] = rnorm(nOut, sd = sigmaNew[isOut] )
             W[isOut] = X[isOut] + eps[isOut]
                isOut = W>K
      }
    }
  }
}
else
{
  if(is.null(n)) n = length(eps.sigma)
  sigmaNew = eps.sigma * error_factor
  if (any(sigmaNew != 0)) {
    eps = extraDistr::rtnorm(n, mean = eps.mean, sd = sigmaNew, a=-Inf, b=(K-theta.true))
  }
  else
    eps = eps.mean
  X = runif(n, min = theta.true, max = K-eps)
  W = X + eps
}
return(list(W=W,eps=sigmaNew))
}

simulate_datasets = function (config) {
  df.length = config$n.trials * config$n.error_factors
  
  dataset.df = list(
    error_factor = rep(config$error_factor, length.out = df.length),
               W = matrix(NA,config$n.samples,df.length),
             eps = matrix(NA,config$n.samples,df.length)
  )
  for(iTrial in 1:df.length)
  {
    dati = simulate_dataset(theta.true = config$theta.true,
                     K = config$K,
                     eps.mean = config$dating_error.mean,
                     eps.sigma = config$fossil.sd,
                     n = config$n.samples,
                     error_factor=dataset.df$error_factor[iTrial])
    dataset.df$W[,iTrial]   = dati$W
    dataset.df$eps[,iTrial] = dati$eps
  }
  return(dataset.df)
}

estimate_extinction = function (W, sd, method, K, dating_error.mean) {
  estimate = NA
  runtime = NA
  start_time = Sys.time()
  estimate = switch(
    method,
    MLE = mle(W),
    `BA-MLE` = ba_mle(W, K),
    Strauss = strauss(W)
  )
  runtime = calculate_tdiff(start_time, Sys.time())
  return(list(point = estimate, point_runtime = runtime))
}

calculate_tdiff = function(start, end) {
  as.numeric(difftime(end, start), units = "secs")
}

estimate_conf_int = function (W,
                              sd,
                              method,
                              alpha,
                              K,
                              dating_error.mean,
                              sd_model) {
  estimate = switch(
    method,
    GRIWM = griwm(
      alpha = alpha,
      dates = W,
      sd = sd,
      K = K,
      p_t = alpha,
      bias_adjusted = F
    ),
    `GRIWM-corrected` = griwm(
      alpha = alpha,
      dates = W,
      sd = sd,
      K = K,
      p_t = 0.5,
      bias_adjusted = T
    ),
    MINMI = do_minmi(
      ages = W,
      sd = sd,
      alpha = alpha,
      K = K
    ),
    UNci = do_UNci(theta.true,
      ages = W,
      sd = sd,
      alpha = alpha,
      K = K,
      wald=FALSE
    ),
    UNwald = do_UNci(theta.true,
                   ages = W,
                   sd = sd,
                   alpha = alpha,
                   K = K,
                   wald=TRUE
    )
    ,
    mleInv = do_mleInv(ages = W,
                     sd = sd,
                     alpha = alpha,
                     K = K)
    ,
    mleInvS = do_mleInv(ages = W,
                       sd = sd,
                       alpha = alpha,
                       K = K,
                       doMod=TRUE)
    ,
    mleInvST = do_mleInv(ages = W,
                        sd = sd,
                        alpha = alpha,
                        K = K,
                        doMod=TRUE,
                        sd_model=sd_model)
    ,
    mleInv2 = do_mleInv(ages = W,
                       sd = sd,
                       alpha = alpha,
                       K = K,
                       method="rq2")
    ,
    mleInvP = do_mleInv(ages = W,
                        sd = sd,
                        alpha = alpha,
                        K = K,
                        method="prob")
    ,
    mleInvW = do_mleInv(ages = W,
                        sd = sd,
                        alpha = alpha,
                        K = K,
                        method="wrq")
    ,
    mleInvWS = do_mleInv(ages = W,
                        sd = sd,
                        alpha = alpha,
                        K = K,
                        method="wrq",
                        doMod=TRUE)
    ,
    mleInvWST = do_mleInv(ages = W,
                         sd = sd,
                         alpha = alpha,
                         K = K,
                         method="wrq",
                         doMod=TRUE,
                         sd_model=sd_model)
    ,
    mleInvA1 = do_mleInv(ages = W,
                        sd = sd,
                        alpha = alpha,
                        K = K,
                        method="rq",a=1)
  )
  return(estimate)
}

mle = function (W) {
  return(min(W))
}

ba_mle = function (W, K) {
  n = length(W)
  return(min(W) * (n + 1) / n - K / n)
}

strauss = function (W) {
  n = length(W)
  return((n * min(W) - max(W)) / (n - 1))
}

do_minmi = function (ages, sd, K, alpha) {
  start_time = Sys.time()
  
  results = tryCatch(
    {
      rminmi::minmi(ages, sd, K, alpha=alpha)
    },
    error = function(cond) {
      message("Something went wrong with `rminmi`:")
      message(cond)
      # Choose a return value in case of error
      return(list())
    }
  )
  runtime = calculate_tdiff(start_time, Sys.time())
  return(
    list(
      lower = as.numeric(results$theta["lower"]),
      point = as.numeric(results$theta["point"]),
      upper = as.numeric(results$theta["upper"]),
      point_runtime = runtime,
      conf_int_runtime = runtime,
      B.lower = results$B["lower"],
      B.point = results$B["point"],
      B.upper = results$B["upper"]
    )
  )
}

griwm = function(alpha,
                 dates,
                 sd,
                 K,
                 p_t,
                 .iter = 10000,
                 bias_adjusted) {
  df = data.frame(dates = dates, sd = sd)
  start_time = Sys.time()
  results = GRIWM(df, alpha, K, p_t, .iter, bias_adjusted = bias_adjusted)
  runtime = calculate_tdiff(start_time, Sys.time())
  return(
    list(
      lower = as.numeric(results$lower_ci),
      point = as.numeric(results$centroid),
      upper = as.numeric(results$upper_ci),
      point_runtime = runtime,
      conf_int_runtime = runtime
    )
  )
}

do_UNci = function (theta, ages, sd, K, alpha, wald=wald) {
  start_time = Sys.time()
  
  results = tryCatch(
    {
      getUNci(min(ages), ages, sd, K, alpha=alpha, wald=wald)
    },
    error = function(cond) {
      message("Something went wrong with `getUNci`:")
      message(cond)
      # Choose a return value in case of error
      return(list())
    }
  )
  runtime = calculate_tdiff(start_time, Sys.time())
  return(
    list(
      lower = as.numeric(results$theta["lower"]),
      point = as.numeric(results$theta["point"]),
      upper = as.numeric(results$theta["upper"]),
      point_runtime = results$se,
      conf_int_runtime = runtime,
      B.lower = results$B["lower"],
      B.point = results$B["point"],
      B.upper = results$B["upper"]
    )
  )
}

do_mleInv = function(ages, sd, K, alpha, iterMax=500, B=100, trans=trans, method="rq", doMod=FALSE, a=0, sd_model=NULL)
{

  u = matrix(runif(B*length(ages)),ncol=B)

  pt.start_time = Sys.time()
  ft.mle = getTheta(ages=ages, theta=min(ages), eps.sigma=sd, K=K, u=u)
  pt.runtime = calculate_tdiff(pt.start_time, Sys.time())

  ci.start_time = Sys.time()

  if(doMod & any(sd>0))
  {
    if(is.null(sd_model))
      eps.s = glm(sd~ages,family=Gamma("log"))
    else
      eps.s = sd_model
  }
  else
    eps.s=sd
  dat=list(W=ages,sd=sd)
  
  stepSize = max(1/sqrt(-ft.mle$hessian), IQR(ages)*0.1, na.rm=TRUE)
  thetaInits = ft.mle$par + stepSize*seq(-5,5,length=20)
  ft.lo = regInversion(dat,getT=getTh,simulateData=simFn,thetaInits=thetaInits,
                q=alpha/2,iterMax=iterMax,K=K,eps.sigma=eps.s, u=u, method=method, aMean=a, n=length(dat$W))
  ft.hi = regInversion(dat,getT=getTh,simulateData=simFn,thetaInits=thetaInits,
                       q=1-alpha/2, iterMax=iterMax, #stats=ft.lo$stats[1:length(thetaInits),], 
                       K=K,eps.sigma=eps.s, u=u, method=method, aMean=a, n=length(dat$W))
  ci.runtime = calculate_tdiff(ci.start_time, Sys.time())
  return(
    list(
      lower = ft.lo$theta,
      point = ft.mle$par,
      upper = ft.hi$theta,
      point_runtime = pt.runtime,
      conf_int_runtime = ci.runtime,
      B.lower = B,
      B.point = B,
      B.upper = B
    )
  )
}
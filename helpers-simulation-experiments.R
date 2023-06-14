
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

simulate_dataset = function (theta.true,
                             K,
                             eps.mean,
                             eps.sigma,
                             n) {
  
  eps = eps.mean
  if (any(eps.sigma != 0)) {
    eps = extraDistr::rtnorm(n, mean = eps.mean, sd = eps.sigma, a=-Inf, b=(K-theta.true))
  }
  X = runif(n, min = theta.true, max = K-eps)
  
  W = X + eps
  return(W)
}

simulate_datasets = function (config) {
  df.length = config$n.trials * config$n.error_factors
  
  dataset.df = data.frame(
    error_factor = rep(config$error_factor, length.out = df.length)
  )
  
  W = lapply(dataset.df$error_factor, 
             function(e) simulate_dataset(theta.true = config$theta.true,
                                          K = config$K,
                                          eps.mean = config$dating_error.mean,
                                          eps.sigma = e * config$fossil.sd,
                                          n = config$n.samples))
  
  dataset.df$W = W
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
                              dating_error.mean) {
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
      K = K
    )
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

do_UNci = function (theta, ages, sd, K, alpha) {
  start_time = Sys.time()
  
  results = tryCatch(
    {
      getUNci(0.8*theta, ages, sd, K, alpha=alpha)
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
      point_runtime = runtime,
      conf_int_runtime = runtime,
      B.lower = results$B["lower"],
      B.point = results$B["point"],
      B.upper = results$B["upper"]
    )
  )
}

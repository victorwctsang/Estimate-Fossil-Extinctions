# Simulation Experiments - minmi vs UNloglik
# Get estimates using simulated datasets

library(logger)
library(rminmi) # version 0.1.1 or greater
source('helpers-simulation-experiments.R')
source('UNloglik.R')

log_info("Loading synthetic dataset and configuration")

load("data/synthetic-data.RData")
attach(synthetic.data.config)

alpha = 0.05

methods.point_estimates = c(
  # "Strauss"
  # ,"MLE"
  # ,"BA-MLE"
)
methods.conf_int = c(
  "MINMI"
  , "UNci"
  , "UNwald"
  # ,"GRIWM"
  # ,"GRIWM-corrected"
)

RESULTS_PATH = paste0("data/simResults-", format(Sys.Date(), "%Y%m%d"), ".RData")

results = data.frame(
  error_factor=double(),
  method=factor(),
  lower=double(),
  point=double(),
  upper=double(),
  point_runtime=double(),
  conf_int_runtime=double(),
  B.lower=double(),
  B.point=double(),
  B.upper=double()
)

pilot.dates = datasets[1, "W"][[1]]
A = 0.1 * (mean(fossil.sd))

############################################################
# Run Trials
############################################################
log_info("Performing Simulations")

start_time = Sys.time()
#for (i in 1:100) {
for (i in 1:nrow(datasets)) {
  log_info(sprintf("Dataset: %i/%i", i, nrow(datasets)))
  iter = datasets[i, ]
  W = as.numeric(iter$W[[1]])
  sd = as.numeric(iter$error_factor * fossil.sd)
  
  log_info("Getting point estimates")
  for (method in methods.point_estimates) {
    estimation = estimate_extinction(
      W = W,
      sd = sd,
      method = method,
      K = K,
      dating_error.mean = dating_error.mean
    )
    log_info(sprintf("Time taken for %s: %.02f seconds", method, estimation$point_runtime))
    
    results = tibble::add_row(
      results,
      error_factor = iter$error_factor,
      method = method,
      point = estimation$point,
      point_runtime = estimation$point_runtime,
      conf_int_runtime = estimation$conf_int_runtime
    )
  }
  
  log_info("Getting confidence intervals")
  for (method in methods.conf_int) {
    estimation = estimate_conf_int(
      W = W,
      sd = sd,
      method = method,
      alpha = alpha,
      K = K,
      dating_error.mean = dating_error.mean
    )
    log_info(sprintf("Time taken for %s: %.02f seconds", method, estimation$conf_int_runtime))
    results = tibble::add_row(
      results,
      error_factor = iter$error_factor,
      method = method,
      lower = estimation$lower,
      point = estimation$point,
      upper = estimation$upper,
      point_runtime = estimation$point_runtime,
      conf_int_runtime = estimation$conf_int_runtime,
      B.lower = estimation$B.lower,
      B.point = estimation$B.point,
      B.upper = estimation$B.upper
    )
  }
}
sims.time = Sys.time() - start_time
log_info(sims.time)

## Save results
log_info(sprintf("Saving results to %s", RESULTS_PATH))
save(results, file = RESULTS_PATH)
View(results)

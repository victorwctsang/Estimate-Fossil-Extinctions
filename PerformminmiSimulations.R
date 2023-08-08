# Simulation Experiments - minmi vs UNloglik
# Get estimates using simulated datasets

#nSim=5000
library(logger)
library(rminmi) # version 0.1.1 or greater
source('helpers-simulation-experiments.R')

#log_info("Loading synthetic dataset and configuration")

#load("data/synthetic-data-12-20230702.RData")

PerformminmiSimulations = function(whichSims,datasets,synthetic.data.config,alpha=0.05,
                                   methods.point_estimates = c(
  # "Strauss"
  # ,"MLE"
  # ,"BA-MLE"
),
methods.conf_int = c(
  "MINMI"
  ,"mlereginv"
  ,"reginvUNci"
  ,"reginvUNwald"
#  , "UNciA"
#  , "UNwaldA"
#  , "mleInvA"
#  , "mleInvAW"
  # ,"GRIWM"
  # ,"GRIWM-corrected"
)
)
{
  
RESULTS_PATH = paste0("data/simResults-", synthetic.data.config$n.samples, "-", format(Sys.Date(), "%Y%m%d"), "_", max(whichSims),".RData")
print(RESULTS_PATH)

results = data.frame(
  which_sim=double(),
  n.samples=double(),
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
A = 0.1 * (mean(synthetic.data.config$fossil.sd))

############################################################
# Run Trials
############################################################
log_info("Performing Simulations")

start_time = Sys.time()
for (i in whichSims) {
#for (i in 1:nrow(datasets)) {
  log_info(sprintf("Dataset: %i/%i", i, nrow(datasets)))
  iter = datasets[i, ]
  W = as.numeric(iter$W[[1]])
  sd = as.numeric(iter$error_factor * synthetic.data.config$fossil.sd)
  
  log_info("Getting point estimates")
  for (method in methods.point_estimates) {
    estimation = estimate_extinction(
      W = W,
      sd = sd,
      method = method,
      K = synthetic.data.config$K,
      dating_error.mean = synthetic.data.config$dating_error.mean
    )
    log_info(sprintf("Time taken for %s: %.02f seconds", method, estimation$point_runtime))
    
    results = tibble::add_row(
      results,
      which_sim=i,
      n.samples=synthetic.data.config$n.samples,
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
      K = synthetic.data.config$K,
      dating_error.mean = synthetic.data.config$dating_error.mean
    )
    log_info(sprintf("Time taken for %s: %.02f seconds", method, estimation$conf_int_runtime))
    results = tibble::add_row(
      results,
      which_sim=i,
      n.samples=synthetic.data.config$n.samples,
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
#View(results)
invisible(results)
}

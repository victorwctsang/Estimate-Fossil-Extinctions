# Simulation Experiments - minmi vs UNloglik
# Get estimates using simulated datasets

library(logger)
library(rminmi) # version 0.1.1 or greater
source('helpers-simulation-experiments.R')
source('UNloglik.R')

log_info("Loading synthetic dataset and configuration")

#load("data/synthetic-data.RData")
load("data/synthetic-data-100-20230616.RData")
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
<<<<<<< Updated upstream
  # ,"GRIWM"
=======
  , "mleInv"
#  , "mleInv2"
#  , "mleInvS"
, "mleInvST"
, "mleInvW"
  , "mleInvWS"
, "mleInvWST"
#, "mleInvA1"
# ,"GRIWM"
>>>>>>> Stashed changes
  # ,"GRIWM-corrected"
)

n.samples = length(datasets[1,]$W[[1]])
RESULTS_PATH = paste0("data/simResults-", n.samples, "-", format(Sys.Date(), "%Y%m%d"), ".RData")

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

<<<<<<< Updated upstream
pilot.dates = datasets[1, "W"][[1]]
A = 0.1 * (mean(fossil.sd))

=======
>>>>>>> Stashed changes
############################################################
# Run Trials
############################################################
log_info("Performing Simulations")

start_time = Sys.time()
<<<<<<< Updated upstream
#for (i in 1:100) {
for (i in 1:nrow(datasets)) {
  log_info(sprintf("Dataset: %i/%i", i, nrow(datasets)))
  iter = datasets[i, ]
  W = as.numeric(iter$W[[1]])
  sd = as.numeric(iter$error_factor * fossil.sd)
=======
for (iSim in whichSims) {
#for (iSim in 1:length(datasets$error_factor)) {
  log_info(sprintf("Dataset: %i/%i", iSim, length(datasets$error_factor)))
  W = as.numeric(datasets$W[,iSim])
  sd = as.numeric(datasets$eps[,iSim])

  # set sdCurrent as sds or sdModel, depending on inputs
  sdCurrent = synthetic.data.config$fossil.sd
  if(datasets$error_factor[iSim]==0) sdCurrent = rep(0,synthetic.data.config$n.samples)
  if(inherits(sdCurrent,"glm"))
    sdCurrent$coefficients[1] = synthetic.data.config$fossil.sd$coefficients[1] + log(datasets$error_factor[iSim])
>>>>>>> Stashed changes
  
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
<<<<<<< Updated upstream
      error_factor = iter$error_factor,
=======
      which_sim=iSim,
      n.samples=synthetic.data.config$n.samples,
      error_factor = datasets$error_factor[iSim],
>>>>>>> Stashed changes
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
<<<<<<< Updated upstream
      K = K,
      dating_error.mean = dating_error.mean
=======
      K = synthetic.data.config$K,
      dating_error.mean = synthetic.data.config$dating_error.mean,
      sd_model = sdCurrent
>>>>>>> Stashed changes
    )
    log_info(sprintf("Time taken for %s: %.02f seconds", method, estimation$conf_int_runtime))
    results = tibble::add_row(
      results,
<<<<<<< Updated upstream
      error_factor = iter$error_factor,
=======
      which_sim=iSim,
      n.samples=synthetic.data.config$n.samples,
      error_factor = datasets$error_factor[iSim],
>>>>>>> Stashed changes
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

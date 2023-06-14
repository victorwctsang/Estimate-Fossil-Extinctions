load("synthetic-data.RData")

# Select the problematic dataset (this is only one of many)
idx <- 310
W = as.numeric(datasets[idx, ]$W[[1]])
error_factor = as.numeric(datasets[idx, ]$error_factor[[1]])
# NB: We use the same measurement error standard errors for all the datasets, but apply different error multipliers.
s = synthetic.data.config$fossil.sd * error_factor

rminmi::minmi(ages = W,
              sd = s,
              K = synthetic.data.config$K,
              alpha = 0.05
              # Don't give any additional tuning parameters for the Monte Carlo estimates
              # B = NA, A = NA, .B_init = NA
)


thetas = seq(10000,25000,length=100)
eeqns = rep(NA,length(thetas))
for (iTheta in 1:length(thetas))
  eeqns[iTheta] = estimating_eqn(thetas[iTheta], q, K, u, m, n, eps.sigma)
plot(thetas,eeqns,type="l")
abline(a=0,b=0,col="red")

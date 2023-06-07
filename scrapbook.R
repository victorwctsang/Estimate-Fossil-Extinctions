
source('helpers-simulation-experiments.R')

load("data/synthetic-data.RData")
attach(synthetic.data.config)

alpha = 0.05

problem_idx <- 310
problem_row <- datasets[problem_idx, ]
 
estimate <- rminmi::minmi(
  ages = as.numeric(problem_row$W[[1]]),
  sd = as.numeric(problem_row$error_factor * fossil.sd),
  K,
  alpha
)

print(estimate)
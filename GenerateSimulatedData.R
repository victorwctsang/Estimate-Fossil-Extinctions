# Simulation Experiments
# Generate synthetic datasets for simulation experiments


#### Setup
library(tidyverse)
library(readxl)

source("helpers-simulation-experiments.R")

options(nwarnings = 10000)

# Configure synthetic data generation
synthetic.data.config <- list(
  seed = 60,
  theta.true = 15000,
  K = 20000,
  n.trials = 1000,
  n.samples = 60,
  error_factors = c(0, 0.5, 1, 2, 4),
  dating_error.mean = 0
)

SYNTH_DATA_PATH = paste0("data/synthetic-data-", synthetic.data.config$n.samples, "-", format(Sys.Date(), "%Y%m%d"), ".RData")

#set.seed(12) # for n=12
#set.seed(24) # for n=24
#set.seed(36) # for n=36
#set.seed(48) # for n=48
set.seed(60) # for n=60

synthetic.data.config$n.error_factors <- length(synthetic.data.config$error_factors)

mammoth_data = read_excel(path='data/fossildata.xlsx',
                          sheet="Mammoths Eurasian", 
                          range="M3:N204", 
                          col_names=c("age", "sd"), 
                          col_types=c('numeric', 'numeric'))

if(synthetic.data.config$n.samples<=length(mammoth_data$sd))
  synthetic.data.config$fossil.sd <- mammoth_data$sd[1:synthetic.data.config$n.samples]
if(synthetic.data.config$n.samples>length(mammoth_data$sd))
  synthetic.data.config$fossil.sd <- sample(mammoth_data$sd,synthetic.data.config$n.samples,replace=TRUE)

attach(synthetic.data.config)

# Generate synthetic data
datasets <- simulate_datasets(synthetic.data.config)

# Export config and datasets
save(synthetic.data.config, datasets, file = SYNTH_DATA_PATH)

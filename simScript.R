nSim=2000

start.all.sims=Sys.time()

load("data/synthetic-data-12-20230616.RData")
attach(synthetic.data.config)
source("PerformminmiSimulations.R")

load("data/synthetic-data-25-20230616.RData")
attach(synthetic.data.config)
source("PerformminmiSimulations.R")

load("data/synthetic-data-50-20230616.RData")
attach(synthetic.data.config)
source("PerformminmiSimulations.R")

end.all.sims = Sys.time()
tot.all.sims = end.all.sims-start.all.sims
print("Total sim time:")
print(tot.all.sims)


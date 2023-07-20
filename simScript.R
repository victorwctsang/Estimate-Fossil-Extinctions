nSim=5000

start.all.sims=Sys.time()

load("data/synthetic-data-48-20230630.RData")
source("PerformminmiSimulations.R")

load("data/synthetic-data-25-20230616.RData")
source("PerformminmiSimulations.R")

load("data/synthetic-data-50-20230616.RData")
source("PerformminmiSimulations.R")

end.all.sims = Sys.time()
tot.all.sims = end.all.sims-start.all.sims
print("Total sim time:")
print(tot.all.sims)

nSim=5
load("data/synthetic-data-60-20230704.RData")

# things to do:
# repeat sims of 48-20230630 which had good results, see if now bad
# change the thing where I take lower and use that towards upper, maybe increase sim limit or something
# try a non-linear fit (method=rq2)
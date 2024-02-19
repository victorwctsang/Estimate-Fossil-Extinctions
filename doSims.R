source("PerformminmiSimulations.R")

do.sim = function(whichSims,nSamp=c(6,9,12,18,24,36,48,72,96),date)
{
  for(iSamp in nSamp)
  {
    filename=paste0("data/synthetic-data-",iSamp,"-",date,".RData")
    print(filename)
    load(filename)
    PerformminmiSimulations(whichSims,datasets,synthetic.data.config)
  }
}
whichSims=8001:8005
system.time(do.sim(whichSims, nSamp=c(6), date="20230901"))

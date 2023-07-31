source("PerformminmiSimulations.R")


do.sim = function(whichSims,nSamp=c(12,24,36,48,60),date)
{
  for(iSamp in nSamp)
  {
    filename=paste0("data/synthetic-data-",iSamp,"-",date,".RData")
    print(filename)
    load(filename)
    PerformminmiSimulations(whichSims,datasets,synthetic.data.config)
  }
}
#whichSims=6:10
#system.time(do.sim(whichSims, date="20230725"))

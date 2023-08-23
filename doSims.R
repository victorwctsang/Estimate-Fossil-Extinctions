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
#whichSims=5901:5910
#system.time(do.sim(whichSims, nSamp=48,date="20230823"))

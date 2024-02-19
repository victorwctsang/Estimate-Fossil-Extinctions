# Take results from cluster and stick them all together in a big fat file yay

date="20230926"
for(iSamp in c(6,9,12,18,24,36,48,72,96))
#  for(iSamp in c(96))
#    for(iSamp in c(18,24,36,48))
    {
  fileNameAll = paste0("data/simResults-",iSamp,"-",date,".RData")
  fileName = paste0("data/simResults-",iSamp,"-",date,"_30.RData")
  load(fileName)
  resultsAll = results
  for (iRun in 1:299)
  {
    fileNamei = paste0("data/simResults-",iSamp,"-",date,"_",30*(iRun+1),".RData")
    tf=try(load(fileNamei))
    if(inherits(tf,"try-error")==FALSE)
      resultsAll = rbind(resultsAll,results)
  }
  results=resultsAll
  save(results,file=fileNameAll)
}

# Take results from cluster and stick them all together in a big fat file yay

for(iSamp in c(12,24,36,48,60))
{
  fileNameAll = paste0("data/simResults-",iSamp,"-20230822.RData")
  fileName = paste0("data/simResults-",iSamp,"-20230822_100.RData")
  load(fileName)
  resultsAll = results
  for (iRun in 1:59)
  {
    fileNamei = paste0("data/simResults-",iSamp,"-20230822_",100*(iRun+1),".RData")
    tf=try(load(fileNamei))
    if(inherits(tf,"try-error")==FALSE)
      resultsAll = rbind(resultsAll,results)
  }
  results=resultsAll
  save(results,file=fileNameAll)
}

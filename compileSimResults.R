# Take results from cluster and stick them all together in a big fat file yay

for(iSamp in c(12,24,36,48,60))
{
  fileName = paste0("data/simResults-",iSamp,"-20230714_40.RData")
  load(fileName)
  resultsAll = results
  for (iRun in 1:124)
  {
    fileNamei = paste0("data/simResults-",iSamp,"-20230714_",40*(iRun+1),".RData")
    tf=try(load(fileNamei))
    if(inherits(tf,"try-error")==FALSE)
      resultsAll = rbind(resultsAll,results)
  }
  fileNameAll = paste0("data/simResults-",iSamp,"-20230714.RData")
  results=resultsAll
  save(results,file=fileNameAll)
}

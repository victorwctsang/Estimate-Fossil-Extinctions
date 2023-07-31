# Take results from cluster and stick them all together in a big fat file yay

fileName = paste0("data/sdFit/simResults-12-20230726_50.RData")
load(fileName)
for(iSamp in c(12,24,36,48,60))
{
  resultsAll = results[0,]
  for (iRun in 1:100)
  {
    fileNamei = paste0("data/sdFit/simResults-",iSamp,"-20230726_",50*iRun,".RData")
    tf=try(load(fileNamei))
    if(inherits(tf,"try-error")==FALSE)
      resultsAll = rbind(resultsAll,results)
  }
  fileNameAll = paste0("data/simResults-sdFit-",iSamp,"-20230726.RData")
  results=resultsAll
  save(results,file=fileNameAll)
}

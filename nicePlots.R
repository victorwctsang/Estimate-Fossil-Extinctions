library(dplyr)

load("data/synthetic-data-6-20230901.RData")
attach(synthetic.data.config)

### do a plot of bias at sd=500, as a function of interval length

RESULTS_PATH <- 'data/simResults-24-20230926.RData'
load(RESULTS_PATH)

performance.point <- results %>%
  filter(!is.na(point),method%in%c("Strauss","UTbias")) %>%
  group_by(error_factor, method) %>%
  summarise(MSE_000 = mean((point - theta.true)^2,na.rm=TRUE)/1000,
            bias = mean(point,na.rm=TRUE)-theta.true,
            variance_000 = var(point,na.rm=TRUE)/1000,
            avg_runtime = round(mean(point_runtime,na.rm=TRUE), 5))

performance.rel = performance.point %>% filter(error_factor>0.5) %>% mutate(biasAdj = bias*5/error_factor + 10000, interval=50000/error_factor)

png(filename="biasPlot.png",width=3000,height=2500,res=500)
{
  par(mar=c(3,4,0.5,0.5),mgp=c(1.75,0.75,0),las=1)
  plot(biasAdj~interval,data=performance.rel%>%filter(method=="Strauss"),type="l",ylim=range(performance.rel$biasAdj),log="x",xaxt="n",col="blue",ylab="",xlab="")
  mtext(expression(theta),2,las=1,line=3,cex=1.5)
  mtext(expression(K-theta),1,las=1,line=2,cex=1.5)
  axis(1,c(1000,2000,5000,10000,20000,50000))
  #points(biasAdj~interval,data=performance.rel%>%filter(method=="UNci"),type="l",col="darkgreen")
  #points(biasAdj~interval,data=performance.rel%>%filter(method=="UTci"),type="l",col="red")
  abline(h=10000,col="grey80")
  points(biasAdj~interval,data=performance.rel%>%filter(method=="UTbias"),type="l",col="purple")
  text(performance.rel[performance.rel$method=="Strauss","interval"][5,1],performance.rel[performance.rel$method=="Strauss","biasAdj"][5,1],"Strauss",pos=4,col="blue")
  text(performance.rel[performance.rel$method=="UTbias","interval"][5,1],performance.rel[performance.rel$method=="UTbias","biasAdj"][5,1],"bias-corrected",pos=1,col="purple")
}
dev.off()

### do a plot of CI coverage against n at error=2
error_fac_to_plot = 2 #change this value to look at plots for a different error_factor

whichDate = "20230926" #change this value to date of sim run

RESULTS_PATH <- paste0('data/simResults-6-',whichDate,'.RData')
load(RESULTS_PATH)

head(results)
all_results=results
all_results$n.samples = synthetic.data.config$n.samples
all_results=all_results[0,]

n.samples=c(6,9,12,18,24,36,48,72,96)
for (iSample in 1:length(n.samples))
{
  RESULTS_PATH <- paste0("data/simResults-",n.samples[iSample],"-",whichDate,".RData")
  load(RESULTS_PATH)
  all_results=tibble::add_row(
    all_results,
    error_factor = results$error_factor,
    method=results$method,
    lower=results$lower,
    point=results$point,
    upper=results$upper,
    point_runtime=results$point_runtime,
    conf_int_runtime=results$conf_int_runtime,
    B.lower=results$B.lower,
    B.point=results$B.point,
    B.upper=results$B.upper,
    n.samples=n.samples[iSample]
  )
}
performance.CI <- all_results %>%
  filter(!is.na(conf_int_runtime), method%in%c("UTci","reginvUT","GRIWM"), error_factor == error_fac_to_plot) %>%
  mutate(width = upper - lower,
         contains_theta = ifelse(theta.true > lower & theta.true < upper, 1, 0)) %>%
  group_by(n.samples, method) %>%
  summarise(Coverage = mean(contains_theta, na.rm=TRUE) * 100,
            `Average Width` = mean(width, na.rm=TRUE),
            `Trimmed Mean Width` = mean(width, trim=0.05, na.rm=TRUE),
            `Average Runtime` = mean(conf_int_runtime, na.rm=TRUE)) %>%
  ungroup() %>%
  arrange(method, n.samples)


png(filename="CIcover.png",width=3000,height=2500,res=500)
{
  par(mar=c(3,3,0.5,0.5),mgp=c(1.75,0.75,0),las=1)
  plot(Coverage~n.samples,data=performance.CI%>%filter(method=="GRIWM"),type="l",col="blue",ylim=range(performance.CI$Coverage),ylab="",xlab="")
  #points(biasAdj~interval,data=performance.rel%>%filter(method=="UNci"),type="l",col="darkgreen")
  #points(biasAdj~interval,data=performance.rel%>%filter(method=="UTci"),type="l",col="red")
  mtext("Coverage Probability (%)",2,line=1.75,las=0,cex=1.25)
  mtext("Sample size",1,line=1.75,las=0,cex=1.25)
  abline(h=10000,col="grey80")
  xRange=range(performance.CI$n.samples)
  yRange = qbinom(c(0.025,0.975),1000,0.95)/10
  polygon(xRange[c(1,2,2,1)],yRange[c(1,1,2,2)],col="grey95",border=NA)
    points(Coverage~n.samples,data=performance.CI%>%filter(method=="reginvUT"),type="l",col="purple")
  text(performance.CI[performance.CI$method=="GRIWM","n.samples"][5,1],performance.CI[performance.CI$method=="GRIWM","Coverage"][5,1],"GRIWM",pos=4,col="blue")
  text(performance.CI[performance.CI$method=="reginvUT","n.samples"][5,1],performance.CI[performance.CI$method=="reginvUT","Coverage"][5,1],"est_cutt",pos=1,col="purple")
}
dev.off()

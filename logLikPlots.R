
load("data/synthetic-data-1000-20230630.RData")

sds=c(0,125,250,375,500)

cols=colorRampPalette(c("darkblue","pink"))(length(sds))

i=1
iter = datasets[i, ]
W = as.numeric(iter$W[[1]])
#sd = as.numeric(iter$error_factor * fossil.sd)
sd=rep(sds[i],length(W))

source("UNloglik.R")
thetas=seq(13000,16000,length=100)
B=100
u = matrix(runif(B*length(W)),ncol=B)
lls = scores = list()
lls[[1]]=scores[[1]] = matrix(NA,length(W),length(thetas))
for(iTheta in 1:length(thetas))
{
  for(iData in 1:length(W))
  {
    lls[[1]][iData,iTheta] = UNaltloglik(thetas[iTheta],W[iData],sd[iData],K=synthetic.data.config$K)
    scores[[1]][iData,iTheta] = UNaltscore(thetas[iTheta],W[iData],sd[iData],K=synthetic.data.config$K)
  }
}

for(i in 2:5)
{
#  iter = datasets[i, ]
#  W = as.numeric(iter$W[[1]])
  sd=rep(sds[i],length(W))
  #  sd = as.numeric(iter$error_factor * fossil.sd)

  thetas=seq(13000,16000,length=100)
  B=100
  u = matrix(runif(B*length(W)),ncol=B)
  lls[[i]] = scores[[i]] = matrix(NA,length(W),length(thetas))
  for(iTheta in 1:length(thetas))
  {
    for(iData in 1:length(W))
    {
      lls[[i]][iData,iTheta] = UNaltloglik(thetas[iTheta],W[iData],sd[iData],K=synthetic.data.config$K)
      scores[[i]][iData,iTheta] = UNaltscore(thetas[iTheta],W[iData],sd[iData],K=synthetic.data.config$K)
    }
  }
}



library(latex2exp)
label=TeX("$\\frac{\\sigma}{\\mu}exp(-\\pi \\phi)$")


whichTheta = which(thetas== synthetic.data.config$theta.true)

#par(mfrow=c(3,2),mar=c(3,3,1,1),mgp=c(1.75,0.75,0))
#for(i in 1:5)
#  boxplot(lls[[i]][,whichTheta])

sc=data.frame(cbind(scores[[1]][,whichTheta],scores[[2]][,whichTheta],scores[[3]][,whichTheta],scores[[4]][,whichTheta],scores[[5]][,whichTheta]))
names(sc) = paste0(sds/50,"%")

# get the plots

png("likelihoods.png",width=3000,height=3000,res=500)
{
  par(mfrow=c(2,1),mar=c(3,3,1,1),mgp=c(1.5,0.75,0))

llSumm = cbind(apply(lls[[1]],2,sum),apply(lls[[2]],2,sum), apply(lls[[3]],2,sum), apply(lls[[4]],2,sum), apply(lls[[5]],2,sum) )
plot(thetas,llSumm[,1],type="l",xlab=TeX("$\\theta$"),ylab=TeX("$l(\\theta,y)$"),col=cols[1],xaxt="n",yaxt="n",ylim=range(llSumm[llSumm[,1]>-Inf,1]),cex.lab=1.25,lwd=1.25)
for(i in 2:5)
  lines(thetas,apply(lls[[i]],2,sum),type="l",col=cols[i],lwd=1.25)
xLab=seq(13000,16000,by=500)
yLab=seq(min(llSumm[llSumm[,1]>-Inf,1]),max(llSumm[,1]),length=4)
yLab = 2*max(llSumm[,1])
axis(2,at=yLab,labels=rep("",length(yLab)))
legend("topleft",legend=names(sc),lty=1,col=cols)
W = as.numeric(datasets$W[[1]])
abline(v=min(W),col="grey70",lty=1,lwd=0.5)
xLab=min(W)
axis(1,at=xLab,labels=rep("min(W)",length(xLab)),col.axis="grey70")


boxplot(sc,col=cols,xlab="",ylab=TeX("$l\\prime(\\theta,y_i)$"),outcol=cols,yaxt="n",cex.lab=1.25)
abline(h=0,col="grey70",lty=1,lwd=0.5)
boxplot(sc,col=cols,xlab="",ylab=TeX("$l\\prime(\\theta,y_i)$"),outcol=cols,yaxt="n",add=TRUE)
mtext(expression(sigma),1,line=1.75)
step=diff(range(sc))*4
y2Lab=step*c(-2,-1,0,1,2,3)
axis(2,at=y2Lab,labels=c("","",0,"","",""),las=1,col.axis="grey70")
}
dev.off()

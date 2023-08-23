iSamp=600
date="20230822"
filename=paste0("data/synthetic-data-",iSamp,"-",date,".RData")
K = synthetic.data.config$K
df = synthetic.data.config$df

load(filename)
i=108
iter=datasets[i,]
W = as.numeric(iter$W[[1]])
sd = as.numeric(iter$error_factor * synthetic.data.config$fossil.sd)

library(reginv)
mle_fossil(W,sd,K,df)
reginv_fossil(W,sd,K,df,q=0.975) #every now and then returns 20000 (!)

ft.mle = mle_fossil(W, sd=sd, K, df, alpha=NULL)
stepSize = ifelse( any(sd==0), IQR(W)*0.1, ft.mle$se )
paramInits = ft.mle$theta + stepSize*seq(-5,5,length=20)

df=4
ft=reginv(W, getT=reginv:::getThMLE, simulateData=reginv:::simFn_fossil, paramInits=paramInits,q=0.975,K=K,sd=sd,df=4,n=length(W))
plot(T~theta,data=ft$stats)

Wi=reginv:::simFn_fossil(theta=15000,K,sd,df=4,n=length(W))
reginv:::getThMLE(Wi,K,sd,df=4,n=length(Wi))

qfossil(0.1,45000,50000,10000,df)

ts=seq(10000,50000,length=100)
plot(ts,pfossil(ts,45000,50000,1000,2),type="l")

pfossil(10000,45000,50000,1000,4)

qt(0.001,2)

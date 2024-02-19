# testing estimators of df
require(reginv)
df=exp(seq(log(2.8),log(50),length=100))

ages=rcutt(100,10000,25000,1000,df=2100)
th = mle_cutt(ages,sd=1000,K=25000,df=2100); print(th)

lls = profs= rep(NA,length=length(df))
for(iDF in 1:length(df))
{
  lls[iDF] = reginv:::cutt_LogLikT(1/df[iDF],ages=ages,sd=1000,theta=th$theta,K=25000) 
  thI = mle_cutt(ages,sd=1000,K=25000,df=df[iDF])
  profs[iDF] = reginv:::cutt_LogLikT(1/df[iDF],ages=ages,sd=1000,theta=thI$theta,K=25000) 
  
}
plot(lls~df,type="l",ylim=range(c(lls,profs)))
lines(profs~df,col="red")

thJoint = mle_cutt(ages,sd=1000,K=25000,df=NULL); print(thJoint); print(thJoint$df)

n=10000
ps=(1:n)/(n+1)
qs=reginv::qcutt(ps,10000,20000,3000,df=3)
plot(qs~ps,type="l")

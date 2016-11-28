#Question 4#

rm(list=ls())
library(MASS)
#input simdathw2#
SimDatHW2=function(p.alt=10, # the number of differentially expressed features ... 
                   p.null=90, # the number of non-differentially expressed features 
                   n=20, # the number of samples in each of the treatment and control groups
                   rho.alt=0.2, # correlation of the alt hyp variables
                   rho.null=0.1, # correlation of the null hyp variables
                   delta=2,# the mean of the p.alt features in the "treatment" group
                   sdC=0 # the variance of the sample specific centering error
){
  p=p.alt+p.null
  Sigma=array(rep(0,p^2),dim=c(p,p))
  Sigma[1:p.alt,1:p.alt]=rho.alt
  Sigma[(p.alt+(1:p.null)),(p.alt+(1:p.null))]=rho.null
  diag(Sigma)=1
  Xc=mvrnorm(n,mu=rep(0,p),Sigma=Sigma)
  Xt=mvrnorm(n,mu=c(rep(delta,p.alt),rep(0,p.null)),Sigma=Sigma)
  x=t(rbind(Xc,Xt))
  colnames(x)=c(paste("cntrl",1:n,sep=""),paste("trt",1:n,sep=""))
  rownames(x)=c(paste("DEgene",1:(p.alt)),paste("nonDEgene",1:p.null,sep=""))
  if(sdC>0){
    CentError=rnorm(2*n,sd=sdC)
    for(j in 1:(2*n)) x[,j]=x[,j]+CentError[j]
  }
  return(x)
}


#First, we generate data without and with centering error from SimDatHW2 
#and named them X1 and X2a respectively
X1= SimDatHW2(p.alt=25,p.null=50,n=20,rho.alt=.3,rho.null=0.0,delta=2,sdC=0)
X2a= SimDatHW2(p.alt=25,p.null=50,n=20,rho.alt=.3,rho.null=0.0,delta=2,sdC=1)
#Create another data samw as X2a but re-centered about each sample (i.e., column) mean. 
X2b= X2a 
for(j in 1:40) X2b[,j]=X2b[,j]-mean(X2b[,j])


#Question 4.1
#Since all the null targets were generated with zero correlation
#Now we check all the pairwise correlations between null targets 
#and see if their pairwise correlation distributions are centered about zero or not
#The number of pairwise correlation in null target now should be 
choose(50,2)
#1225

#We generate the pairwise correlation of null targets in each group
#and also plot a histogram to see their distribution
#For Data X1
RhoNullX1=cor(t(X1[26:75,]))
hist(RhoNullX1[upper.tri(RhoNullX1)])
#For Data X2a
RhoNullX2a=cor(t(X2a[26:75,]))
hist(RhoNullX2a[upper.tri(RhoNullX2a)])
#For Data X2b
RhoNullX2b=cor(t(X2b[26:75,]))
hist(RhoNullX2b[upper.tri(RhoNullX2b)])


#Question 4.2
#Check if all the nulls were generated such that the p-values 
#for null targets have an uniform distribution
myTX1=function(i) t.test(X1[i,1:20],X1[i,21:40],var.equal=T)$p.value
pvalueX1=mapply(myTX1,26:dim(X1)[1])
plot(ecdf(pvalueX1));abline(0,1,col=2)
ks.test(pvalueX1,'punif')
#
myTX2a=function(i) t.test(X2a[i,1:20],X2a[i,21:40],var.equal=T)$p.value
pvalueX2a=mapply(myTX2a,26:dim(X2a)[1])
plot(ecdf(pvalueX2a));abline(0,1,col=2)
ks.test(pvalueX2a,'punif')
#
myTX2b=function(i) t.test(X2b[i,1:20],X2b[i,21:40],var.equal=T)$p.value
pvalueX2b=mapply(myTX2b,26:dim(X2b)[1])
plot(ecdf(pvalueX2b));abline(0,1,col=2)
ks.test(pvalueX2b,'punif')


# let's look at the correlation plot of null data without centering error
plot(RhoNullX1[1,],RhoNullX1[2,],main=paste("Sample Correlation =",round(cor(RhoNullX1[1,],RhoNullX1[2,]),4)))

# Now,let's look at the plot of data with centering error added to each column 
plot(RhoNullX2a[1,],RhoNullX2a[2,],main=paste("Sample Correlation =",round(cor(RhoNullX2a[1,],RhoNullX2a[2,]),4)))

xylim=range(c(RhoNullX1[1,],RhoNullX2a[1,],RhoNullX1[2,],RhoNullX2a[2,])) 

plot(c(RhoNullX1[1,],RhoNullX2a[1,]),c(RhoNullX1[2,],RhoNullX2a[2,]),
     xlim=xylim,ylim=xylim,type="n",xlab="Target 1", ylab="Target 2") 
abline(0,1,lty=2,lwd=3,col="gray67") 
# now,let's add original points 
points(RhoNullX1[1,],RhoNullX1[2,],col="gray47") 
# now, let's add the points with centering error and connect with lines 
for(j in 1:25){ 
        lines(c(RhoNullX1[1,j],RhoNullX2a[1,j]),c(RhoNullX1[2,j],RhoNullX2a[2,j]),col="red")
        lines(c(RhoNullX1[1,j],RhoNullX2b[1,j]),c(RhoNullX1[2,j],RhoNullX2b[2,j]),col="blue") 
        points(RhoNullX2a[1,],RhoNullX2a[2,],pch=10, col="red") 
        points(RhoNullX2b[1,],RhoNullX2b[2,],pch=10, col="blue") 
}


plot(c(RhoNullX1[1,],RhoNullX2b[1,]),c(RhoNullX1[2,],RhoNullX2b[2,]),
     xlim=xylim,ylim=xylim,type="n",xlab="Target 1", ylab="Target 2") 
abline(0,1,lty=2,lwd=3,col="gray67") 
# now,let's add original points 
points(RhoNullX1[1,],RhoNullX1[2,],col="gray47") 
# now, let's add the points with centering error and connect with lines 
for(j in 1:25){ 
        lines(c(RhoNullX1[1,j],RhoNullX2b[1,j]),c(RhoNullX1[2,j],RhoNullX2b[2,j]),col=2) 
        points(RhoNullX2b[1,],RhoNullX2b[2,],pch=10) 
}





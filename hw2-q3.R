#hw2-q3
#install.packages("pwr")
library(pwr)

# 200 targets across two patient populations (group A and B; n=10 ea.)
# 1 target is DE, 199 not DE
# target exprs values dist N(0,1) for 199
# target exprs value  mean=d and sd = 1 for 1 
# null targets have multivariate distn with a pairwise correlation of rho for all null target pairings

## ANALYSIS PLAN:
## perform two sample equal variance t-tests and then use maxT to adjust the pvalue
## Question: why might we prefer the maxT adjustment to that of Bonferroni


# null hypothesis is no assoication between the expression levels and the responses or covariates
# since several thousand genes simulantenously ... extreme multiple testing problem
#two types of errors can occurs: false positive (Type I error) -- when a gene is declared to be a DEG and it is not
# OR false negative (type II error) -- when test fails to identify truly DEG
#special probelms arise from the multiplicity aspect include defining an appropriate type I error rate, and devising
#multiple testing procedures which control this error rate and incorporate the joint distn of the test statistics
# 

#POWER -- probability of correctly rejecting Ho, and avoiding a type II error -- 0.80 is standard

#Effect Size -- how large the effect of one variable on another variable needs to be. The degree to which Ho is false

#Significance level (alpha) - the probability of rejecting a true Ho which means making a type I error
#0.05 is the minimum and most used

#Power = 1-beta
#Power --- likelihood of rejecting Ho when you should
#ability of test to effect an effect if an effect really exists
#if power is lower ... risk of type II error
#rule of thumb: power increases as sample size increases 

#Cohens guidelines ... SD = 0.2 (small), SD = 0.5 (medium effect), SD = 0.8 (large effect)
#effect size -- the degree to which Ho is false. It measures how large the effect of one variable needs to be another needs to be
#type I and type II are inversely related .. if you decrease one you increase the other
#type I error is more serious ... if you go to smaller 0.01 ... it will increase type II error

#alpha values are determined before the study
#if P < alpha ... results are statistically significant


library(multtest)
library(MASS)
?pwr.t.test

# calculate power for n=10, sig.level=0.05, and d=1
#delta = difference between mu and xbar
pwr.t.test(power = 0.8, sig.level = 0.05, d=1, type = "two.sample")

# what if we wanted to know the minimum detectable effect size for 0.80 power
pwr.t.test(n=10,sig.level=0.05,power=0.8)

# now, let's calculate power for a collection of d values
dvals=seq(0,3.5,length=101)
pwr=rep(0,length(dvals))
for(i in 1:(length(dvals))) pwr[i]=pwr.t.test(n=10,d=dvals[i],sig.level=0.05)$power


# make figure
plot(dvals,pwr,type="l",xlab=expression(delta),ylab="power",lwd=2,ylim=c(0,1))
abline(h=0.05,lty=2)


#if the test of truly DE target was conducted as part of the set of 200 tests across all features
## then we would control the FWER using a Bonferroni correction
# we could conduct each of the 200 tests at an alpha of 0.05/200
# and FWER would be controlled no greater than 0.05

######
## QUESTION: What would be the power to reject the null hypothesis if the test was conducted at the Bonferroni correct alpha 0.05/200
######

pwrBon=rep(0,length(dvals))
for(i in 1:(length(dvals))) pwrBon[i]=pwr.t.test(n=10,d=dvals[i],sig.level=0.05/200)$power

# black curve corresponds to alpha=0.05 and blue curve corresponds to alpha=0.05/200


## So what would the power curve look like if we adjusted for multiplicity using maxT instead of bonferroni?
# REMEMBER: Bonferroni correction is reasonable if the set of tests are independent but it can be very conservative if the tests are correlated
#maxT approach -- adjusts for multiplicity while accounting for the underlying correlation structure
#therefore, it is reasonable to expect that the power curve for a maxT approach will lie someone between the two power curves

#==========================================================#
#
# Simple Multivariate Normal Simulation
#
#==========================================================#
# A simulation function to generate a matrix of data:
# the simulated data will be a matrix that has:
#
# 2*n columns with the first n columns corresponding to n samples from a 'control' group
# and the second n columns correspnding to n samples from a 'treatment' group
#
# p.alt+p.null rows with the first p.alt rows corresponding to the the p.alt features (e.g. genes)
# that are diferrentially expressed and the last p.null rows corresponding to
# p.null features that aren't differentiall expressed
#
# the matrix values will all be generated from a multivariate normal distribution with covariance
# matrix:
# | 1 rho.alt rho.alt ... 0 0 0 |
# | rho.alt 1 rho.alt ... 0 0 0 |
# | : : : : : |
# | 0 0 0 ... 1 rho.nul rho.nul |
# | 0 0 0 ... rho.nul 1 rho.nul |
# | 0 0 0 ... rho.nul rho.nul 1 |
#
# The p.null rows will all have mean zero while the p.alt rows will have mean zero for
# the first n columns and mean delta for the last n columns
# (i.e., the columns corresponding to 'treatment').
#
# Additionally, the function allows for sample specific centering error
# from a normal distribution with mean 0 and sd = sdC


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

## for example we can run the following simulation


if(F){ # I ran this as if(T) the first time, then changed it
        # so I don't have to rerun it every time
        nreps <- 8
        f <- factor(c(rep("cntrl", 10), rep("trt", 10)))
        PhiVec <- rep(NA,nreps)
        for(k in 1:nreps){
                X=SimDatHW2(p.alt=1,p.null=199,n=10,rho.alt=0.0,rho.null=0.75,delta=2)
                # calculate minP adjusted p-values
                resT <- mt.maxT(X, f, B = 2500)
                PhiVec[k]=as.integer(resT$adjp[resT$index==1]<=0.05)
        }
        save(file="PhiVec.RData",PhiVec)
} # end if(F)
load(file="PhiVec.RData") # provides PhiVec


plot(dvals,pwr,type="l",xlab=expression(delta),ylab="power",ylim=c(0,1))
abline(h=0.05,lty=2)
lines(dvals,pwrBon,col=4,lwd=2)
abline(h=0.05/200,lty=3,col=4)
# add our simulation point
points(2,mean(PhiVec),pch="+",col="darkgreen", cex=2)


## Run a simulation and add some more points to Figure 1. Specifically, perform simulations for rho is element of the set 0.25, 0.75
# and delta [0.75, 1.50, 2.25, 3.0]

PhiVecs <- function(rho.null, delta){
        f <- factor(c(rep("cntrl", 10), rep("trt", 10)))
        nreps <- 1000
        PhiVec <- rep(NA,nreps)
        for(k in 1:nreps){
                X <- SimDatHW2(p.alt=1,p.null=199,n=10,rho.alt=0,rho.null=rho.null,delta=delta)
                # calculate minP adjusted p-values
                resT <- mt.maxT(X, f, B = 2500)
                PhiVec[k] <- as.integer(resT$adjp[resT$index==1]<=0.05)
        }
        PhiVec
}

deltas <- c(0.75, 1.50, 2.25, 3.0)
rho.null <- c(0.25, 0.75)

#phivecs.delta.rho25 <- mapply(PhiVecs, rho.null[1], deltas[1:4])
#phivecs.delta.rho75 <- mapply(PhiVecs, rho.null[2], deltas[1:4])

#save(file= "phivecs_delta.rho25.RData",phivecs.delta.rho25)
#save(file= "phivecs.delta.rho75.RData",phivecs.delta.rho75)
load("phivecs_deltarho25.RData")
load("phivecs.delta.rho75.RData")

for(i in 1:length(deltas)){
        #colpal <- c('#e41a1c','#377eb8','#4daf4a','#984ea3')
        colpal <- c('#1b9e77','#d95f02','#7570b3','#e7298a')
        points(deltas[i], mean(phivecs.delta.rho75[,i]), pch="*", col=colpal[i], cex=3)
        #legend("topleft", pch="+", paste("delta =", deltas, sep = " "), col=colpal)
}

#Interpret your results
#what can you say about the gains in power when using max instead of Bonferroni
#Do the results for different values of rho make sense?

# Give a short lecture on power and transition into talking about effect sizes (quickR power package)
#Cohen's Effect Size - based on sum of less squares of some distances

#take power curve and ... look at power detecting 1 ... feature selection
# 1 of 1000, 2 of 1000, 3 of 1000 ....
#Power to detect 1 of 1, 1 of 2, 1 of 3 .... 1-1-p^1 ... 1-1-p^2 .... 1-1-p^3 ... etc


#power in terms of when people compare different testing methods
#simulate another power curve for wilcoxon-rank sum test (will have lower power than t-test)
#more robust to outliers ... 
#do a miracle test (a t-test with a level of 0.1) ... it will be more powerful but is it meaningful?





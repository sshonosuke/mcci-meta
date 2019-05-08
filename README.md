# mcci-meta: A Unified Method for Improved Inference in Random-effects Meta-analysis via Monte Carlo Conditioning   
This package implements the unified method for improved inference via Monte Carlo Conditioning (MC) in random-effects meta-analysis, as proposed by Sugasawa and Noma (2017).

Functions are implemented in *MC-function.R*.
Install R functions and example datasets.
```{r}
load("Example-Data.RData") 
source("MC-function.R")
```

# Univariate meta-analysis
```{r}
library(metafor)
mag=magnesium
```

Apply the proposed method.
```{r}
set.seed(1)
EX=MA.ECI(mag$y,mag$v,mc=1000)
EX.CI=exp(EX)
```

Apply other existing methods.
```{r}
# Darsimonian-Laird method
DL=rma.uni(mag$y,mag$v,method="DL")
DL.CI=exp(c(DL$ci.lb,DL$ci.ub))

# Likelihood ratio method
LR=MA.LRCI(mag$y,mag$v)
LR.CI=exp(LR)

# Maximum likelihood method
ML=rma.uni(mag$y,mag$v,method="ML")

# Restricted maximum likelihood method
REML=rma.uni(mag$y,mag$v,method="REML")
REML.CI=exp(c(REML$ci.lb,REML$ci.ub))

# Knapp-Hurtung method
KNHA=rma.uni(mag$y,mag$v,test="knha")
KNHA.CI=exp(c(KNHA$ci.lb,KNHA$ci.ub))

# Peto method
Y=mag$y; W=1/mag$v
hmu=sum(W*Y)/sum(W); sd=sqrt(1/sum(W))
PTE.CI=exp(c(hmu+qnorm(0.025)*sd,hmu+qnorm(0.975)*sd))
```

Summarize the results.
```{r}
EST=c(exp(ML$b),exp(DL$b),exp(ML$b),exp(REML$b),exp(hmu))
CI=rbind(EX.CI,DL.CI,LR.CI,REML.CI,PTE.CI)
result=cbind(EST,CI)
dimnames(result)[[1]]=c("EX","DL","LR","REML","PTE")
dimnames(result)[[2]]=c("Estimate","CI-lower","CI-upper")
round(result,3)
```

# Bivariate meta-analysis
```{r}
set.seed(1)
library(mada)
data(AuditC)
dd=madad(AuditC)

logit=function(x){ log(x/(1-x)) }
invlogit=function(x){ exp(x)/(1+exp(x)) }

pp=cbind(dd$sens$sens,1-dd$spec$spec)
y=logit(pp)
ni=cbind(AuditC$TP+AuditC$FN,AuditC$FP+AuditC$TN)
S=1/(pp*(1-pp)*ni)

# Proposed method
ECR=ECR.BMA(y,S,mc=100,R=5,grid=200,alpha=0.05) 

MA=function(x,n){ 
  a=(n-1)/2; m=length(x)
  xx=c(x[(m-a+1):m],x,x[1:a])
  as.vector(filter(xx,rep(1,n))/n)[(a+1):(m+a)]
}

m=7
sECR=cbind(MA(ECR[,1],m),MA(ECR[,2],m))

# Reitsuma's method
fits=reitsma(AuditC,method="reml")

# Plot
plot(fits,sroclwd=2,main ="",
     ylim=c(0.5,1),xlim=c(0,0.52))
points(invlogit(sECR)[,c(2,1)],type="l",lty=2)
points(invlogit(y)[,c(2,1)],pch=2)
legend("bottomright",c("Data point","Summary estimate","SROC","Approximate CR","Proposed CR"),
       pch=c(2,1,NA,NA,NA),lwd=c(NA,NA,2,1,1),lty=c(NA,NA,1,1,2))
```


# Network meta-analysis
```{r}
set.seed(1)
Y=SCZ[[1]]; X=SCZ[[2]]; Z=SCZ[[3]]; S=SCZ[[4]]
p=dim(X)[2]
MC=1000
alpha=0.05

IC=function(cc){
  C=matrix(0,p,p); C[1,]=cc
  ind=(1:p)[cc!=1]
  for(j in 2:p){
    C[j,ind[j-1]]=1
  }
  solve(C)
}

cc.set=c()
for(j in 1:p){
  cc=rep(0,p); cc[j]=1
  cc.set=rbind(cc.set,cc)
}

Lab=dimnames(X)[[2]]

ECI=matrix(NA,p,2); RMCI=matrix(NA,p,2)
MLCI=matrix(NA,p,4)
dimnames(ECI)[[2]]=dimnames(RMCI)[[2]]=c("Lower","Upper")
dimnames(MLCI)[[2]]=c("Lower","Upper","Lower","Upper")
dimnames(ECI)[[1]]=dimnames(RMCI)[[1]]=dimnames(MLCI)[[1]]=Lab
Mu.ML=c(); Mu.RML=c()


# REML method
res=NMA.RML(Y,X,Z,S)
rml.beta=res[[1]]; ACV=res[[2]]
ml.beta=NMA(Y,X,Z,S)[1:p]

for(i in 1:p){
  invC=IC(cc.set[i,])
  ECI[i,]=NMA.CI(Y,X%*%invC,Z,S,mc=MC,alpha=alpha)[[1]]  # Proposed method
  Mu.ML[i]=t(cc.set[i,])%*%ml.beta
  Mu.RML[i]=t(cc.set[i,])%*%rml.beta
  s=sqrt(t(cc.set[i,])%*%ACV%*%cc.set[i,])
  RMCI[i,]=c(Mu.RML[i]+s*qnorm(alpha/2),Mu.RML[i]+s*qnorm(1-alpha/2))
  MLCI[i,]=NMA.LRCI(Y,X%*%invC,Z,S,alpha=alpha)
  print(i)
}

# Results
CI=cbind(ECI,MLCI[,1:2],RMCI)
dimnames(CI)[[2]]=c("ECI-low","ECI-up","LR-low","LR-up","RML-low","RML-up")
exp(CI)


```

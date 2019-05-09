# mcci-meta: A Unified Method for Improved Inference in Random-effects Meta-analysis via Monte Carlo Conditioning   
This package implements the unified method for improved inference via Monte Carlo Conditioning (MC) in random-effects meta-analysis, as proposed by the following papers.

Sugasawa, S. and Noma, H. (2017). A Unified Method for Improved Inference in Random-effects Meta-analysis.  https://arxiv.org/abs/1711.06393

Functions are implemented in *MC-function.R*.
This tutorial replicates the results in Sections 3.4, 4.4 and 5.4 in the paper.
First, install R functions and example datasets.
```{r}
load("Example-Data.RData") 
source("MC-function.R")
```

# Univariate meta-analysis
The dataset used here consists 7 trials of the treatment of suspected acute myocardial infarction with intravenous magnesium. 

Apply the proposed method.
- `y`: vector of estimates
- `Di`: vector of estimated variances
- `alpha`: significance level
- `mc`: number of Monte Calro samples
```{r}
set.seed(1)
EX=MA.ECI(y=mag$y,Di=mag$v,alpha=0.05,mc=1000)
EX.CI=exp(EX)
```

Apply other existing methods.
```{r}
library(metafor)

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
We analyze the dataset including 14 studies regarding a short screening test for alcohol problems.
The dataset is available in R package *mada*.

Compute estimates of log-odds ratios and their associated standard errors.
```{r}
library(mada)
data(AuditC)
dd=madad(AuditC)

logit=function(x){ log(x/(1-x)) }
invlogit=function(x){ exp(x)/(1+exp(x)) }

pp=cbind(dd$sens$sens,1-dd$spec$spec)
y=logit(pp)
ni=cbind(AuditC$TP+AuditC$FN,AuditC$FP+AuditC$TN)
S=1/(pp*(1-pp)*ni)
```

Apply the proposed method using `ECR.BMA`. (It would take much time when large values of *mc*, *R* and *grid* are used.)
- `y`: (k,2) matrix of estimates (k: number of studies)
- `S`: (k,2) matrix of estimated variances
- `mc`: number of Monte Carlo samples
- `R`: number of iterations in the bisectional method
- `grid`: number of points to approximate confidence boundary
- `alpha`: significance level
```{r}
set.seed(2)
ECR=ECR.BMA(y=y,S=S,mc=100,R=5,grid=200,alpha=0.05) 
```

Smoothing the resutling confidence limits via moving average.
```{r}
MA=function(x,n){ 
  a=(n-1)/2; m=length(x)
  xx=c(x[(m-a+1):m],x,x[1:a])
  as.vector(filter(xx,rep(1,n))/n)[(a+1):(m+a)]
}

m=7  # number of points to compute local average
sECR=cbind(MA(ECR[,1],m),MA(ECR[,2],m))
```

Apply the Reitsuma's method.
```{r}
fits=reitsma(AuditC,method="reml")
```

Make a figure summarizing the results.
```{r}
plot(fits,sroclwd=2,main ="",
     ylim=c(0.5,1),xlim=c(0,0.52))
points(invlogit(sECR)[,c(2,1)],type="l",lty=2)
points(invlogit(y)[,c(2,1)],pch=2)
legend("bottomright",c("Data point","Summary estimate","SROC","Approximate CR","Proposed CR"),
       pch=c(2,1,NA,NA,NA),lwd=c(NA,NA,2,1,1),lty=c(NA,NA,1,1,2))
```


# Network meta-analysis
We consider network meta-analysis of antipsychotic medication for prevention of relapse of schizophrenia. This analysis includes 15 trials comparing eight treatments with placebo.
The dataset `SCZ` and list of matrices `invC` are incldued in `Example-Data.RData`.

```{r}
Y=SCZ[[1]]; X=SCZ[[2]]; Z=SCZ[[3]]; S=SCZ[[4]]
p=dim(X)[2]

Lab=dimnames(X)[[2]]
ECI=matrix(NA,p,2); RMCI=matrix(NA,p,2)
dimnames(ECI)[[2]]=dimnames(RMCI)[[2]]=c("Lower","Upper")
dimnames(ECI)[[1]]=dimnames(RMCI)[[1]]=Lab
```

Apply the REML method.
```{r}
res=NMA.RML(Y,X,Z,S)
Mu.RML=res[[1]]; ss=sqrt(diag(res[[2]]))
RMCI[,1]=Mu.RML+qnorm(0.05/2)*ss
RMCI[,2]=Mu.RML+qnorm(1-0.05/2)*ss
```

Apply the proposed method using the function `NMA.CI`.
- `Y`: vector of estimates
- `X`: (N,p) design matrix for fixed effects
- `Z`: (N,N*p) design matrix for random effects
- `S`: vector of estimated varinces
- `mc`: number of Monte Carlo samples
- `alpha`: significance level
```{r}
set.seed(3)
for(i in 1:p){
  ECI[i,]=NMA.CI(Y=Y,X=X%*%invC[[i]],Z=Z,S=S,mc=100,alpha=0.05)[[1]] 
  print(i)
}
```

Summarizing the results.
```{r}
CI=cbind(ECI,RMCI)
dimnames(CI)[[2]]=c("ECI-low","ECI-up","RML-low","RML-up")
exp(CI)
```

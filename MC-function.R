library(MASS)
library(mvmeta)



#-------------------------------------------------#
#          Univariate Meta-analysis               #
#-------------------------------------------------#

###  MC p-value 
# y: vector of estimates
# Di: vector of estimated variances
# mu0: value of treatment effect to be tested
# mc: number of Monte Carlo samples
MA.ET=function(y,Di,mu0=0,mc=1000){
  n=length(y)
  like=function(y,para){
    mu=para[1]; A=para[2]
    sum(log(A+Di))+sum((y-mu)^2/(A+Di))
  }
  
  stat=function(obs){
    opt1=function(para){ like(obs,para) }
    out1=optim(par=c(1,1),fn=opt1,method="L-BFGS-B",lower=c(-100,0),upper=c(100,100))
    opt2=function(A){ like(obs,c(mu0,A)) }
    out2=optim(par=1,fn=opt2,method="L-BFGS-B",lower=0,upper=100)
    return(-out1$value+out2$value)
  }
  
  cml=function(A){ like(y,c(mu0,A)) }
  hA=optim(par=1,fn=cml,method="L-BFGS-B",lower=0,upper=100)$par
  origin=stat(y)
  
  num=c(); denom=c()
  for(k in 1:mc){
    u=rnorm(n)
    a1=sum(u^2/(hA+Di)^2); a2=sum(1/(hA+Di))-sum(Di*u^2/(hA+Di)^2)
    tA=max(a2/a1,0)
    b1=sum((2*(tA+Di)*u^2-(hA+Di))/(hA+Di)^3); b2=sum(u^2/(hA+Di)^2)
    W=abs(b1)/b2
    yt=mu0+sqrt(tA+Di)*u; st=stat(yt)
    num[k]=ifelse(st>origin,1,0)*W; denom[k]=W
  }
  
  p.val=mean(num)/mean(denom)
  return(p.val)
}



### Confidece interval based on MC-pvalue
# y: vector of estimates
# Di: vector of estimated variances
# alpha: significance level
# mc: number of Monte Carlo samples
MA.ECI=function(y,Di,alpha=0.05,mc=1000){
  n=length(y)
  like=function(para){
    mu=para[1]; A=para[2]
    sum(log(A+Di))+sum((y-mu)^2/(A+Di))
  }
  out=optim(par=c(1,1),fn=like,method="L-BFGS-B",lower=c(-100,0),upper=c(100,100))$par
  hmu=out[1]; hA=out[2]
  et=function(mu){ MA.ET(y,Di,mu0=mu,mc=mc) }
  s=sqrt(1/sum(1/(hA+Di)))
  
  ll=hmu-4*s
  while(et(ll)>alpha){ ll=ll-s }
  L.int=c(ll,hmu)
  ul=hmu+4*s
  while(et(ul)>alpha){ ul=ul+s }
  U.int=c(hmu,ul)
  
  L=mean(L.int); d=1
  while(d>0.05){
    pv=et(L)
    if(pv<alpha){ new.L=(L+L.int[2])/2; L.int[1]=L }
    if(pv>alpha){ new.L=(L+L.int[1])/2; L.int[2]=L }
    if(pv==alpha){ new.L=L }
    d=abs(new.L-L)/abs(L+0.00001)*100
    L=new.L
  }
  
  U=mean(U.int); d=1
  while(d>0.05){
    pv=et(U)
    if(pv<alpha){ new.U=(U+U.int[1])/2; U.int[2]=U }
    if(pv>alpha){ new.U=(U+U.int[2])/2; U.int[1]=U }
    if(pv==alpha){ new.U=U }
    d=abs(new.U-U)/abs(U+0.00001)*100
    U=new.U
  }
  return(c(L,U))
}



### Confidece interval based on LRT
# y: vector of estimates
# Di: vector of estimated variances
# alpha: significance level
MA.LRCI=function(y,Di,alpha=0.05){
  n=length(y)
  like=function(para){
    mu=para[1]; A=para[2]
    sum(log(A+Di))+sum((y-mu)^2/(A+Di))
  }
  
  LRT=function(mu0=0){
    out1=optim(par=c(1,1),fn=like,method="L-BFGS-B",lower=c(-100,0),upper=c(100,100))
    opt2=function(A){ like(c(mu0,A)) }
    out2=optim(par=1,fn=opt2,method="L-BFGS-B",lower=0,upper=100)
    stat=-out1$value+out2$value
    p.val=1-pchisq(stat,1)
    return(p.val)
  }
  
  out=optim(par=c(1,1),fn=like,method="L-BFGS-B",lower=c(-100,0),upper=c(100,100))$par
  hmu=out[1]; hA=out[2]
  PV=function(mu){ LRT(mu0=mu) }
  
  s=sqrt(1/sum(1/(hA+Di)))
  L.int=c(hmu-4*s,hmu)
  U.int=c(hmu,hmu+4*s)
  
  L=mean(L.int); d=1
  while(d>0.01){
    pv=PV(L)
    if(pv<alpha){ new.L=(L+L.int[2])/2; L.int[1]=L }
    if(pv>alpha){ new.L=(L+L.int[1])/2; L.int[2]=L }
    if(pv==alpha){ new.L=L }
    d=abs(new.L-L)/abs(L+0.00001)*100
    L=new.L
  }
  
  U=mean(U.int); d=1
  while(d>0.01){
    pv=PV(U)
    if(pv<alpha){ new.U=(U+U.int[1])/2; U.int[2]=U }
    if(pv>alpha){ new.U=(U+U.int[2])/2; U.int[1]=U }
    if(pv==alpha){ new.U=U }
    d=abs(new.U-U)/abs(U+0.00001)*100
    U=new.U
  }
  return(c(L,U))
}




#-------------------------------------------------#
#           Bivariate Meta-analysis               #
#-------------------------------------------------#

###  Fiting bivariate model
# y: (k,2) matrix of estimates (k: number of studies)
# S: (k,2) matrix of estimated variances
BMA=function(y,S){
  res=mvmeta(y,S,method="ml")
  mu=as.vector(res$coefficients)
  vv=res$Psi
  tauA=sqrt(vv[1,1]); tauB=sqrt(vv[2,2])
  rho=vv[1,2]/(tauA*tauB)
  est=c(mu,tauA,tauB,rho)
  names(est)=c("mu-A","mu-B","tau-A","tau-B","rho")
  return(est)
}


###   MC p-value of LRT
# y: (k,2) matrix of estimates (k: number of studies)
# S: (k,2) matrix of estimated variances
# mu0: vector of values to be tested
# mc: number of Monte Carlo samples
EP.BMA=function(y,S,mu0,mc=100){
  init=BMA(y,S)[3:5]
  hSig=matrix(c(init[1]^2,prod(init),prod(init),init[2]^2),2,2)
  
  rBMA=function(y,S,mu0){
    n=dim(y)[1]
    Sig=hSig
    diff=1; R=1
    while(diff>0.01){
      val2=array(NA,c(2,2,n))
      for(i in 1:n){
        invCi=diag(1/S[i,]); invSig=ginv(Sig)
        mat=solve(invCi+invSig)
        tht=mat%*%as.vector(t(y[i,])%*%invCi+t(mu0)%*%invSig)
        val2[,,i]=mat+(tht-mu0)%*%t(tht-mu0)
      }
      new.Sig=apply(val2,c(1,2),mean)
      diff=100*sum((new.Sig-Sig)^2)/sum(Sig^2)
      R=R+1; if(R>20){ diff=0 }
      Sig=new.Sig
    }
    tauA=sqrt(Sig[1,1]); tauB=sqrt(Sig[2,2]); rho=Sig[1,2]/(tauA*tauB)
    rest=c(tauA,tauB,rho)
    names(rest)=c("tau-A","tau-B","rho")
    return(rest)
  }
  
  LR.stat=function(y,S,mu0){
    n=dim(y)[1]
    nloglike=function(yy,S,para){
      mu=para[1:2]; tauA=para[3]; tauB=para[4]; rho=para[5]
      Sig=matrix(c(tauA^2,tauA*tauB*rho,tauA*tauB*rho,tauB^2),2,2)
      val=c()
      for(i in 1:n){
        mat=Sig+diag(S[i,])
        val[i]=log(det(mat))+(yy[i,]-mu)%*%solve(mat)%*%(yy[i,]-mu)
      }
      return(sum(val))
    }
    ML=BMA(y,S); rML=rBMA(y,S,mu0=mu0)
    q1=nloglike(y,S,ML); q2=nloglike(y,S,c(mu0,rML))
    return(q2-q1)
  }
  
  n=dim(y)[1]
  ml=BMA(y,S); rml=rBMA(y,S,mu0)
  origin=LR.stat(y,S,mu0)
  
  tr=function(mat){ sum(diag(mat)) }
  
  J1=matrix(c(1,0,0,0),2,2)
  J2=matrix(c(0,1,1,0),2,2)
  J3=matrix(c(0,0,0,1),2,2)
  e1=c(1,0,0); e2=c(0,1,0); e3=c(0,0,1)
  z=0.001; ad=1.1
  
  EE=function(U,vv,hvv){
    tauA=vv[1]; tauB=vv[2]; rho=vv[3]
    Sig=matrix(c(tauA^2,tauA*tauB*rho,tauA*tauB*rho,tauB^2),2,2)
    htauA=hvv[1]; htauB=hvv[2]; hrho=hvv[3]
    hSig=matrix(c(htauA^2,htauA*htauB*hrho,htauA*htauB*hrho,htauB^2),2,2)
    val=matrix(NA,n,3)
    for(i in 1:n){
      A=solve(hSig+diag(S[i,]))
      B=t(chol(Sig+diag(S[i,])))
      a=A%*%B%*%U[i,]
      val[i,1]=A[1,1]-a[1]^2
      val[i,2]=A[2,2]-a[2]^2
      val[i,3]=A[1,2]-a[1]*a[2]
    }
    eq=apply(val,2,mean)
    return(sum(eq^2))
  }
  
  boot.vv=function(U,hvv){
    OB=function(x){ EE(U,x,hvv) }
    optim(par=hvv,OB,lower=c(ad*z,ad*z,-1+ad*z),upper=c(100,100,1-ad*z),method="L-BFGS-B")$par
  }
  
  est.vv=function(U,vv){
    tauA=vv[1]; tauB=vv[2]; rho=vv[3]
    Sig=matrix(c(tauA^2,tauA*tauB*rho,tauA*tauB*rho,tauB^2),2,2)
    yy=matrix(NA,n,2)
    for(i in 1:n){
      yy[i,]=mu0+t(chol(Sig+diag(S[i,])))%*%U[i,]
    }
    return(rBMA(yy,S,mu0))
  }
  
  weight=c(); st=c()
  for(k in 1:mc){
    U=matrix(rnorm(2*n),n,2)
    ast.vv=boot.vv(U,rml)
    dd=matrix(NA,3,3)
    dd[,1]=(est.vv(U,ast.vv+e1*z)-est.vv(U,ast.vv-e1*z))/(2*z)
    dd[,2]=(est.vv(U,ast.vv+e2*z)-est.vv(U,ast.vv-e2*z))/(2*z)
    dd[,3]=(est.vv(U,ast.vv+e3*z)-est.vv(U,ast.vv-e3*z))/(2*z)
    weight[k]=1/abs(det(dd))
    1/abs(det(dd))
    
    atauA=ast.vv[1]; atauB=ast.vv[2]; arho=ast.vv[3]
    ast.Sig=matrix(c(atauA^2,atauA*atauB*arho,atauA*atauB*arho,atauB^2),2,2)
    boot.y=matrix(NA,n,2)
    for(i in 1:n){
      boot.y[i,]=mu0+t(chol(ast.Sig+diag(S[i,])))%*%U[i,]
    }
    st[k]=LR.stat(boot.y,S,mu0)
  }
  weight[weight>5]=5
  pp=mean(ifelse(st>origin,1,0)*weight)/mean(weight)
  return(pp)
}



###  Exact confidece region
# y: (k,2) matrix of estimates (k: number of studies)
# S: (k,2) matrix of estimated variances
# mc: number of Monte Carlo samples
# R: number of itertaions to compute boundary for each grid 
# grid: number of grid to approximate region boundary
# alpha: significance level
ECR.BMA=function(y,S,mc=200,R=5,grid=50,alpha=0.05){
  et=function(mu){ EP.BMA(y,S,mu0=mu,mc=mc) }
  t=seq(0,2*pi,length=grid+1)[-(grid+1)]
  
  res=mvmeta(y,S)
  hmu=as.vector(res$coefficients)
  vv=res$Psi
  tauA=sqrt(vv[1,1]); tauB=sqrt(vv[2,2])
  hrho=vv[1,2]/(tauA*tauB)
  AC=res$vcov
  dist=sqrt(diag(AC))*qchisq(1-alpha,2)*1.2
  
  rlen=c()
  for(k in 1:grid){
    dir=c(cos(t[k]),sin(t[k]))
    dd=dist*c(cos(t[k]),cos(t[k]+acos(hrho)))
    r1=0; r2=sqrt(sum(dd^2))
    root=c(); root[1]=(r1+r2)/2
    for(r in 1:R){
      ppc=et(hmu+root[r]*dir)
      if(ppc<alpha){ root[r+1]=(root[r]+r1)/2; r2=root[r] }
      if(ppc>alpha){ root[r+1]=(root[r]+r2)/2; r1=root[r] }
      if(ppc==alpha){ root[r+1]=root[r] }
    }
    rlen[k]=root[R]
    print(k)
  }
  
  muA=hmu[1]+rlen*cos(t); muB=hmu[2]+rlen*sin(t)
  return(cbind(muA,muB))
}




#-------------------------------------------------#
#           Network Meta-analysis                 #
#-------------------------------------------------#

### Fitting network model
# Y: vector of estimates
# X: (N,p) design matrix for fixed effects
# Z: (N,N*p) design matrix for random effects
# S: vector of estimated varinces
NMA=function(Y,X,Z,S){
  p=dim(X)[2]; n=dim(Z)[2]/p
  P=matrix(0.5,p,p); diag(P)=1
  Q=Z%*%(diag(n)%x%P)%*%t(Z)
  
  ML=function(beta,tau){
    mat=tau^2*Q+S
    resid=as.vector(Y-X%*%beta)
    val=log(det(mat))+t(resid)%*%solve(mat)%*%resid
    return(as.vector(val))
  }
  
  GLS=function(tau){
    mat=tau^2*Q+S; invmat=solve(mat)
    hbeta=as.vector(solve(t(X)%*%invmat%*%X)%*%t(X)%*%invmat%*%Y)
    return(hbeta)
  }
  
  PL=function(tau){ ML(GLS(tau),tau) }
  htau=optim(par=1,fn=PL,method="L-BFGS-B",lower=0,upper=100)$par
  hbeta=GLS(htau)
  Est=c(hbeta,htau)
  names(Est)=c(paste0("beta",1:p),"tau")
  return(Est)
}


###   Fitting network model with fixed beta1
# Y: vector of estimates
# X: (N,p) design matrix for fixed effects
# Z: (N,N*p) design matrix for random effects
# S: vector of estimated varinces
# beta1: value of fixed effect to be tested
rNMA=function(Y,X,Z,S,beta1=0){
  p=dim(X)[2]; n=dim(Z)[2]/p
  P=matrix(0.5,p,p); diag(P)=1
  Q=Z%*%(diag(n)%x%P)%*%t(Z)
  W1=X[,1]; W2=X[,-1]
  
  ML=function(om,tau){
    mat=tau^2*Q+S
    resid=as.vector(Y-W1*beta1-W2%*%om)
    val=log(det(mat))+t(resid)%*%solve(mat)%*%resid
    return(as.vector(val))
  }
  
  GLS=function(tau){
    mat=tau^2*Q+S; invmat=solve(mat)
    hbeta=as.vector(solve(t(W2)%*%invmat%*%W2)%*%t(W2)%*%invmat%*%(Y-W1*beta1))
    return(hbeta)
  }
  
  PL=function(tau){ ML(GLS(tau),tau) }
  htau=optim(par=1,fn=PL,method="L-BFGS-B",lower=10^(-10),upper=100)$par
  hom=GLS(htau)
  Est=c(hom,htau)
  names(Est)=c(paste0("beta",2:p),"tau")
  return(Est)
}




###   Likelihood ratio test statistics
# Y: vector of estimates
# X: (N,p) design matrix for fixed effects
# Z: (N,N*p) design matrix for random effects
# S: vector of estimated varinces
# beta1: value of fixed effect to be tested
LRT=function(Y,X,Z,S,beta1=0){
  p=dim(X)[2]; n=dim(Z)[2]/p
  P=matrix(0.5,p,p); diag(P)=1
  Q=Z%*%(diag(n)%x%P)%*%t(Z)
  
  LL=function(Est){
    beta=Est[1:p]; tau=Est[p+1]
    mat=tau^2*Q+S
    resid=as.vector(Y-X%*%beta)
    val=-0.5*log(det(mat))-0.5*t(resid)%*%solve(mat)%*%resid
    return(as.vector(val))
  }
  est1=rNMA(Y,X,Z,S,beta1=beta1); est2=NMA(Y,X,Z,S)
  stat=-2*(LL(c(beta1,est1))-LL(est2))
  return(stat)
}


###  MC p-value for LRT
# Y: vector of estimates
# X: (N,p) design matrix for fixed effects
# Z: (N,N*p) design matrix for random effects
# S: vector of estimated varinces
# beta1: value of fixed effect to be tested
# mc: number of Monte Carlo samples
NMA.EP=function(Y,X,Z,S,beta1=0,mc=200){
  p=dim(X)[2]; n=dim(Z)[2]/p
  P=matrix(0.5,p,p); diag(P)=1
  Q=Z%*%(diag(n)%x%P)%*%t(Z)
  W1=X[,1]; W2=X[,-1]
  N=length(Y)
  origin=LRT(Y,X,Z,S,beta1=beta1)
  
  tr=function(mat){ sum(diag(mat)) }
  
  rest=rNMA(Y,X,Z,S,beta1=beta1)
  omh=rest[1:(p-1)]; tauh=rest[p]
  Vh=tauh^2*Q+S; IV=solve(Vh)
  Bh=diag(N)-W2%*%solve(t(W2)%*%IV%*%W2)%*%t(W2)%*%IV
  z=0.001; th=10
  lower=max(z+0.001,tauh^2/2)
  
  St=c(); weight=c()
  for(k in 1:mc){
    u=rnorm(N)
    Eq=function(tau2){
      A=t(chol(tau2*Q+S))
      tr(IV%*%Q)-t(u)%*%t(A)%*%Bh%*%IV%*%Q%*%IV%*%Bh%*%A%*%u
    }
    Op=function(x){ Eq(x)^2 }
    tau2.ast=optim(par=tauh^2,fn=Op,method="L-BFGS-B",lower=lower,upper=100)$par
    
    A.ast=t(chol(tau2.ast*Q+S))
    om.ast=as.vector(omh-solve(t(W2)%*%IV%*%W2)%*%t(W2)%*%IV%*%A.ast%*%u)
    Y.ast=as.vector(W1*beta1+W2%*%om.ast+A.ast%*%u)
    St[k]=LRT(Y.ast,X,Z,S,beta1=beta1)
    
    Mat1=matrix(NA,p,p); Mat2=matrix(NA,p,p)
    Mat1[1:(p-1),1:(p-1)]=t(W2)%*%IV%*%W2
    Mat1[p,1:(p-1)]=-2*t(W2)%*%IV%*%Q%*%IV%*%(W2%*%(om.ast-omh)+A.ast%*%u)
    Mat2[1:(p-1),1:(p-1)]=-Mat1[1:(p-1),1:(p-1)]
    Mat2[p,1:(p-1)]=-Mat1[p,1:(p-1)]
    
    A1=t(chol((tau2.ast+z)*Q+S))
    A2=t(chol((tau2.ast-z)*Q+S))
    r1=W2%*%(om.ast-omh)+A1%*%u
    r2=W2%*%(om.ast-omh)+A2%*%u
    Ea1=t(W2)%*%IV%*%r1
    Ea2=t(W2)%*%IV%*%r2
    Mat1[1:(p-1),p]=(Ea1-Ea2)/(2*z)
    Eb1=tr(IV%*%Q)-t(r1)%*%IV%*%Q%*%IV%*%r1
    Eb2=tr(IV%*%Q)-t(r2)%*%IV%*%Q%*%IV%*%r2
    Mat1[p,p]=(Eb1-Eb2)/(2*z)
    
    rr=W2%*%(om.ast-omh)+A.ast%*%u
    IV1=solve((tauh^2+z)*Q+S)
    IV2=solve((tauh^2-z)*Q+S)
    Ea1=t(W2)%*%IV1%*%rr
    Ea2=t(W2)%*%IV2%*%rr
    Mat2[1:(p-1),p]=(Ea1-Ea2)/(2*z)
    Eb1=tr(IV1%*%Q)-t(rr)%*%IV1%*%Q%*%IV1%*%rr
    Eb2=tr(IV2%*%Q)-t(rr)%*%IV2%*%Q%*%IV2%*%rr
    Mat2[p,p]=(Eb1-Eb2)/(2*z)
   
    m1=abs(eigen(Mat1)$values)
    m2=abs(eigen(Mat2)$values)
    exp(sum(log(m2)-log(m1)))
    
    weight[k]=exp(sum(log(m2)-log(m1)))
  }
  weight[weight>th]=th
  pp=mean(ifelse(St>origin,1,0)*weight)/mean(weight)
  return(pp)
}



###  Confidece interval based on MC p-value
# Y: vector of estimates
# X: (N,p) design matrix for fixed effects
# Z: (N,N*p) design matrix for random effects
# S: vector of estimated varinces
# mc: number of Monte Carlo samples
# alpha: significance level
NMA.CI=function(Y,X,Z,S,mc=200,alpha=0.05){
  p=dim(X)[2]; n=dim(Z)[2]/p
  P=matrix(0.5,p,p); diag(P)=1
  Q=Z%*%(diag(n)%x%P)%*%t(Z)
  Est=NMA(Y,X,Z,S)
  hbeta1=Est[1]; htau=Est[p+1]
  Vh=htau^2*Q+S; IV=solve(Vh)
  ACV=solve(t(X)%*%IV%*%X)
  s=sqrt(ACV[1,1])
  NCI=c(hbeta1+s*qnorm(alpha/2),hbeta1+s*qnorm(1-alpha/2))
  
  et=function(bb){ NMA.EP(Y,X,Z,S,beta1=bb,mc=mc) }
  L.int=c(hbeta1-4*s,hbeta1-s)
  U.int=c(hbeta1+s,hbeta1+4*s)
  
  L=mean(L.int); d=1
  while(d>0.01){
    pv=et(L)
    if(pv<alpha){ new.L=(L+L.int[2])/2; L.int[1]=L }
    if(pv>alpha){ new.L=(L+L.int[1])/2; L.int[2]=L }
    if(pv==alpha){ new.L=L }
    d=abs(new.L-L)/abs(L+0.00001)*100
    L=new.L
  }
  
  U=mean(U.int); d=1
  while(d>0.1){
    pv=et(U)
    if(pv<alpha){ new.U=(U+U.int[1])/2; U.int[2]=U }
    if(pv>alpha){ new.U=(U+U.int[2])/2; U.int[1]=U }
    if(pv==alpha){ new.U=U }
    d=abs(new.U-U)/abs(U+0.00001)*100
    U=new.U
  }
  
  ECI=c(L,U)
  Result=list(ECI,NCI); names(Result)=c("ECI","NCI")
  return(Result)
}



###  Fiting network model via REML
# Y: vector of estimates
# X: (N,p) design matrix for fixed effects
# Z: (N,N*p) design matrix for random effects
# S: vector of estimated varinces
NMA.RML=function(Y,X,Z,S){
  p=dim(X)[2]; n=dim(Z)[2]/p
  P=matrix(0.5,p,p); diag(P)=1
  Q=Z%*%(diag(n)%x%P)%*%t(Z)
  
  RML=function(tau){
    V=tau^2*Q+S; IV=solve(V)
    mat=t(X)%*%IV%*%X
    P=IV-IV%*%X%*%solve(mat)%*%t(X)%*%IV
    log(det(mat))+log(det(V))+t(Y)%*%P%*%Y
  }
  
  GLS=function(tau){
    mat=tau^2*Q+S; invmat=solve(mat)
    hbeta=as.vector(solve(t(X)%*%invmat%*%X)%*%t(X)%*%invmat%*%Y)
    return(hbeta)
  }
  
  PL=function(tau){ RML(tau) }
  htau=optim(par=1,fn=PL,method="L-BFGS-B",lower=0,upper=1000)$par
  
  hbeta=GLS(htau)
  Vh=htau^2*Q+S; IV=solve(Vh)
  ACV=solve(t(X)%*%IV%*%X)
  return(list(hbeta,ACV,htau))
}



###  Confidece interval based on likelihood ratio test 
# Y: vector of estimates
# X: (N,p) design matrix for fixed effects
# Z: (N,N*p) design matrix for random effects
# S: vector of estimated varinces
# alpha: significance level
NMA.LRCI=function(Y,X,Z,S,alpha=0.05){
  p=dim(X)[2]; n=dim(Z)[2]/p
  P=matrix(0.5,p,p); diag(P)=1
  Q=Z%*%(diag(n)%x%P)%*%t(Z)
  Est=NMA(Y,X,Z,S)
  hbeta1=Est[1]; htau=Est[p+1]
  Vh=htau^2*Q+S; IV=solve(Vh)
  ACV=solve(t(X)%*%IV%*%X)
  s=sqrt(ACV[1,1])
  NCI=c(hbeta1+s*qnorm(alpha/2),hbeta1+s*qnorm(1-alpha/2))
  
  et=function(bb){ 
    st=LRT(Y,X,Z,S,beta1=bb) 
    1-pchisq(st,1)
  }
  
  L.int=c(hbeta1-4*s,hbeta1)
  U.int=c(hbeta1,hbeta1+4*s)
  L=mean(L.int); d=1
  while(d>0.01){
    pv=et(L)
    if(pv<alpha){ new.L=(L+L.int[2])/2; L.int[1]=L }
    if(pv>alpha){ new.L=(L+L.int[1])/2; L.int[2]=L }
    if(pv==alpha){ new.L=L }
    d=abs(new.L-L)/abs(L+0.00001)*100
    L=new.L
  }
  U=mean(U.int); d=1
  while(d>0.1){
    pv=et(U)
    if(pv<alpha){ new.U=(U+U.int[1])/2; U.int[2]=U }
    if(pv>alpha){ new.U=(U+U.int[2])/2; U.int[1]=U }
    if(pv==alpha){ new.U=U }
    d=abs(new.U-U)/abs(U+0.00001)*100
    U=new.U
  }
  LRCI=c(L,U)
  Result=c(LRCI,NCI)
  return(Result)
}





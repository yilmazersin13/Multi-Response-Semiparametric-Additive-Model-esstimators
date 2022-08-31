library(condSURV)
library(psych)
#DATA GENERATION----------------------------------------------------------------
n <- 100
no.cov <- 2
if(no.cov == 2){
  t<-0
  p<-2
  x<-matrix(c(runif(n),runif(n)),n,2)
  for (j in 1:n){
    t[j]<-(j-0.5)/n
  }
  z<-t+seq(-2,2,length.out = n)
  f1<-1-48*t+218*t^2-315*t^3+145*t^4
  f2<-sin(2*z)+2*exp(-16*z^2)
  beta<-c(-0.5,1)
  error<-rnorm(n,sd=0.5)
  y<-x%*%beta+f1+f2+error
  nc<-matrix(c(t,z),n,2)
  allf<-matrix(c(f1,f2),n,2)
  dat<-new.env()
  dat$x<-x
  dat$y<-y
  dat$nc<-nc
  dat$beta<-beta
  dat$allf<-allf
} 
if(no.cov == 4){
  t2<-0
  t<-0
  x<-matrix(c(runif(n),runif(n)),n,2)
  for (j in 1:n){
    t[j] <-(j-0.5)/n
  }
  t2 <- 5*seq(0,1,length.out = n)
  z<-t+seq(-2,2,length.out = n)
  z2<-seq(-2,2,length.out = n)+(t2/10)
  f1<-1-48*t+218*t^2-315*t^3+145*t^4
  f2<-sin(2*z)+2*exp(-16*z^2)
  f3<-t2*(sin(t2))^2
  f4<-z2+2*exp(-15*z2^2)
  beta<-c(-1,2)
  error<-rnorm(n,sd=0.5)
  y<-x%*%beta+f1+f2+f3+f4+error
  nc<-matrix(c(t,t2,z,z2),n,4)
  allf<-matrix(c(f1,f2,f3,f4),n,4)
}
#-------------------------------------------------------------------------------
z1  <- t       #Nonparametric cov. #1
z2  <- z       #Nonparametric cov. #2
#-------------------------------------------------------------------------------
#BACKFITTING PROCEDURE----------------------------------------------------------
#Initialization-----------------------------------------------------------------
f01 <- fitted(lm(y~z1))
f02 <- fitted(lm(y~z2))
#SMOOTHING MATRIX FOR KS--------------------------------------------------------
selectionKS <- function(x,z,y){
  aic <- 0
  gcv <- 0
  cp  <- 0
  aiccfunc<-function(y,yhat,H){
    p <- tr(H)
    y<-matrix(c(y))
    yhat<-matrix(c(yhat))
    n<-length(y)
    score<-1+log((norm(((H-diag(n))%*%y))^2)/n)+(2*(tr(H)+1)/(n-tr(H)-2))
    return(score)
  }
  #---------------GCV-----------------
  gcvfunc<-function(y,yhat,H){
    y<-matrix(c(y))
    yhat<-matrix(c(yhat))
    n<-length(y)
    score<-(1/n)*(norm(y-yhat)^2)/(((1/n)*tr(diag(n)-H))^2)
    return(score)
  }
  #---------------GCV-----------------
  cpfunc<-function(y,yhat,H){
    y<-matrix(c(y))
    yhat<-matrix(c(yhat))
    n<-length(y)
    MSE<-norm(((diag(n)-H)%*%y)^2)
    ssq<-var(y-yhat)
    DF <- tr(t(diag(n)-H)%*%(diag(n)-H))
    score<-MSE+(2*ssq*DF)
    return(score)
  }
  #-----------------------------------------------------------------------------
  n      <- length(y)
  index  <-seq(min(z)-0.1,max(z)+0.1,length.out=n)
  tp_seq <-seq(0.01,0.5,length.out=50)
  W      <-matrix(0,n,n) 
  for (i in 1:50){
    for (j in 1:n){
      W[,j] <- NWW(z,index[j],bw=tp_seq[i])
    }
    xtil <- (diag(n)-W)%*%x
    ytil <- (diag(n)-W)%*%y
    
    beta <- solve(t(xtil)%*%xtil)%*%t(xtil)%*%ytil
    fhat <- W%*%(y-x%*%beta)
    yhat <- x%*%beta+fhat
    H    <- xtil%*%solve(t(xtil)%*%xtil)%*%t(xtil)
    aic[i] <- aiccfunc(y,yhat,H)
    gcv[i] <- gcvfunc(y,yhat,H)
    cp[i]  <- cpfunc(y,yhat,H)
  }
  for (i2 in 1:50){
    if (aic[i2]==min(aic)){
      lam_aic <- tp_seq[i2]
    }
    if (gcv[i2]==min(gcv)){
      lam_gcv <- tp_seq[i2]
    }
    if (cp[i2]==min(cp)){
      lam_cp <- tp_seq[i2]
    }
  }
  res <- new.env()
  res$lam.aic <- lam_aic
  res$lam.gcv <- lam_gcv
  res$lam.cp  <- lam_cp
  
    par(mfrow=c(1,3))
    plot(aic,type="l")
    plot(gcv,type="l")
    plot(cp,type="l")
  
  return(res)
}
Smat <- function(z,tp){
  n      <- length(z)
  index  <-seq(min(z)-0.1,max(z)+0.1,length.out=n)
  W      <-matrix(0,n,n) 
  for (j in 1:n){
    W[,j] <- NWW(z,index[j],bw=tp)
  }
  return(W)
}
alpha0 <- mean(y)
#-------------------------------------------------------------------------------
fhat1_aic <- f01 #matrix(0,n,1)
fhat2_aic <- f02 #matrix(0,n,1)

fhat1_gcv <- f01
fhat2_gcv <- f02

fhat1_cp <- f01
fhat2_cp <- f02

fhat_aic <- matrix(0,n,2)
fhat_gcv <- matrix(0,n,2)
fhat_cp  <- matrix(0,n,2)

iter <- 100
fh1 <- matrix(0,n,iter)
fh2 <- matrix(0,n,iter)
for (k in 1:2){
  tp <- selectionKS(x,nc[,k],(y))
  tp.aic <- tp$lam.aic
  tp.gcv <- tp$lam.gcv
  tp.cp  <- tp$lam.cp
  S_aic  <- Smat(nc[,k],tp.aic) 
  S_gcv  <- Smat(nc[,k],tp.gcv)
  S_cp   <- Smat(nc[,k],tp.cp)
  tol    <- 0.05
  ctol   <- 99
  i <- 1
while (ctol>=tol){
    
    betai_aic   <- solve(t(x)%*%x)%*%t(x)%*%(y-alpha0-(fhat_aic[,k]+fhat_aic[,-k])) 
    fhat_aic[,k]<- S_aic%*%(y-alpha0-(x%*%betai_aic)-fhat_aic[,-k])
    
    if (k==1){
      fhat1_aic <- fhat_aic[,k]
      fh1[,i]   <- fhat1_aic
    }
    if (k==2){
      fhat2_aic <- fhat_aic[,k]
      fh2[,i]   <- fhat2_aic
    }
    if (i>1){
    tol2      <-(mean(abs(fh2[,(i-1)]-fhat2_aic)))
    tol1      <-(mean(abs(fh1[,(i-1)]-fhat1_aic)))
    ctol <- (tol1+tol2)/2
    }
    i <- i+1
    if (i==iter){
      break
    }
}
  
}
for (i3 in 1:(i-1)){
  plot(fh1[,i3],type="l", ylim=c(min(f1),max(f1)),col = rgb(red = 0, green = 0, blue = 0, alpha = 0.2))
  par(new=TRUE)
  plot(f1,type="l",col=2)
  par(new=TRUE)
}
par(new=FALSE)
for (i3 in 1:(i-1)){
  plot(fh2[,i3],type="l", ylim=c(min(f2)-2,max(f2)),col = rgb(red = 0, green = 0, blue = 0, alpha = 0.2))
  par(new=TRUE)
  plot(f2,type="l",col=3, ylim=c(min(f2)-2,max(f2)))
  par(new=TRUE)
}
par(new=FALSE)

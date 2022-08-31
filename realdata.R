#SIMULATION STUDY FOR THE ADDITIVE MULTI RESPONSE REGRESSION MODEL--------------
library(svMisc)
library(pspline)
library(splines)
library(psych)
library(splines2)
library(Rlab)
library(pracma)
library(psych)
stdf <- function(x){
  stdx <- (x-min(x))/(max(x)-min(x))
  return(stdx)
}
#---------------GCV-------------------------------------------------------------
gcvfunc<-function(y,yhat,p){
  y<-matrix(c(y))
  yhat<-matrix(c(yhat))
  n<-length(y)
  score<-(1/n)*(norm(y-yhat)^2)/((1/n)*p)
  return(score)
}
################################################################################
#SMOOTHING FUNCTION WITH LOCAL METHODS------------------------------------------
MRLP<-function(x1,x2,y1,y2,nc1,nc2){
  library(optR)
  library(condSURV)
  #x: matrix of parametric covariates
  #y: response variable(s) all
  #nc: Matrix of nonparametric covariate
  stdf <- function(x){
    stdx <- (x-min(x))/(max(x)-min(x))
    return(stdx)
  }
  n1 <- length(y1)
  n2 <- length(y2)
  r  <- 2
  p  <- 2
  #Prepare data-----------------------------------------------------------------
  t11 <- sort(stdf(nc1[,1]))
  t12 <- sort(stdf(nc1[,2]))
  t21 <- sort(stdf(nc2[,1]))
  t22 <- sort(stdf(nc2[,2]))
  y1  <- stdf(y1)
  y2  <- stdf(y2)
  x1  <- stdf(x1) 
  x2  <- stdf(x2)
  
  #-----------------------------------------------------------------------------
  #GCV function-------------------------------------------------------------------
  gcvfunc<-function(y,yhat,p){
    y<-matrix(c(y))
    yhat<-matrix(c(yhat))
    n<-length(y)
    score<-(1/n)*(norm(y-yhat)^2)/((1/n)*p)
    return(score)
  }
  #Preparation of estimation-----------------------------------------------
  ones1   <-matrix(1,n1,1) 
  ones2   <-matrix(1,n2,1) 
  W11     <-matrix(0,n1,n1)
  W21     <-matrix(0,n1,n1)
  W12     <-matrix(0,n2,n2)
  W22     <-matrix(0,n2,n2)
  GCV_g11 <- 0
  GCV_g21 <- 0
  GCV_g12 <- 0
  GCV_g22 <- 0
  
  LO     <-50 
  index11 <-seq(min(t11),max(t11),length.out=n1)
  index12 <-seq(min(t12),max(t12),length.out=n1)
  index21 <-seq(min(t21),max(t21),length.out=n2)
  index22 <-seq(min(t22),max(t22),length.out=n2)
  
  T11     <-matrix(c(ones1,(t11-index11),(t11-index11)^2),n1,3)
  T21     <-matrix(c(ones1,(t12-index12),(t12-index12)^2),n1,3)
  T12     <-matrix(c(ones2,(t21-index21),(t21-index21)^2),n2,3)
  T22     <-matrix(c(ones2,(t22-index22),(t22-index22)^2),n2,3)
  bw_seq <-seq(0.08,0.15,length.out=LO)
 
    for (i in 1:LO){
      for (j in 1:n1){
        W11[,j] <- NWW(t11,index11[j],bw=bw_seq[i])
        W21[,j] <- NWW(t12,index12[j],bw=bw_seq[i])
      }
      for (j in 1:n2){
        W12[,j] <- NWW(t21,index21[j],bw=bw_seq[i])
        W22[,j] <- NWW(t22,index22[j],bw=bw_seq[i])
      }
      
      
      S_g11  <- T11%*%solve(t(T11)%*%W11%*%T11,tol=1e-200)%*%t(T11)%*%W11
      S_g21  <- T21%*%solve(t(T21)%*%W21%*%T21,tol=1e-200)%*%t(T21)%*%W21
      
      S_g12  <- T12%*%solve(t(T12)%*%W12%*%T12,tol=1e-200)%*%t(T12)%*%W12
      S_g22  <- T22%*%solve(t(T22)%*%W22%*%T22,tol=1e-200)%*%t(T22)%*%W22
      
      xtil_g11  <- (diag(n1)-S_g11)%*%x1
      ytil_g11  <- (diag(n1)-S_g11)%*%y1
    
      xtil_g21  <- (diag(n1)-S_g21)%*%x1
      ytil_g21  <- (diag(n1)-S_g21)%*%y1
      
      xtil_g12  <- (diag(n2)-S_g12)%*%x2
      ytil_g12  <- (diag(n2)-S_g12)%*%y2
      
      xtil_g22  <- (diag(n2)-S_g22)%*%x2
      ytil_g22  <- (diag(n2)-S_g22)%*%y2
      
      betahat1_g1 <-solve(t(xtil_g11)%*%xtil_g11)%*%t(xtil_g11)%*%ytil_g11 
      betahat1_g2 <-solve(t(xtil_g21)%*%xtil_g21)%*%t(xtil_g21)%*%ytil_g21
      
      betahat2_g1 <-solve(t(xtil_g12)%*%xtil_g12)%*%t(xtil_g12)%*%ytil_g12 
      betahat2_g2 <-solve(t(xtil_g22)%*%xtil_g22)%*%t(xtil_g22)%*%ytil_g22
      
      beta1hat <- rowMeans(betahat1_g1,betahat1_g2)
      beta2hat <- rowMeans(betahat2_g1,betahat2_g2)
      
      g11 <- W11%*%(y1-x1%*%beta1hat)
      g21 <- W21%*%(y1-x1%*%beta1hat)
      
      g12 <- W12%*%(y2-x2%*%beta2hat)
      g22 <- W22%*%(y2-x2%*%beta2hat)
        
      library(psych)
      p0 <- tr(S_g11)
      p1 <- tr(S_g12)
      
      GCV_g11[i]  <-gcvfunc(y1,g11,p0)
      GCV_g21[i]  <-gcvfunc(y1,g21,p0)
      
      GCV_g12[i]  <-gcvfunc(y2,g12,p1)
      GCV_g22[i]  <-gcvfunc(y2,g22,p1)
    }
    ming <- min(scale(GCV_g11))
    maxg <- max(scale(GCV_g21))
    plot(bw_seq,scale(GCV_g11),type="l",ylab="GCV score",xlab="bandwidth",ylim=c(ming,maxg),col=2,lty=1,main="Bandwidth Selection for g(t1)")
    par(new=TRUE)
    plot(bw_seq,scale(GCV_g21),type="l",ylab="GCV score",xlab="bandwidth",ylim=c(ming,maxg),col=3,lty=2)
    grid()
    legend(x = "topright",legend = c("GCV(g1) for y1", "GCV(g1) for y2"), lty = c(1, 2),col = c(2, 3),lwd = 1)
    
    ming <- min(scale(GCV_g12))
    maxg <- max(scale(GCV_g22))
    plot(bw_seq,scale(GCV_g12),type="l",ylab="GCV score",xlab="bandwidth",ylim=c(ming,maxg),col=2,lty=1,main="Bandwidth Selection for g(t2)")
    par(new=TRUE)
    plot(bw_seq,scale(GCV_g22),type="l",ylab="GCV score",xlab="bandwidth",ylim=c(ming,maxg),col=3,lty=2)
    grid()
    legend(x = "topright",legend = c("GCV(g2) for y1", "GCV(g2) for y2"), lty = c(1, 2),col = c(2,3),lwd = 1)
    
    for (i2 in 1:LO){
      if (GCV_g11[i2]==min(GCV_g11)){
        bwg11 <- bw_seq[i2]
      }
      if (GCV_g21[i2]==min(GCV_g21)){
        bwg21 <- bw_seq[i2]
      }
      
      
      if (GCV_g12[i2]==min(GCV_g12)){
        bwg12 <- bw_seq[i2]
      }
      if (GCV_g22[i2]==min(GCV_g22)){
        bwg22 <- bw_seq[i2]
      }
    }
  #-------------------------------------------------------------------------------
  #ESTIMATION PART---------------------------------------------------------------
    ones1   <-matrix(1,n1,1) 
    ones2   <-matrix(1,n2,1) 
    W11     <-matrix(0,n1,n1)
    W21     <-matrix(0,n1,n1)
    W12     <-matrix(0,n2,n2)
    W22     <-matrix(0,n2,n2)
    
    index11 <-seq(min(t11)-0.1,max(t11)+0.1,length.out=n1)
    index12 <-seq(min(t12)-0.5,max(t12)+0.5,length.out=n1)
    index21 <-seq(min(t21)-0.1,max(t21)+0.1,length.out=n2)
    index22 <-seq(min(t22)-0.5,max(t22)+0.5,length.out=n2)
    
    T11     <-matrix(c(ones1,(t11-index11),(t11-index11)^2),n1,3)
    T21     <-matrix(c(ones1,(t12-index12),(t12-index12)^2),n1,3)
    T12     <-matrix(c(ones2,(t21-index21),(t21-index21)^2),n2,3)
    T22     <-matrix(c(ones2,(t22-index22),(t22-index22)^2),n2,3)
    
    
 
    for (j in 1:n1){
      W11[,j] <- NWW(t11,index11[j],bw=bwg11)
      W21[,j] <- NWW(t12,index12[j],bw=bwg21)
    }
    for (j in 1:n2){
      W12[,j] <- NWW(t21,index21[j],bw=bwg12)
      W22[,j] <- NWW(t22,index22[j],bw=bwg22)
    }
    S_g11  <- T11%*%solve(t(T11)%*%W11%*%T11)%*%t(T11)%*%W11
    S_g21  <- T21%*%solve(t(T21)%*%W21%*%T21)%*%t(T21)%*%W21
    
    S_g12  <- T12%*%solve(t(T12)%*%W12%*%T12)%*%t(T12)%*%W12
    S_g22  <- T22%*%solve(t(T22)%*%W22%*%T22)%*%t(T22)%*%W22
    
    xtil_g11  <- (diag(n1)-S_g11)%*%x1
    ytil_g11  <- (diag(n1)-S_g11)%*%y1
    
    xtil_g21  <- (diag(n1)-S_g21)%*%x1
    ytil_g21  <- (diag(n1)-S_g21)%*%y1
    
    xtil_g12  <- (diag(n2)-S_g12)%*%x2
    ytil_g12  <- (diag(n2)-S_g12)%*%y2
    
    xtil_g22  <- (diag(n2)-S_g22)%*%x2
    ytil_g22  <- (diag(n2)-S_g22)%*%y2
    
    betahat1_g1 <-solve(t(xtil_g11)%*%xtil_g11)%*%t(xtil_g11)%*%ytil_g11 
    betahat1_g2 <-solve(t(xtil_g21)%*%xtil_g21)%*%t(xtil_g21)%*%ytil_g21
    
    betahat2_g1 <-solve(t(xtil_g12)%*%xtil_g12)%*%t(xtil_g12)%*%ytil_g12 
    betahat2_g2 <-solve(t(xtil_g22)%*%xtil_g22)%*%t(xtil_g22)%*%ytil_g22
    
    beta1hat <- rowMeans(betahat1_g1,betahat1_g2)
    beta2hat <- rowMeans(betahat2_g1,betahat2_g2)
    
    g11 <- W11%*%(y1-x1%*%beta1hat)
    g21 <- W21%*%(y1-x1%*%beta1hat-g11)
    
    g12 <- W12%*%(y2-x2%*%beta2hat)
    g22 <- W22%*%(y2-x2%*%beta2hat-g12)
    
    
    y1hat <- x1%*%beta1hat+g11+g21
    y2hat <- x2%*%beta2hat+g12+g22
   
    
    res <- new.env()
    
    res$beta1hat <- beta1hat
    res$beta2hat <- beta2hat
    #------------------------------------------------------------------------------    
    res$g11 <- g11
    res$g21 <- g21
    
    res$g12 <- g12
    res$g22 <- g22
    #------------------------------------------------------------------------------
    res$y1hat <- y1hat
    res$y2hat <- y2hat
    res$t11 <- index11
    res$t12 <- index12
    res$t21 <- index21
    res$t22 <- index22
  
  #-------------------------------------------------------------------------------
  return(res)
  #-----------------------------------------------------------------------------
  #-----------------------------------------------------------------------------
}
MRLL<-function(x1,x2,y1,y2,nc1,nc2){
  library(optR)
  library(condSURV)
  #x: matrix of parametric covariates
  #y: response variable(s) all
  #nc: Matrix of nonparametric covariate
  #Prepare data-------------------------------------------------------------------
  stdf <- function(x){
    stdx <- (x-min(x))/(max(x)-min(x))
    return(stdx)
  }
  n1 <- length(y1)
  n2 <- length(y2)
  r  <- 2
  p  <- 2
  #Prepare data-----------------------------------------------------------------
  t11 <- sort(stdf(nc1[,1]))
  t12 <- sort(stdf(nc1[,2]))
  t21 <- sort(stdf(nc2[,1]))
  t22 <- sort(stdf(nc2[,2]))
  y1  <- stdf(y1)
  y2  <- stdf(y2)
  x1  <- stdf(x1) 
  x2  <- stdf(x2)
  #-----------------------------------------------------------------------------
  
  #-------------------------------------------------------------------------------
  #GCV function-------------------------------------------------------------------
  gcvfunc<-function(y,yhat,p){
    y<-matrix(c(y))
    yhat<-matrix(c(yhat))
    n<-length(y)
    score<-(1/n)*(norm(y-yhat)^2)/((1/n)*p)
    return(score)
  }
  #Preparation of estimation-----------------------------------------------
  ones1   <-matrix(1,n1,1) 
  ones2   <-matrix(1,n2,1) 
  W11     <-matrix(0,n1,n1)
  W21     <-matrix(0,n1,n1)
  W12     <-matrix(0,n2,n2)
  W22     <-matrix(0,n2,n2)
  GCV_g11 <- 0
  GCV_g21 <- 0
  GCV_g12 <- 0
  GCV_g22 <- 0
  
  LO     <-50 
  index11 <-seq(min(t11)-0.1,max(t11)+0.1,length.out=n1)
  index12 <-seq(min(t12)-0.5,max(t12)+0.5,length.out=n1)
  index21 <-seq(min(t21)-0.1,max(t21)+0.1,length.out=n2)
  index22 <-seq(min(t22)-0.5,max(t22)+0.5,length.out=n2)
  
  T11     <-matrix(c(ones1,(t11-index11)),n1,2)
  T21     <-matrix(c(ones1,(t12-index12)),n1,2)
  T12     <-matrix(c(ones2,(t21-index21)),n2,2)
  T22     <-matrix(c(ones2,(t22-index22)),n2,2)
  bw_seq <-seq(0.08,0.15,length.out=LO)
  
  for (i in 1:LO){
    for (j in 1:n1){
      W11[,j] <- NWW(t11,index11[j],bw=bw_seq[i])
      W21[,j] <- NWW(t12,index12[j],bw=bw_seq[i])
    }
    for (j in 1:n2){
      W12[,j] <- NWW(t21,index21[j],bw=bw_seq[i])
      W22[,j] <- NWW(t22,index22[j],bw=bw_seq[i])
    }
    
    S_g11  <- T11%*%solve(t(T11)%*%W11%*%T11)%*%t(T11)%*%W11
    S_g21  <- T21%*%solve(t(T21)%*%W21%*%T21)%*%t(T21)%*%W21
    
    S_g12  <- T12%*%solve(t(T12)%*%W12%*%T12)%*%t(T12)%*%W12
    S_g22  <- T22%*%solve(t(T22)%*%W22%*%T22)%*%t(T22)%*%W22
    
    xtil_g11  <- (diag(n1)-S_g11)%*%x1
    ytil_g11  <- (diag(n1)-S_g11)%*%y1
    
    xtil_g21  <- (diag(n1)-S_g21)%*%x1
    ytil_g21  <- (diag(n1)-S_g21)%*%y1
    
    xtil_g12  <- (diag(n2)-S_g12)%*%x2
    ytil_g12  <- (diag(n2)-S_g12)%*%y2

    xtil_g22  <- (diag(n2)-S_g22)%*%x2
    ytil_g22  <- (diag(n2)-S_g22)%*%y2
    
    betahat1_g1 <-solve(t(xtil_g11)%*%xtil_g11)%*%t(xtil_g11)%*%ytil_g11 
    betahat1_g2 <-solve(t(xtil_g21)%*%xtil_g21)%*%t(xtil_g21)%*%ytil_g21
    
    betahat2_g1 <-solve(t(xtil_g12)%*%xtil_g12)%*%t(xtil_g12)%*%ytil_g12 
    betahat2_g2 <-solve(t(xtil_g22)%*%xtil_g22)%*%t(xtil_g22)%*%ytil_g22
    
    beta1hat <- rowMeans(betahat1_g1,betahat1_g2)
    beta2hat <- rowMeans(betahat2_g1,betahat2_g2)
    
    g11 <- W11%*%(y1-x1%*%beta1hat)
    g21 <- W21%*%(y1-x1%*%beta1hat)
    g11 <- W11%*%(y1-x1%*%beta1hat-g21)
    g21 <- W21%*%(y1-x1%*%beta1hat-g11)
    
    g12 <- W12%*%(y2-x2%*%beta2hat)
    g22 <- W22%*%(y2-x2%*%beta2hat-g12)
    g12 <- W12%*%(y2-x2%*%beta2hat-g22)
    g22 <- W22%*%(y2-x2%*%beta2hat-g12)
    
    library(psych)
    p0 <- tr(S_g11)
    p1 <- tr(S_g12)
    
    GCV_g11[i]  <-gcvfunc(y1,g11,p0)
    GCV_g21[i]  <-gcvfunc(y1,g21,p0)
    
    GCV_g12[i]  <-gcvfunc(y2,g12,p1)
    GCV_g22[i]  <-gcvfunc(y2,g22,p1)
  }
  ming <- min(scale(GCV_g11))
  maxg <- max(scale(GCV_g21))
  plot(bw_seq,scale(GCV_g11),type="l",ylab="GCV score",xlab="bandwidth",ylim=c(ming,maxg),col=2,lty=1,main="Bandwidth Selection for g(t1)")
  par(new=TRUE)
  plot(bw_seq,scale(GCV_g21),type="l",ylab="GCV score",xlab="bandwidth",ylim=c(ming,maxg),col=3,lty=2)
  grid()
  legend(x = "topright",legend = c("GCV(g1) for y1", "GCV(g1) for y2"), lty = c(1, 2),col = c(2, 3),lwd = 1)
  
  ming <- min(scale(GCV_g12))
  maxg <- max(scale(GCV_g22))
  plot(bw_seq,scale(GCV_g12),type="l",ylab="GCV score",xlab="bandwidth",ylim=c(ming,maxg),col=2,lty=1,main="Bandwidth Selection for g(t2)")
  par(new=TRUE)
  plot(bw_seq,scale(GCV_g22),type="l",ylab="GCV score",xlab="bandwidth",ylim=c(ming,maxg),col=3,lty=2)
  grid()
  legend(x = "topright",legend = c("GCV(g2) for y1", "GCV(g2) for y2"), lty = c(1, 2),col = c(2,3),lwd = 1)
  
  for (i2 in 1:LO){
    if (GCV_g11[i2]==min(GCV_g11)){
      bwg11 <- bw_seq[i2]
    }
    if (GCV_g21[i2]==min(GCV_g21)){
      bwg21 <- bw_seq[i2]
    }
    
    
    if (GCV_g12[i2]==min(GCV_g12)){
      bwg12 <- bw_seq[i2]
    }
    if (GCV_g22[i2]==min(GCV_g22)){
      bwg22 <- bw_seq[i2]
    }
  }
  #-------------------------------------------------------------------------------
  #ESTIMATION PART---------------------------------------------------------------
  ones1   <-matrix(1,n1,1) 
  ones2   <-matrix(1,n2,1) 
  W11     <-matrix(0,n1,n1)
  W21     <-matrix(0,n1,n1)
  W12     <-matrix(0,n2,n2)
  W22     <-matrix(0,n2,n2)
  
  index11 <-seq(min(t11)-0.1,max(t11)+0.1,length.out=n1)
  index12 <-seq(min(t12)-0.5,max(t12)+0.5,length.out=n1)
  index21 <-seq(min(t21)-0.1,max(t21)+0.1,length.out=n2)
  index22 <-seq(min(t22)-0.5,max(t22)+0.5,length.out=n2)
  
  T11     <-matrix(c(ones1,(t11-index11)),n1,2)
  T21     <-matrix(c(ones1,(t12-index12)),n1,2)
  T12     <-matrix(c(ones2,(t21-index21)),n2,2)
  T22     <-matrix(c(ones2,(t22-index22)),n2,2)
  
  
  for (j in 1:n1){
    W11[,j] <- NWW(t11,index11[j],bw=bwg11)
    W21[,j] <- NWW(t12,index12[j],bw=bwg21)
  }
  for (j in 1:n2){
    W12[,j] <- NWW(t21,index21[j],bw=bwg12)
    W22[,j] <- NWW(t22,index22[j],bw=bwg22)
  }
  S_g11  <- T11%*%solve(t(T11)%*%W11%*%T11)%*%t(T11)%*%W11
  S_g21  <- T21%*%solve(t(T21)%*%W21%*%T21)%*%t(T21)%*%W21
  
  S_g12  <- T12%*%solve(t(T12)%*%W12%*%T12)%*%t(T12)%*%W12
  S_g22  <- T22%*%solve(t(T22)%*%W22%*%T22)%*%t(T22)%*%W22
  
  xtil_g11  <- (diag(n1)-S_g11)%*%x1
  ytil_g11  <- (diag(n1)-S_g11)%*%y1
  
  xtil_g21  <- (diag(n1)-S_g21)%*%x1
  ytil_g21  <- (diag(n1)-S_g21)%*%y1
  
  xtil_g12  <- (diag(n2)-S_g12)%*%x2
  ytil_g12  <- (diag(n2)-S_g12)%*%y2
  
  xtil_g22  <- (diag(n2)-S_g22)%*%x2
  ytil_g22  <- (diag(n2)-S_g22)%*%y2
  
  betahat1_g1 <-solve(t(xtil_g11)%*%xtil_g11)%*%t(xtil_g11)%*%ytil_g11 
  betahat1_g2 <-solve(t(xtil_g21)%*%xtil_g21)%*%t(xtil_g21)%*%ytil_g21
  
  betahat2_g1 <-solve(t(xtil_g12)%*%xtil_g12)%*%t(xtil_g12)%*%ytil_g12 
  betahat2_g2 <-solve(t(xtil_g22)%*%xtil_g22)%*%t(xtil_g22)%*%ytil_g22
  
  beta1hat <- rowMeans(betahat1_g1,betahat1_g2)
  beta2hat <- rowMeans(betahat2_g1,betahat2_g2)
  
  g11 <- W11%*%(y1-x1%*%beta1hat)
  g21 <- W21%*%(y1-x1%*%beta1hat-g11)
  
  g12 <- W12%*%(y2-x2%*%beta2hat)
  g22 <- W22%*%(y2-x2%*%beta2hat-g12)
  
  y1hat <- x1%*%beta1hat+g11+g21
  y2hat <- x2%*%beta2hat+g12+g22
  
  
  res <- new.env()
  
  res$beta1hat <- beta1hat
  res$beta2hat <- beta2hat
  #------------------------------------------------------------------------------    
  res$g11 <- g11
  res$g21 <- g21
  
  res$g12 <- g12
  res$g22 <- g22
  #------------------------------------------------------------------------------
  res$y1hat <- y1hat
  res$y2hat <- y2hat
  
  res$t11 <- index11
  res$t12 <- index12
  res$t21 <- index21
  res$t22 <- index22
  
  
  #-------------------------------------------------------------------------------
  return(res)
  #-----------------------------------------------------------------------------
  #-----------------------------------------------------------------------------
}
MRLC<-function(x1,x2,y1,y2,nc1,nc2){
  library(optR)
  library(condSURV)
  #x: matrix of parametric covariates
  #y: response variable(s) all
  #nc: Matrix of nonparametric covariate
  #Prepare data-------------------------------------------------------------------
  stdf <- function(x){
    stdx <- (x-min(x))/(max(x)-min(x))
    return(stdx)
  }
  n1 <- length(y1)
  n2 <- length(y2)
  r  <- 2
  p  <- 2
  #Prepare data-----------------------------------------------------------------
  t11 <- sort(stdf(nc1[,1]))
  t12 <- sort(stdf(nc1[,2]))
  t21 <- sort(stdf(nc2[,1]))
  t22 <- sort(stdf(nc2[,2]))
  y1  <- stdf(y1)
  y2  <- stdf(y2)
  x1  <- stdf(x1) 
  x2  <- stdf(x2)
  #-----------------------------------------------------------------------------
  
  #-------------------------------------------------------------------------------
  #GCV function-------------------------------------------------------------------
  gcvfunc<-function(y,yhat,p){
    y<-matrix(c(y))
    yhat<-matrix(c(yhat))
    n<-length(y)
    score<-(1/n)*(norm(y-yhat)^2)/((1/n)*p)
    return(score)
  }
  #Preparation of estimation-----------------------------------------------
  ones1   <-matrix(1,n1,1) 
  ones2   <-matrix(1,n2,1) 
  W11     <-matrix(0,n1,n1)
  W21     <-matrix(0,n1,n1)
  W12     <-matrix(0,n2,n2)
  W22     <-matrix(0,n2,n2)
  GCV_g11 <- 0
  GCV_g21 <- 0
  GCV_g12 <- 0
  GCV_g22 <- 0
  
  LO     <-50 
  index11 <-seq(min(t11)-0.1,max(t11)+0.1,length.out=n1)
  index12 <-seq(min(t12)-0.5,max(t12)+0.5,length.out=n1)
  index21 <-seq(min(t21)-0.1,max(t21)+0.1,length.out=n2)
  index22 <-seq(min(t22)-0.5,max(t22)+0.5,length.out=n2)
  
  T11     <-matrix(c(ones1,(t11-index11),n1,2))
  T21     <-matrix(c(ones1,(t12-index12),n1,2))
  T12     <-matrix(c(ones2,(t21-index21),n2,2))
  T22     <-matrix(c(ones2,(t22-index22),n2,2))
  bw_seq <-seq(0.08,0.15,length.out=LO)
  
  for (i in 1:LO){
    for (j in 1:n1){
      W11[,j] <- NWW(t11,index11[j],bw=bw_seq[i])
      W21[,j] <- NWW(t12,index12[j],bw=bw_seq[i])
    }
    for (j in 1:n2){
      W12[,j] <- NWW(t21,index21[j],bw=bw_seq[i])
      W22[,j] <- NWW(t22,index22[j],bw=bw_seq[i])
    }
    
    S_g11  <- W11
    S_g21  <- W21
    
    S_g12  <- W12
    S_g22  <- W22
    
    xtil_g11  <- (diag(n1)-S_g11)%*%x1
    ytil_g11  <- (diag(n1)-S_g11)%*%y1
    
    xtil_g21  <- (diag(n1)-S_g21)%*%x1
    ytil_g21  <- (diag(n1)-S_g21)%*%y1
    
    xtil_g12  <- (diag(n2)-S_g12)%*%x2
    ytil_g12  <- (diag(n2)-S_g12)%*%y2
    
    xtil_g22  <- (diag(n2)-S_g22)%*%x2
    ytil_g22  <- (diag(n2)-S_g22)%*%y2
    
    betahat1_g1 <-solve(t(xtil_g11)%*%xtil_g11)%*%t(xtil_g11)%*%ytil_g11 
    betahat1_g2 <-solve(t(xtil_g21)%*%xtil_g21)%*%t(xtil_g21)%*%ytil_g21
    
    betahat2_g1 <-solve(t(xtil_g12)%*%xtil_g12)%*%t(xtil_g12)%*%ytil_g12 
    betahat2_g2 <-solve(t(xtil_g22)%*%xtil_g22)%*%t(xtil_g22)%*%ytil_g22
    
    beta1hat <- rowMeans(betahat1_g1,betahat1_g2)
    beta2hat <- rowMeans(betahat2_g1,betahat2_g2)
    
    g11 <- W11%*%(y1-x1%*%beta1hat)
    g21 <- W21%*%(y1-x1%*%beta1hat)
    
    g12 <- W12%*%(y2-x2%*%beta2hat)
    g22 <- W22%*%(y2-x2%*%beta2hat)
    
    library(psych)
    p0 <- tr(S_g11)
    p1 <- tr(S_g12)
    
    GCV_g11[i]  <-gcvfunc(y1,g11,p0)
    GCV_g21[i]  <-gcvfunc(y1,g21,p0)
    
    GCV_g12[i]  <-gcvfunc(y2,g12,p1)
    GCV_g22[i]  <-gcvfunc(y2,g22,p1)
  }
  ming <- min(scale(GCV_g11))
  maxg <- max(scale(GCV_g21))
  plot(bw_seq,scale(GCV_g11),type="l",ylab="GCV score",xlab="bandwidth",ylim=c(ming,maxg),col=2,lty=1,main="Bandwidth Selection for g(t1)")
  par(new=TRUE)
  plot(bw_seq,scale(GCV_g21),type="l",ylab="GCV score",xlab="bandwidth",ylim=c(ming,maxg),col=3,lty=2)
  grid()
  legend(x = "topright",legend = c("GCV(g1) for y1", "GCV(g1) for y2"), lty = c(1, 2),col = c(2, 3),lwd = 1)
  
  ming <- min(scale(GCV_g12))
  maxg <- max(scale(GCV_g22))
  plot(bw_seq,scale(GCV_g12),type="l",ylab="GCV score",xlab="bandwidth",ylim=c(ming,maxg),col=2,lty=1,main="Bandwidth Selection for g(t2)")
  par(new=TRUE)
  plot(bw_seq,scale(GCV_g22),type="l",ylab="GCV score",xlab="bandwidth",ylim=c(ming,maxg),col=3,lty=2)
  grid()
  legend(x = "topright",legend = c("GCV(g2) for y1", "GCV(g2) for y2"), lty = c(1, 2),col = c(2,3),lwd = 1)
  
  for (i2 in 1:LO){
    if (GCV_g11[i2]==min(GCV_g11)){
      bwg11 <- bw_seq[i2]
    }
    if (GCV_g21[i2]==min(GCV_g21)){
      bwg21 <- bw_seq[i2]
    }
    
    
    if (GCV_g12[i2]==min(GCV_g12)){
      bwg12 <- bw_seq[i2]
    }
    if (GCV_g22[i2]==min(GCV_g22)){
      bwg22 <- bw_seq[i2]
    }
  }
  #-------------------------------------------------------------------------------
  #ESTIMATION PART---------------------------------------------------------------
  ones1   <-matrix(1,n1,1) 
  ones2   <-matrix(1,n2,1) 
  W11     <-matrix(0,n1,n1)
  W21     <-matrix(0,n1,n1)
  W12     <-matrix(0,n2,n2)
  W22     <-matrix(0,n2,n2)
  
  index11 <-seq(min(t11)-0.1,max(t11)+0.1,length.out=n1)
  index12 <-seq(min(t12)-0.5,max(t12)+0.5,length.out=n1)
  index21 <-seq(min(t21)-0.1,max(t21)+0.1,length.out=n2)
  index22 <-seq(min(t22)-0.5,max(t22)+0.5,length.out=n2)
  
  T11     <-matrix(c(ones1,(t11-index11),n1,2))
  T21     <-matrix(c(ones1,(t12-index12),n1,2))
  T12     <-matrix(c(ones2,(t21-index21),n2,2))
  T22     <-matrix(c(ones2,(t22-index22),n2,2))
  
  
  for (j in 1:n1){
    W11[,j] <- NWW(t11,index11[j],bw=bwg11)
    W21[,j] <- NWW(t12,index12[j],bw=bwg21)
  }
  for (j in 1:n2){
    W12[,j] <- NWW(t21,index21[j],bw=bwg12)
    W22[,j] <- NWW(t22,index22[j],bw=bwg22)
  }
  
  S_g11  <- W11
  S_g21  <- W21
  
  S_g12  <- W12
  S_g22  <- W22
  
  xtil_g11  <- (diag(n1)-S_g11)%*%x1
  ytil_g11  <- (diag(n1)-S_g11)%*%y1
  
  xtil_g21  <- (diag(n1)-S_g21)%*%x1
  ytil_g21  <- (diag(n1)-S_g21)%*%y1
  
  xtil_g12  <- (diag(n2)-S_g12)%*%x2
  ytil_g12  <- (diag(n2)-S_g12)%*%y2
  
  xtil_g22  <- (diag(n2)-S_g22)%*%x2
  ytil_g22  <- (diag(n2)-S_g22)%*%y2
  
  betahat1_g1 <-solve(t(xtil_g11)%*%xtil_g11)%*%t(xtil_g11)%*%ytil_g11 
  betahat1_g2 <-solve(t(xtil_g21)%*%xtil_g21)%*%t(xtil_g21)%*%ytil_g21
  
  betahat2_g1 <-solve(t(xtil_g12)%*%xtil_g12)%*%t(xtil_g12)%*%ytil_g12 
  betahat2_g2 <-solve(t(xtil_g22)%*%xtil_g22)%*%t(xtil_g22)%*%ytil_g22
  
  beta1hat <- rowMeans(betahat1_g1,betahat1_g2)
  beta2hat <- rowMeans(betahat2_g1,betahat2_g2)
  
  g11 <- W11%*%(y1-x1%*%beta1hat)
  g21 <- W21%*%(y1-x1%*%beta1hat-g11)
  
  g12 <- W12%*%(y2-x2%*%beta2hat)
  g22 <- W22%*%(y2-x2%*%beta2hat-g12)
  
  
  y1hat <- x1%*%beta1hat+g11+g21
  y2hat <- x2%*%beta2hat+g12+g22
  
  
  res <- new.env()
  
  res$beta1hat <- beta1hat
  res$beta2hat <- beta2hat
  #------------------------------------------------------------------------------    
  res$g11 <- g11
  res$g21 <- g21
  
  res$g12 <- g12
  res$g22 <- g22
  #------------------------------------------------------------------------------
  res$y1hat <- y1hat
  res$y2hat <- y2hat
  
  res$t11 <- index11
  res$t12 <- index12
  res$t21 <- index21
  res$t22 <- index22
  
  #-------------------------------------------------------------------------------
  return(res)
  #-----------------------------------------------------------------------------
  #-----------------------------------------------------------------------------
}
syndata<-function(y,delta){
  library(pracma)
  # where y: right-censored observations, delta: censorship indicator
  n<-length(y)
  M<-0
  yg<-0
  M<-0
  delta2<-0
  #Synthetic data transformation
  y1<-cbind(y,delta)
  y2<-sortrows(y1,k=2)
  delta2<-y2[,2]
  delta2[n]<-1
  sy<-sort(y)
  for (i1 in 1:n){
    Ma<-1
    for (i2 in 1:n){
      mGt<-((n-i2)/(n-i2+1))^((delta2[i2]==0)*(sy[i2]<=sy[i1]))
      Ma<-Ma*mGt
    }
    M[i1]=Ma
    yg[i1]<-(delta[i1]*y[i1])/M[i1]
  }
  return(yg)
}
#-------------------------------------------------------------------------------

data  <- read.table("cancer.txt",header=TRUE)
data <- as.matrix(data)
n     <- length(data[,1])
data1 <- matrix(0,n,7)
data2 <- matrix(0,n,7)
a <- 1
b <- 1
for (i in 1:n){
  if (data[i,1]>=61 & data[i,1]<=71){
    data1[a,] <-data[i,] 
    a <- a+1
  }
  else{
    data2[b,] <-data[i,]
    b <- b+1
  }
}
#-------------------------------------------------------------------------------
data1           <- as.matrix(data.frame(data1[1:(a-1),]))
data1[,4]       <-data1[,4]/100
colnames(data1) <- colnames(data)
data2           <- as.matrix(data.frame(data2[1:(b-1),]))
data2[,4]       <- data2[,4]/100
colnames(data2) <- colnames(data)
data1 <- data.frame(data1)
data2 <- data.frame(data2)
n1              <-a-1
n2              <-b-1 
#FOR DATA 1---------------------------------------------------------------------
x1     <- matrix(c(scale(data1$optime),data1$tumorlev),n1,2)
delta1 <- data1$delta
y1     <- syndata(scale(data1$stime),delta1) 
t11    <- (scale(data1$alb))
t12    <- (scale(data1$tblood)+rnorm(n1,sd=0.1))
nc1    <- matrix(c(t11,t12),n1,2)
CL1    <- (n1-sum(delta1))/n1
#FOR DATA 2---------------------------------------------------------------------
x2     <- matrix(c(scale(data2$optime),data2$tumorlev),n2,2)
delta2 <- data2$delta
y2     <- syndata(scale(data2$stime),delta2) 
t21    <- (scale(data2$alb))
t22    <- (scale(data2$tblood)+rnorm(n2,sd=0.1))
nc2    <- matrix(c(t21,t22),n2,2)
CL2    <- (n2-sum(delta2))/n2
#-------------------------------------------------------------------------------
library(GGally)
ggpairs(data1[,2:6])
ggpairs(data2[,2:6])
summary(data1[,2:6])
summary(data2[,2:6])
#-------------------------------------------------------------------------------
LPobj <- MRLP(x1,x2,y1,y2,nc1,nc2)
LLobj <- MRLL(x1,x2,y1,y2,nc1,nc2)
LCobj <- MRLC(x1,x2,y1,y2,nc1,nc2)

LP.var1        <- var(y1-LPobj$y1hat)
LP.var2        <- var(y2-LPobj$y2hat)
LP.RMSE.model1 <- sqrt(mean((y1-LPobj$y1hat)^2))
LP.RMSE.model2 <- sqrt(mean((y2-LPobj$y2hat)^2))
data.frame(LP.var1,LP.var2,LP.RMSE.model1,LP.RMSE.model2)

LP.MSE1  <- mean(((y1-x1%*%LPobj$beta1hat)-LPobj$g11)^2)
LP.MSE2  <- mean(((y2-x2%*%LPobj$beta2hat)-LPobj$g12)^2)
gall1    <- LPobj$g11+LPobj$g21 
gall2    <- LPobj$g22+LPobj$g12 
LP.TMSE1 <- mean(((y1-x1%*%LPobj$beta1hat)-gall1)^2)
LP.TMSE2 <- mean(((y2-x2%*%LPobj$beta2hat)-gall2)^2)

data.frame(LP.MSE1,LP.MSE2,LP.TMSE1,LP.TMSE2)

LL.var1        <- var(y1-LLobj$y1hat)
LL.var2        <- var(y2-LLobj$y2hat)
LL.RMSE.model1 <- sqrt(mean((y1-LLobj$y1hat)^2))
LL.RMSE.model2 <- sqrt(mean((y2-LLobj$y2hat)^2))

data.frame(LL.var1,LL.var2,LL.RMSE.model1,LL.RMSE.model2)

LL.MSE1  <- mean(((y1-x1%*%LLobj$beta1hat)-LLobj$g11)^2)
LL.MSE2  <- mean(((y2-x2%*%LLobj$beta2hat)-LLobj$g12)^2)
LL.gall1    <- LLobj$g11+LLobj$g21 
LL.gall2    <- LLobj$g22+LLobj$g12 
LL.TMSE1 <- mean(((y1-x1%*%LLobj$beta1hat)-LL.gall1)^2)
LL.TMSE2 <- mean(((y2-x2%*%LLobj$beta2hat)-LL.gall2)^2)

data.frame(LL.MSE1,LL.MSE2,LL.TMSE1,LL.TMSE2)

LC.var1        <- var(y1-LCobj$y1hat)
LC.var2        <- var(y2-LCobj$y2hat)
LC.RMSE.model1 <- sqrt(mean((y1-LCobj$y1hat)^2))
LC.RMSE.model2 <- sqrt(mean((y2-LCobj$y2hat)^2))
data.frame(LC.var1,LC.var2,LC.RMSE.model1,LC.RMSE.model2)

LC.MSE1  <- mean(((y1-x1%*%LCobj$beta1hat)-LCobj$g11)^2)
LC.MSE2  <- mean(((y2-x2%*%LCobj$beta2hat)-LCobj$g12)^2)
LC.gall1 <- LCobj$g11+LCobj$g21 
LC.gall2 <- LCobj$g22+LCobj$g12 
LC.TMSE1 <- mean(((y1-x1%*%LCobj$beta1hat)-LC.gall1)^2)
LC.TMSE2 <- mean(((y2-x2%*%LCobj$beta2hat)-LC.gall2)^2)

data.frame(LC.MSE1,LC.MSE2,LC.TMSE1,LC.TMSE2)

LPobj$beta1hat
LPobj$beta2hat

LLobj$beta1hat
LLobj$beta2hat

LCobj$beta1hat
LCobj$beta2hat

#-------------------------------------------------------------------------------
#This plot for the g1(t1) for y1 model (stime_1)
par(mar=c(4,4,2,1))
points1 <- stdf(y1-(x1%*%LPobj$beta1hat-LPobj$g21))        
ming1 <- min(points1)
maxg1 <- max(points1)
plot(t11,points1,pch=19,ylim=c(ming1,maxg1),xlab="Alb",ylab="g(Alb)",cex=0.7,main="(a)")
par(new=TRUE)
plot(sort(t11),LPobj$g11,type="l",ylim=c(ming1,maxg1),xlab="Alb",ylab="g(Alb)",col=1,lty=1,lwd=2)
par(new=TRUE)
plot(sort(t11),LLobj$g11,type="l",ylim=c(ming1,maxg1),xlab="Alb",ylab="g(Alb)",col=2,lty=2,lwd=2)
par(new=TRUE)
plot(sort(t11),LCobj$g11,type="l",ylim=c(ming1,maxg1),xlab="Alb",ylab="g(Alb)",col=3,lty=3,lwd=2)
grid()
legend(x = "topright",legend = c("LP.g(Alb)","LL.g(Alb)","LC.g(Alb)"), lty = c(1, 2,3),col = c(1, 2,3),lwd = 2)
par(new=FALSE)
#This plot for the g2(t2) for y1 model (stime_1)
points1.2 <- stdf(y1-(x1%*%LPobj$beta1hat-LPobj$g11))
ming1.2 <- min(LPobj$g21+0.4)
maxg1.2 <- max(LPobj$g21+0.4)
plot(t12,points1.2,pch=19,ylim=c(ming1.2,maxg1.2),xlim=c(-1,2),xlab="tx",ylab="g(tx)",cex=0.7,main="(b)")
par(new=TRUE)
plot(sort(t12),LPobj$g21+0.2,type="l",ylim=c(ming1.2,maxg1.2),xlim=c(-1,2),xlab="tx",ylab="g(tx)",col=1,lty=1,lwd=2)
par(new=TRUE)
plot(sort(t12),LLobj$g21+0.2,type="l",ylim=c(ming1.2,maxg1.2),xlim=c(-1,2),xlab="tx",ylab="g(tx)",col=2,lty=2,lwd=2)
par(new=TRUE)
plot(sort(t12),LCobj$g21+0.3,type="l",ylim=c(ming1.2,maxg1.2),xlim=c(-1,2),xlab="tx",ylab="g(tx)",col=3,lty=3,lwd=2)
grid()
legend(x = "topright",legend = c("LP.g(tx)","LL.g(tx)","LC.g(tx)"), lty = c(1, 2,3),col = c(1, 2,3),lwd = 2)
par(new=FALSE)

#-------------------------------------------------------------------------------
#This plot for the g1(t1) for y2 model (stime_2)
points2 <- stdf(y2-(x2%*%LPobj$beta2hat-LPobj$g22))
ming2 <- min(points2)
maxg2 <- max(points2)
plot(t21,points2,pch=19,ylim=c(ming2,maxg2),xlab="Alb",ylab="g(Alb)",cex=0.7,main="(c)")
par(new=TRUE)
plot(sort(t21) ,LPobj$g12,type="l",ylim=c(ming2,maxg2),xlab="Alb",ylab="g(Alb)",col=1,lty=1,lwd=2)
par(new=TRUE)
plot(sort(t21),LLobj$g12,type="l",ylim=c(ming2,maxg2),xlab="Alb",ylab="g(Alb)",col=2,lty=2,lwd=2)
par(new=TRUE)
plot(sort(t21),LCobj$g12,type="l",ylim=c(ming2,maxg2),xlab="Alb",ylab="g(Alb)",col=3,lty=3,lwd=2)
grid()
legend(x = "bottomleft",legend = c("LP.g(Alb)","LL.g(Alb)","LC.g(Alb)"), lty = c(1, 2,3),col = c(1, 2,3),lwd = 2)
par(new=FALSE)
#This plot for the g2(t2) for y2 model (stime_2)

points2.2 <- stdf(y2-(x2%*%LPobj$beta2hat)-LPobj$g12)
ming2.2 <- min(points2.2)
maxg2.2 <- max(points2.2)
plot(t22,points2.2,pch=19,ylim=c(ming2.2,maxg2.2),xlab="tx",ylab="g(tx)",cex=0.7,main="(d)")
par(new=TRUE)
plot(sort(t22),LPobj$g22+0.4,type="l",ylim=c(ming2.2,maxg2.2),xlab="tx",ylab="g(tx)",col=1,lty=1,lwd=2)
par(new=TRUE)
plot(sort(t22),LLobj$g22+0.4,type="l",ylim=c(ming2.2,maxg2.2),xlab="tx",ylab="g(tx)",col=2,lty=2,lwd=2)
par(new=TRUE)
plot(sort(t22),LCobj$g22+0.4,type="l",ylim=c(ming2.2,maxg2.2),xlab="tx",ylab="g(tx)",col=3,lty=3,lwd=2)
grid()
legend(x = "topright",legend = c("LP.g(tx)","LL.g(tx)","LC.g(tx)"), lty = c(1,2,3),col = c(1,2,3),lwd = 2)
par(new=FALSE)



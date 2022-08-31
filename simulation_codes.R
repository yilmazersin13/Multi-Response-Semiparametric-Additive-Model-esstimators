#SIMULATION STUDY FOR THE ADDITIVE MULTI RESPONSE REGRESSION MODEL--------------
library(svMisc)
library(pspline)
library(splines)
library(psych)
library(splines2)
library(Rlab)
library(pracma)
library(psych)
#Data Generation----------------------------------------------------------------
#DATA GENERATION----------------------------------------------------------------
#n=(30, 50, 100)
datagen <- function(n,CL,no.z){
  #n: smaple size
  #CL: Censoring level
  #no.z: number of responses (3 or 5)
#-----------------------------------------------------------------------------
  if (no.z==3){
    x<-matrix(c(runif(n),runif(n),runif(n)),n,3)  
    beta1<-c(1,-0.5,2)
    beta2 <- c(5,3,-2)
    beta3 <- c(-1,3,1)
    allbeta <- cbind(beta1,beta2,beta3)
  }
  if (no.z==5){
    x<-matrix(c(runif(n),runif(n),runif(n)),n,3)
    beta1<-c(1,-0.5,2)
    beta2 <- c(5,3,-2)
    beta3 <- c(-1,3,1)
    beta4 <- c(-2,0.1,1)
    beta5 <- c(3,-1,0.5)
    allbeta <- cbind(beta1,beta2,beta3,beta4,beta5)
  }
#-----------------------------------------------------------------------------  
  t1<-0
  for (j in 1:n){
    t1[j]<-2.5*(j-0.5)/n
  }
  t2<-seq(-2,2,length.out = n)
#------------------------------------------------------------------------------
  g1   <- -t1*(sin(-t1^2))
  g2   <-sin(2*t2)+2*exp(-16*t2^2)
#------------------------------------------------------------------------------ 
  error<-rnorm(n,sd=0.5)
  if (no.z==3){
    z1<-x%*%beta1+g1+g2+error
    z2<-x%*%beta2+g1+g2+error
    z3<-x%*%beta3+g1+g2+error
    Z <- matrix(c(z1,z2,z3),n,3)  
  }
  if (no.z==5){
    z1<-x%*%beta1+g1+g2+error
    z2<-x%*%beta2+g1+g2+error
    z3<-x%*%beta3+g1+g2+error
    z4<-x%*%beta4+g1+g2+error
    z5<-x%*%beta5+g1+g2+error
    Z <- matrix(c(z1,z2,z3,z4,z5),n,5)  
  }
  nc  <-matrix(c(t1,t2),n,2)
  allf<-matrix(c(g1,g2),n,2)
#------------------------------------------------------------------------------  
#Generation of the right-censored data-----------------------------------------
censoring <- function(n,CL,z){
  delta <-1-rbern(n,CL)     
  c <- matrix(0,n,1)
  y <- 0
  for (i in 1:n){
    if (delta[i]==0){
      while (Z[i]<=c[i]){
        c[i]<-rnorm(1,mean(z),sd=sd(z))
      }
    }
    else{
        c[i]<-z[i]
    }
  }                          
  for (j in 1:n){
    if (z[j]<=c[j]){
      y[j]<-z[j]
    }
    else{
      y[j]<-c[j]
    }
  }
  return(y)
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
delta <-1-rbern(n,CL) 
if (no.z==3){
  y1     <- censoring(n,CL,z1)
  syn.y1 <- syndata(y1,delta) 
  y2     <- censoring(n,CL,z2)
  syn.y2 <- syndata(y2,delta)
  y3     <- censoring(n,CL,z3)
  syn.y3 <- syndata(y3,delta)
  syn.Y  <- matrix(c(syn.y1,syn.y2,syn.y3),n,3) 
  Y      <- matrix(c(y1,y2,y3),n,3)
  }
if (no.z==5){
    y1     <- censoring(n,CL,z1)
    syn.y1 <- syndata(y1,delta)
    y2     <- censoring(n,CL,z2)
    syn.y2 <- syndata(y2,delta)
    y3     <- censoring(n,CL,z3)
    syn.y3 <- syndata(y3,delta)
    y4     <- censoring(n,CL,z4)
    syn.y4 <- syndata(y4,delta)
    y5     <- censoring(n,CL,z5)
    syn.y5 <- syndata(y5,delta)
    syn.Y  <- matrix(c(syn.y1,syn.y2,syn.y3,syn.y4,syn.y5),n,5)
    Y      <- matrix(c(y1,y2,y3,y4,y5),n,5)
}            

  dat         <-new.env()
  dat$x       <-x
  dat$z       <-Z
  dat$Y       <- Y
  dat$syn.Y   <-syn.Y 
  dat$nc      <-nc
  dat$allbeta <-allbeta
  dat$allf    <-allf
  return(dat)
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
#SMOOTHING FUNCTION WITH REG.SPLINE---------------------------------------------
MRLP<-function(x,y,nc,allf){
  library(optR)
  library(condSURV)
  #x: matrix of parametric covariates
  #y: response variable(s) all
  #nc: Matrix of nonparametric covariate
#Prepare data-------------------------------------------------------------------
  t1 <- nc[,1]
  t2 <- nc[,2]
  g1 <- allf[,1]
  g2 <- allf[,2]
  dimY  <- dim(y)
  no.y <- dimY[2]
  if (no.y==3){
    y1 <- y[,1]
    y2 <- y[,2]
    y3 <- y[,3]  
  }
  if (no.y==5){
    y1 <- y[,1]
    y2 <- y[,2]
    y3 <- y[,3]
    y4 <- y[,4]
    y5 <- y[,5]
  }
  n <- length(y1)
  r <- 2
  p <- 3
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
  ones   <-matrix(1,n,1) 
  W1     <-matrix(0,n,n)
  W2     <-matrix(0,n,n)
  GCV_g1h1 <- 0
  GCV_g1h2 <- 0
  GCV_g1h3 <- 0
  GCV_g1h4 <- 0
  GCV_g1h5 <- 0
  
  GCV_g2h1 <- 0
  GCV_g2h2 <- 0
  GCV_g2h3 <- 0
  GCV_g2h4 <- 0
  GCV_g2h5 <- 0
  
  LO     <-30 
  index1 <-seq(min(t1)-0.1,max(t1)+0.1,length.out=n)
  index2 <-seq(min(t2)-0.5,max(t2)+0.5,length.out=n)
  T1     <-matrix(c(ones,(t1-index1),(t1-index1)^2),n,3)
  T2     <-matrix(c(ones,(t2-index2),(t2-index2)^2),n,3)
  bw_seq <-seq(0.015,0.9,length.out=LO)
  if(no.y==3){
  for (i in 1:LO){
  for (j in 1:n){
    W1[,j] <- NWW(t1,index1[j],bw=bw_seq[i])
    W2[,j] <- NWW(t2,index2[j],kernel="gaussian",bw=bw_seq[i])
  }
    zeros <- matrix(0,n,1)
    
    S_g1  <- T1%*%solve(t(T1)%*%W1%*%T1)%*%t(T1)%*%W1
    S_g2  <- T2%*%solve(t(T2)%*%W2%*%T2)%*%t(T2)%*%W2
    
      xtil_g1  <- (diag(n)-S_g1)%*%x
      ytil1_g1 <- (diag(n)-S_g1)%*%y1
      ytil2_g1 <- (diag(n)-S_g1)%*%y2
      ytil3_g1 <- (diag(n)-S_g1)%*%y3
      
      xtil_g2 <- (diag(n)-S_g2)%*%x
      ytil1_g2 <- (diag(n)-S_g2)%*%y1
      ytil2_g2 <- (diag(n)-S_g2)%*%y2
      ytil3_g2 <- (diag(n)-S_g2)%*%y3
      
      betahat1_g1 <-solve(t(xtil_g1)%*%xtil_g1)%*%t(xtil_g1)%*%ytil1_g1 
      betahat2_g1 <-solve(t(xtil_g1)%*%xtil_g1)%*%t(xtil_g1)%*%ytil2_g1
      betahat3_g1 <-solve(t(xtil_g1)%*%xtil_g1)%*%t(xtil_g1)%*%ytil3_g1
      
      betahat1_g2 <-solve(t(xtil_g2)%*%xtil_g2)%*%t(xtil_g2)%*%ytil1_g2 
      betahat2_g2 <-solve(t(xtil_g2)%*%xtil_g2)%*%t(xtil_g2)%*%ytil2_g2
      betahat3_g2 <-solve(t(xtil_g2)%*%xtil_g2)%*%t(xtil_g2)%*%ytil3_g2
      
      beta1hat <- rowMeans(betahat1_g1,betahat1_g2)
      beta2hat <- rowMeans(betahat2_g1,betahat2_g2)
      beta3hat <- rowMeans(betahat3_g1,betahat3_g2)
      
      g1hat1 <- W1%*%(y1-x%*%beta1hat)
      g1hat2 <- W1%*%(y2-x%*%beta2hat)
      g1hat3 <- W1%*%(y3-x%*%beta3hat)
      
      g2hat1 <- W2%*%(y1-x%*%beta1hat)
      g2hat2 <- W2%*%(y2-x%*%beta2hat)
      g2hat3 <- W2%*%(y3-x%*%beta3hat)
      
      p0 <- tr(S_g1)
      p1 <- tr(S_g2)
      
      GCV_g1h1[i] <-gcvfunc(g1,g1hat1,p0)
      GCV_g1h2[i] <-gcvfunc(g1,g1hat2,p0)
      GCV_g1h3[i] <-gcvfunc(g1,g1hat3,p0)
      
      GCV_g2h1[i] <-gcvfunc(g2,g2hat1,p1)
      GCV_g2h2[i] <-gcvfunc(g2,g2hat2,p1)
      GCV_g2h3[i] <-gcvfunc(g2,g2hat3,p1)
  }
  ming <- min(scale(GCV_g1h1))
  maxg <- max(scale(GCV_g1h1))
  plot(bw_seq,scale(GCV_g1h1),type="l",ylab="GCV score",xlab="bandwidth",ylim=c(ming,maxg),col=2,lty=1,main="Bandwidth Selection for g(t1)")
  par(new=TRUE)
  plot(bw_seq,scale(GCV_g1h2),type="l",ylab="GCV score",xlab="bandwidth",ylim=c(ming,maxg),col=3,lty=2)
  par(new=TRUE)
  plot(bw_seq,scale(GCV_g1h3),type="l",ylab="GCV score",xlab="bandwidth",ylim=c(ming,maxg),col=4,lty=3)
  grid()
  legend(x = "topright",legend = c("GCV(g1) for y1", "GCV(g1) for y2","GCV(g1) for y3"), lty = c(1, 2,3),col = c(2, 3,4),lwd = 1)
  
  ming <- min(scale(GCV_g2h1))
  maxg <- max(scale(GCV_g2h1))
  plot(bw_seq,scale(GCV_g2h1),type="l",ylab="GCV score",xlab="bandwidth",ylim=c(ming,maxg),col=2,lty=1,main="Bandwidth Selection for g(t2)")
  par(new=TRUE)
  plot(bw_seq,scale(GCV_g2h2),type="l",ylab="GCV score",xlab="bandwidth",ylim=c(ming,maxg),col=3,lty=2)
  par(new=TRUE)
  plot(bw_seq,scale(GCV_g2h3),type="l",ylab="GCV score",xlab="bandwidth",ylim=c(ming,maxg),col=4,lty=3)
  grid()
  legend(x = "topright",legend = c("GCV(g2) for y1", "GCV(g2) for y2","GCV(g2) for y3"), lty = c(1, 2,3),col = c(2, 3,4),lwd = 1)
  
  for (i2 in 1:LO){
    if (GCV_g1h1[i2]==min(GCV_g1h1)){
      bwg1h1 <- bw_seq[i2]
    }
    if (GCV_g1h2[i2]==min(GCV_g1h2)){
      bwg1h2 <- bw_seq[i2]
    }
    if (GCV_g1h3[i2]==min(GCV_g1h3)){
      bwg1h3 <- bw_seq[i2]
    }
    
    if (GCV_g2h1[i2]==min(GCV_g2h1)){
      bwg2h1 <- bw_seq[i2]
    }
    if (GCV_g2h2[i2]==min(GCV_g2h2)){
      bwg2h2 <- bw_seq[i2]
    }
    if (GCV_g2h3[i2]==min(GCV_g2h3)){
      bwg2h3 <- bw_seq[i2]
    }
  }
  }
#-------------------------------------------------------------------------------
    if(no.y==5){
      for (i in 1:LO){
      for (j in 1:n){
          W1[,j] <- NWW(t1,index1[j],bw=bw_seq[i])
          W2[,j] <- NWW(t2,index2[j],kernel="gaussian",bw=bw_seq[i])
        }
        zeros <- matrix(0,n,1)
        
      S_g1  <- T1%*%solve(t(T1)%*%W1%*%T1)%*%t(T1)%*%W1
      S_g2  <- T2%*%solve(t(T2)%*%W2%*%T2)%*%t(T2)%*%W2
      
      xtil_g1  <- (diag(n)-S_g1)%*%x
      ytil1_g1 <- (diag(n)-S_g1)%*%y1
      ytil2_g1 <- (diag(n)-S_g1)%*%y2
      ytil3_g1 <- (diag(n)-S_g1)%*%y3
      ytil4_g1 <- (diag(n)-S_g1)%*%y4
      ytil5_g1 <- (diag(n)-S_g1)%*%y5
      
      xtil_g2 <- (diag(n)-S_g2)%*%x
      ytil1_g2 <- (diag(n)-S_g2)%*%y1
      ytil2_g2 <- (diag(n)-S_g2)%*%y2
      ytil3_g2 <- (diag(n)-S_g2)%*%y3
      ytil4_g2 <- (diag(n)-S_g2)%*%y4
      ytil5_g2 <- (diag(n)-S_g2)%*%y5
      
      betahat1_g1 <-solve(t(xtil_g1)%*%xtil_g1)%*%t(xtil_g1)%*%ytil1_g1 
      betahat2_g1 <-solve(t(xtil_g1)%*%xtil_g1)%*%t(xtil_g1)%*%ytil2_g1
      betahat3_g1 <-solve(t(xtil_g1)%*%xtil_g1)%*%t(xtil_g1)%*%ytil3_g1
      betahat4_g1 <-solve(t(xtil_g1)%*%xtil_g1)%*%t(xtil_g1)%*%ytil4_g1
      betahat5_g1 <-solve(t(xtil_g1)%*%xtil_g1)%*%t(xtil_g1)%*%ytil5_g1
      
      betahat1_g2 <-solve(t(xtil_g2)%*%xtil_g2)%*%t(xtil_g2)%*%ytil1_g2 
      betahat2_g2 <-solve(t(xtil_g2)%*%xtil_g2)%*%t(xtil_g2)%*%ytil2_g2
      betahat3_g2 <-solve(t(xtil_g2)%*%xtil_g2)%*%t(xtil_g2)%*%ytil3_g2
      betahat4_g2 <-solve(t(xtil_g2)%*%xtil_g2)%*%t(xtil_g2)%*%ytil4_g2
      betahat5_g2 <-solve(t(xtil_g2)%*%xtil_g2)%*%t(xtil_g2)%*%ytil5_g2
      
      beta1hat <- rowMeans(betahat1_g1,betahat1_g2)
      beta2hat <- rowMeans(betahat2_g1,betahat2_g2)
      beta3hat <- rowMeans(betahat3_g1,betahat3_g2)
      beta4hat <- rowMeans(betahat4_g1,betahat4_g2)
      beta5hat <- rowMeans(betahat5_g1,betahat5_g2)
      
      g1hat1 <- W1%*%(y1-x%*%beta1hat)
      g1hat2 <- W1%*%(y2-x%*%beta2hat)
      g1hat3 <- W1%*%(y3-x%*%beta3hat)
      g1hat4 <- W1%*%(y4-x%*%beta4hat)
      g1hat5 <- W1%*%(y5-x%*%beta5hat)
      
      g2hat1 <- W2%*%(y1-x%*%beta1hat)
      g2hat2 <- W2%*%(y2-x%*%beta2hat)
      g2hat3 <- W2%*%(y3-x%*%beta3hat)
      g2hat4 <- W2%*%(y4-x%*%beta4hat)
      g2hat5 <- W2%*%(y5-x%*%beta5hat)
      
      p0 <- tr(S_g1)
      p1 <- tr(S_g2)
      
      GCV_g1h1[i] <-gcvfunc(g1,g1hat1,p0)
      GCV_g1h2[i] <-gcvfunc(g1,g1hat2,p0)
      GCV_g1h3[i] <-gcvfunc(g1,g1hat3,p0)
      GCV_g1h4[i] <-gcvfunc(g1,g1hat4,p0)
      GCV_g1h5[i] <-gcvfunc(g1,g1hat5,p0)
      
      GCV_g2h1[i] <-gcvfunc(g2,g2hat1,p1)
      GCV_g2h2[i] <-gcvfunc(g2,g2hat2,p1)
      GCV_g2h3[i] <-gcvfunc(g2,g2hat3,p1)
      GCV_g2h4[i] <-gcvfunc(g2,g2hat4,p1)
      GCV_g2h5[i] <-gcvfunc(g2,g2hat5,p1)
      }
      
      ming <- min(scale(GCV_g1h1))
      maxg <- max(scale(GCV_g1h1))
      plot(bw_seq,scale(GCV_g1h1),type="l",ylab="GCV score",xlab="bandwidth",ylim=c(ming,maxg),col=2,lty=1,main="Bandwidth Selection for g(t1)")
      par(new=TRUE)
      plot(bw_seq,scale(GCV_g1h2),type="l",ylab="GCV score",xlab="bandwidth",ylim=c(ming,maxg),col=3,lty=2)
      par(new=TRUE)
      plot(bw_seq,scale(GCV_g1h3),type="l",ylab="GCV score",xlab="bandwidth",ylim=c(ming,maxg),col=4,lty=3)
      par(new=TRUE)
      plot(bw_seq,scale(GCV_g1h4),type="l",ylab="GCV score",xlab="bandwidth",ylim=c(ming,maxg),col=5,lty=4)
      par(new=TRUE)
      plot(bw_seq,scale(GCV_g1h5),type="l",ylab="GCV score",xlab="bandwidth",ylim=c(ming,maxg),col=6,lty=5)
      grid()
      legend(x = "topright",legend = c("GCV(g1) for y1", "GCV(g1) for y2","GCV(g1) for y3","GCV(g1) for y4","GCV(g1) for y5"), lty = c(1, 2,3,4,5),col = c(2, 3,4,5,6),lwd = 1)
      
      ming <- min(scale(GCV_g2h1))
      maxg <- max(scale(GCV_g2h1))
      plot(bw_seq,scale(GCV_g2h1),type="l",ylab="GCV score",xlab="bandwidth",ylim=c(ming,maxg),col=2,lty=1,main="Bandwidth Selection for g(t2)")
      par(new=TRUE)
      plot(bw_seq,scale(GCV_g2h2),type="l",ylab="GCV score",xlab="bandwidth",ylim=c(ming,maxg),col=3,lty=2)
      par(new=TRUE)
      plot(bw_seq,scale(GCV_g2h3),type="l",ylab="GCV score",xlab="bandwidth",ylim=c(ming,maxg),col=4,lty=3)
      par(new=TRUE)
      plot(bw_seq,scale(GCV_g2h4),type="l",ylab="GCV score",xlab="bandwidth",ylim=c(ming,maxg),col=5,lty=4)
      par(new=TRUE)
      plot(bw_seq,scale(GCV_g2h5),type="l",ylab="GCV score",xlab="bandwidth",ylim=c(ming,maxg),col=6,lty=5)
      grid()
      legend(x = "topright",legend = c("GCV(g2) for y1", "GCV(g2) for y2","GCV(g2) for y3","GCV(g2) for y4","GCV(g2) for y5"), lty = c(1, 2,3,4,5),col = c(2, 3,4,5,6),lwd = 1)
      
    for (i2 in 1:LO){
      if (GCV_g1h1[i2]==min(GCV_g1h1)){
        bwg1h1 <- bw_seq[i2]
      }
      if (GCV_g1h2[i2]==min(GCV_g1h2)){
        bwg1h2 <- bw_seq[i2]
      }
      if (GCV_g1h3[i2]==min(GCV_g1h3)){
        bwg1h3 <- bw_seq[i2]
      }
      if (GCV_g1h4[i2]==min(GCV_g1h4)){
        bwg1h4 <- bw_seq[i2]
      }
      if (GCV_g1h5[i2]==min(GCV_g1h5)){
        bwg1h5 <- bw_seq[i2]
      }
      
      
      if (GCV_g2h1[i2]==min(GCV_g2h1)){
        bwg2h1 <- bw_seq[i2]
      }
      if (GCV_g2h2[i2]==min(GCV_g2h2)){
        bwg2h2 <- bw_seq[i2]
      }
      if (GCV_g2h3[i2]==min(GCV_g2h3)){
        bwg2h3 <- bw_seq[i2]
      }
      if (GCV_g2h4[i2]==min(GCV_g2h4)){
        bwg2h4 <- bw_seq[i2]
      }
      if (GCV_g2h5[i2]==min(GCV_g2h5)){
        bwg2h5 <- bw_seq[i2]
      }
    }
  }
#ESTIMATION PART---------------------------------------------------------------
    if(no.y==3){
      ones   <-matrix(1,n,1) 
      W1h1     <-matrix(0,n,n)
      W1h2     <-matrix(0,n,n)
      W1h3     <-matrix(0,n,n)
      
      W2h1     <-matrix(0,n,n)
      W2h2     <-matrix(0,n,n)
      W2h3     <-matrix(0,n,n)
        for (j in 1:n){
          W1h1[,j] <- NWW(t1,index1[j],bw=bwg1h1)
          W1h2[,j] <- NWW(t1,index1[j],bw=bwg1h2)
          W1h3[,j] <- NWW(t1,index1[j],bw=bwg1h3)
          
          W2h1[,j] <- NWW(t2,index2[j],bw=bwg2h1)
          W2h2[,j] <- NWW(t2,index2[j],bw=bwg2h2)
          W2h3[,j] <- NWW(t2,index2[j],bw=bwg2h3)
        }
  
  S_g1h1  <- T1%*%solve(t(T1)%*%W1h1%*%T1)%*%t(T1)%*%W1h1
  S_g1h2  <- T1%*%solve(t(T1)%*%W1h2%*%T1)%*%t(T1)%*%W1h2
  S_g1h3  <- T1%*%solve(t(T1)%*%W1h3%*%T1)%*%t(T1)%*%W1h3
  
  S_g2h1  <- T2%*%solve(t(T2)%*%W2h1%*%T2)%*%t(T2)%*%W2h1
  S_g2h2  <- T2%*%solve(t(T2)%*%W2h2%*%T2)%*%t(T2)%*%W2h2
  S_g2h3  <- T2%*%solve(t(T2)%*%W2h3%*%T2)%*%t(T2)%*%W2h3
#------------------------------------------------------------------------------
  xtil1_g1  <- (diag(n)-S_g1h1)%*%x
  xtil2_g1  <- (diag(n)-S_g1h2)%*%x
  xtil3_g1  <- (diag(n)-S_g1h3)%*%x
  
  ytil1_g1 <- (diag(n)-S_g1h1)%*%y1
  ytil2_g1 <- (diag(n)-S_g1h2)%*%y2
  ytil3_g1 <- (diag(n)-S_g1h3)%*%y3
  
  xtil1_g2  <- (diag(n)-S_g2h1)%*%x
  xtil2_g2  <- (diag(n)-S_g2h2)%*%x
  xtil3_g2  <- (diag(n)-S_g2h3)%*%x
  
  ytil1_g2 <- (diag(n)-S_g2h1)%*%y1
  ytil2_g2 <- (diag(n)-S_g2h2)%*%y2
  ytil3_g2 <- (diag(n)-S_g2h3)%*%y3
  
  betahat1_g1 <-solve(t(xtil1_g1)%*%xtil1_g1)%*%t(xtil1_g1)%*%ytil1_g1 
  betahat2_g1 <-solve(t(xtil2_g1)%*%xtil2_g1)%*%t(xtil2_g1)%*%ytil2_g1
  betahat3_g1 <-solve(t(xtil3_g1)%*%xtil3_g1)%*%t(xtil3_g1)%*%ytil3_g1
  
  betahat1_g2 <-solve(t(xtil1_g2)%*%xtil1_g2)%*%t(xtil1_g2)%*%ytil1_g2 
  betahat2_g2 <-solve(t(xtil2_g2)%*%xtil2_g2)%*%t(xtil2_g2)%*%ytil2_g2
  betahat3_g2 <-solve(t(xtil3_g2)%*%xtil3_g2)%*%t(xtil3_g2)%*%ytil3_g2
  
  beta1hat <- rowMeans(betahat1_g1,betahat1_g2)
  beta2hat <- rowMeans(betahat2_g1,betahat2_g2)
  beta3hat <- rowMeans(betahat3_g1,betahat3_g2)
  
  g1hat1 <- W1h1%*%(y1-x%*%beta1hat)
  g1hat2 <- W1h2%*%(y2-x%*%beta2hat)
  g1hat3 <- W1h3%*%(y3-x%*%beta3hat)
  
  g2hat1 <- W2h1%*%(y1-x%*%beta1hat)
  g2hat2 <- W2h2%*%(y2-x%*%beta2hat)
  g2hat3 <- W2h3%*%(y3-x%*%beta3hat)
  
  y1hat <- x%*%beta1hat+g1hat1+g2hat1
  y2hat <- x%*%beta2hat+g1hat2+g2hat2
  y3hat <- x%*%beta3hat+g1hat3+g2hat3

  
  res <- new.env()
  
  res$beta1hat <- beta1hat
  res$beta2hat <- beta2hat
  res$beta3hat <- beta3hat
  #------------------------------------------------------------------------------    
  res$g1hat1 <- g1hat1
  res$g1hat2 <- g1hat2
  res$g1hat3 <- g1hat3
  
  res$g2hat1 <- g2hat1
  res$g2hat2 <- g2hat2
  res$g2hat3 <- g2hat3
  
  #------------------------------------------------------------------------------
  res$y1hat <- y1hat
  res$y2hat <- y2hat
  res$y3hat <- y3hat
  
  }
  if(no.y==5){
    ones     <-matrix(1,n,1) 
    W1h1     <-matrix(0,n,n)
    W1h2     <-matrix(0,n,n)
    W1h3     <-matrix(0,n,n)
    W1h4     <-matrix(0,n,n)
    W1h5     <-matrix(0,n,n)
    
    W2h1     <-matrix(0,n,n)
    W2h2     <-matrix(0,n,n)
    W2h3     <-matrix(0,n,n)
    W2h4     <-matrix(0,n,n)
    W2h5     <-matrix(0,n,n)
    
    for (j in 1:n){
      W1h1[,j] <- NWW(t1,index1[j],bw=bwg1h1)
      W1h2[,j] <- NWW(t1,index1[j],bw=bwg1h2)
      W1h3[,j] <- NWW(t1,index1[j],bw=bwg1h3)
      W1h4[,j] <- NWW(t1,index1[j],bw=bwg1h4)
      W1h5[,j] <- NWW(t1,index1[j],bw=bwg1h5)
      
      W2h1[,j] <- NWW(t2,index2[j],bw=bwg2h1)
      W2h2[,j] <- NWW(t2,index2[j],bw=bwg2h2)
      W2h3[,j] <- NWW(t2,index2[j],bw=bwg2h3)
      W2h4[,j] <- NWW(t2,index2[j],bw=bwg2h4)
      W2h5[,j] <- NWW(t2,index2[j],bw=bwg2h5)
    }
    
    S_g1h1  <- T1%*%solve(t(T1)%*%W1h1%*%T1)%*%t(T1)%*%W1h1
    S_g1h2  <- T1%*%solve(t(T1)%*%W1h2%*%T1)%*%t(T1)%*%W1h2
    S_g1h3  <- T1%*%solve(t(T1)%*%W1h3%*%T1)%*%t(T1)%*%W1h3
    S_g1h4  <- T1%*%solve(t(T1)%*%W1h4%*%T1)%*%t(T1)%*%W1h4
    S_g1h5  <- T1%*%solve(t(T1)%*%W1h5%*%T1)%*%t(T1)%*%W1h5
    
    S_g2h1  <- T2%*%solve(t(T2)%*%W2h1%*%T2)%*%t(T2)%*%W2h1
    S_g2h2  <- T2%*%solve(t(T2)%*%W2h2%*%T2)%*%t(T2)%*%W2h2
    S_g2h3  <- T2%*%solve(t(T2)%*%W2h3%*%T2)%*%t(T2)%*%W2h3
    S_g2h4  <- T2%*%solve(t(T2)%*%W2h4%*%T2)%*%t(T2)%*%W2h4
    S_g2h5  <- T2%*%solve(t(T2)%*%W2h5%*%T2)%*%t(T2)%*%W2h5
    #------------------------------------------------------------------------------
    xtil1_g1  <- (diag(n)-S_g1h1)%*%x
    xtil2_g1  <- (diag(n)-S_g1h2)%*%x
    xtil3_g1  <- (diag(n)-S_g1h3)%*%x
    xtil4_g1  <- (diag(n)-S_g1h4)%*%x
    xtil5_g1  <- (diag(n)-S_g1h5)%*%x
    
    ytil1_g1 <- (diag(n)-S_g1h1)%*%y1
    ytil2_g1 <- (diag(n)-S_g1h2)%*%y2
    ytil3_g1 <- (diag(n)-S_g1h3)%*%y3
    ytil4_g1 <- (diag(n)-S_g1h4)%*%y4
    ytil5_g1 <- (diag(n)-S_g1h5)%*%y5
    
    xtil1_g2  <- (diag(n)-S_g2h1)%*%x
    xtil2_g2  <- (diag(n)-S_g2h2)%*%x
    xtil3_g2  <- (diag(n)-S_g2h3)%*%x
    xtil4_g2  <- (diag(n)-S_g2h4)%*%x
    xtil5_g2  <- (diag(n)-S_g2h5)%*%x
    
    ytil1_g2 <- (diag(n)-S_g2h1)%*%y1
    ytil2_g2 <- (diag(n)-S_g2h2)%*%y2
    ytil3_g2 <- (diag(n)-S_g2h3)%*%y3
    ytil4_g2 <- (diag(n)-S_g2h4)%*%y4
    ytil5_g2 <- (diag(n)-S_g2h5)%*%y5
    
    betahat1_g1 <-solve(t(xtil1_g1)%*%xtil1_g1)%*%t(xtil1_g1)%*%ytil1_g1 
    betahat2_g1 <-solve(t(xtil2_g1)%*%xtil2_g1)%*%t(xtil2_g1)%*%ytil2_g1
    betahat3_g1 <-solve(t(xtil3_g1)%*%xtil3_g1)%*%t(xtil3_g1)%*%ytil3_g1
    betahat4_g1 <-solve(t(xtil4_g1)%*%xtil4_g1)%*%t(xtil4_g1)%*%ytil4_g1
    betahat5_g1 <-solve(t(xtil5_g1)%*%xtil5_g1)%*%t(xtil5_g1)%*%ytil5_g1
    
    betahat1_g2 <-solve(t(xtil1_g2)%*%xtil1_g2)%*%t(xtil1_g2)%*%ytil1_g2 
    betahat2_g2 <-solve(t(xtil2_g2)%*%xtil2_g2)%*%t(xtil2_g2)%*%ytil2_g2
    betahat3_g2 <-solve(t(xtil3_g2)%*%xtil3_g2)%*%t(xtil3_g2)%*%ytil3_g2
    betahat4_g2 <-solve(t(xtil4_g2)%*%xtil4_g2)%*%t(xtil4_g2)%*%ytil4_g2
    betahat5_g2 <-solve(t(xtil5_g2)%*%xtil5_g2)%*%t(xtil5_g2)%*%ytil5_g2
    
    beta1hat <- rowMeans(betahat1_g1,betahat1_g2)
    beta2hat <- rowMeans(betahat2_g1,betahat2_g2)
    beta3hat <- rowMeans(betahat3_g1,betahat3_g2)
    beta4hat <- rowMeans(betahat4_g1,betahat4_g2)
    beta5hat <- rowMeans(betahat5_g1,betahat5_g2)
    
    g1hat1 <- W1h1%*%(y1-x%*%beta1hat)
    g1hat2 <- W1h2%*%(y2-x%*%beta2hat)
    g1hat3 <- W1h3%*%(y3-x%*%beta3hat)
    g1hat4 <- W1h4%*%(y4-x%*%beta4hat)
    g1hat5 <- W1h5%*%(y5-x%*%beta5hat)
    
    g2hat1 <- W2h1%*%(y1-x%*%beta1hat)
    g2hat2 <- W2h2%*%(y2-x%*%beta2hat)
    g2hat3 <- W2h3%*%(y3-x%*%beta3hat)
    g2hat4 <- W2h4%*%(y4-x%*%beta4hat)
    g2hat5 <- W2h5%*%(y5-x%*%beta5hat)
    
    y1hat <- x%*%beta1hat+g1hat1+g2hat1
    y2hat <- x%*%beta2hat+g1hat2+g2hat2
    y3hat <- x%*%beta3hat+g1hat3+g2hat3
    y4hat <- x%*%beta4hat+g1hat4+g2hat4
    y5hat <- x%*%beta5hat+g1hat5+g2hat5
    
    
    res <- new.env()
    
    res$beta1hat <- beta1hat
    res$beta2hat <- beta2hat
    res$beta3hat <- beta3hat
    res$beta4hat <- beta4hat
    res$beta5hat <- beta5hat
#------------------------------------------------------------------------------    
    res$g1hat1 <- g1hat1
    res$g1hat2 <- g1hat2
    res$g1hat3 <- g1hat3
    res$g1hat4 <- g1hat4
    res$g1hat5 <- g1hat5
    res$g2hat1 <- g2hat1
    res$g2hat2 <- g2hat2
    res$g2hat3 <- g2hat3
    res$g2hat4 <- g2hat4
    res$g2hat5 <- g2hat5
#------------------------------------------------------------------------------
    res$y1hat <- y1hat
    res$y2hat <- y2hat
    res$y3hat <- y3hat
    res$y4hat <- y4hat
    res$y5hat <- y5hat
  }
#-------------------------------------------------------------------------------
  return(res)
  #-----------------------------------------------------------------------------
  #-----------------------------------------------------------------------------
}
MRLL<-function(x,y,nc,allf){
  library(optR)
  library(condSURV)
  #x: matrix of parametric covariates
  #y: response variable(s) all
  #nc: Matrix of nonparametric covariate
  #Prepare data-------------------------------------------------------------------
  t1 <- nc[,1]
  t2 <- nc[,2]
  g1 <- allf[,1]
  g2 <- allf[,2]
  dimY  <- dim(y)
  no.y <- dimY[2]
  if (no.y==3){
    y1 <- y[,1]
    y2 <- y[,2]
    y3 <- y[,3]  
  }
  if (no.y==5){
    y1 <- y[,1]
    y2 <- y[,2]
    y3 <- y[,3]
    y4 <- y[,4]
    y5 <- y[,5]
  }
  n <- length(y1)
  r <- 2
  p <- 3
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
  ones   <-matrix(1,n,1) 
  W1     <-matrix(0,n,n)
  W2     <-matrix(0,n,n)
  GCV_g1h1 <- 0
  GCV_g1h2 <- 0
  GCV_g1h3 <- 0
  GCV_g1h4 <- 0
  GCV_g1h5 <- 0
  
  GCV_g2h1 <- 0
  GCV_g2h2 <- 0
  GCV_g2h3 <- 0
  GCV_g2h4 <- 0
  GCV_g2h5 <- 0
  
  LO     <-30 
  index1 <-seq(min(t1)-0.1,max(t1)+0.1,length.out=n)
  index2 <-seq(min(t2)-0.5,max(t2)+0.5,length.out=n)
  T1     <-matrix(c(ones,(t1-index1)),n,2)
  T2     <-matrix(c(ones,(t2-index2)),n,2)
  bw_seq <-seq(0.015,0.9,length.out=LO)
  if(no.y==3){
    for (i in 1:LO){
      for (j in 1:n){
        W1[,j] <- NWW(t1,index1[j],bw=bw_seq[i])
        W2[,j] <- NWW(t2,index2[j],kernel="gaussian",bw=bw_seq[i])
      }
      zeros <- matrix(0,n,1)
      
      S_g1  <- T1%*%solve(t(T1)%*%W1%*%T1)%*%t(T1)%*%W1
      S_g2  <- T2%*%solve(t(T2)%*%W2%*%T2)%*%t(T2)%*%W2
      
      xtil_g1  <- (diag(n)-S_g1)%*%x
      ytil1_g1 <- (diag(n)-S_g1)%*%y1
      ytil2_g1 <- (diag(n)-S_g1)%*%y2
      ytil3_g1 <- (diag(n)-S_g1)%*%y3
      
      xtil_g2 <- (diag(n)-S_g2)%*%x
      ytil1_g2 <- (diag(n)-S_g2)%*%y1
      ytil2_g2 <- (diag(n)-S_g2)%*%y2
      ytil3_g2 <- (diag(n)-S_g2)%*%y3
      
      betahat1_g1 <-solve(t(xtil_g1)%*%xtil_g1)%*%t(xtil_g1)%*%ytil1_g1 
      betahat2_g1 <-solve(t(xtil_g1)%*%xtil_g1)%*%t(xtil_g1)%*%ytil2_g1
      betahat3_g1 <-solve(t(xtil_g1)%*%xtil_g1)%*%t(xtil_g1)%*%ytil3_g1
      
      betahat1_g2 <-solve(t(xtil_g2)%*%xtil_g2)%*%t(xtil_g2)%*%ytil1_g2 
      betahat2_g2 <-solve(t(xtil_g2)%*%xtil_g2)%*%t(xtil_g2)%*%ytil2_g2
      betahat3_g2 <-solve(t(xtil_g2)%*%xtil_g2)%*%t(xtil_g2)%*%ytil3_g2
      
      beta1hat <- rowMeans(betahat1_g1,betahat1_g2)
      beta2hat <- rowMeans(betahat2_g1,betahat2_g2)
      beta3hat <- rowMeans(betahat3_g1,betahat3_g2)
      
      g1hat1 <- W1%*%(y1-x%*%beta1hat)
      g1hat2 <- W1%*%(y2-x%*%beta2hat)
      g1hat3 <- W1%*%(y3-x%*%beta3hat)
      
      g2hat1 <- W2%*%(y1-x%*%beta1hat)
      g2hat2 <- W2%*%(y2-x%*%beta2hat)
      g2hat3 <- W2%*%(y3-x%*%beta3hat)
      
      p0 <- tr(S_g1)
      p1 <- tr(S_g2)
      
      GCV_g1h1[i] <-gcvfunc(g1,g1hat1,p0)
      GCV_g1h2[i] <-gcvfunc(g1,g1hat2,p0)
      GCV_g1h3[i] <-gcvfunc(g1,g1hat3,p0)
      
      GCV_g2h1[i] <-gcvfunc(g2,g2hat1,p1)
      GCV_g2h2[i] <-gcvfunc(g2,g2hat2,p1)
      GCV_g2h3[i] <-gcvfunc(g2,g2hat3,p1)
    }
    ming <- min(scale(GCV_g1h1))
    maxg <- max(scale(GCV_g1h1))
    plot(bw_seq,scale(GCV_g1h1),type="l",ylab="GCV score",xlab="bandwidth",ylim=c(ming,maxg),col=2,lty=1,main="Bandwidth Selection for g(t1)")
    par(new=TRUE)
    plot(bw_seq,scale(GCV_g1h2),type="l",ylab="GCV score",xlab="bandwidth",ylim=c(ming,maxg),col=3,lty=2)
    par(new=TRUE)
    plot(bw_seq,scale(GCV_g1h3),type="l",ylab="GCV score",xlab="bandwidth",ylim=c(ming,maxg),col=4,lty=3)
    grid()
    legend(x = "topright",legend = c("GCV(g1) for y1", "GCV(g1) for y2","GCV(g1) for y3"), lty = c(1, 2,3),col = c(2, 3,4),lwd = 1)
    
    ming <- min(scale(GCV_g2h1))
    maxg <- max(scale(GCV_g2h1))
    plot(bw_seq,scale(GCV_g2h1),type="l",ylab="GCV score",xlab="bandwidth",ylim=c(ming,maxg),col=2,lty=1,main="Bandwidth Selection for g(t2)")
    par(new=TRUE)
    plot(bw_seq,scale(GCV_g2h2),type="l",ylab="GCV score",xlab="bandwidth",ylim=c(ming,maxg),col=3,lty=2)
    par(new=TRUE)
    plot(bw_seq,scale(GCV_g2h3),type="l",ylab="GCV score",xlab="bandwidth",ylim=c(ming,maxg),col=4,lty=3)
    grid()
    legend(x = "topright",legend = c("GCV(g2) for y1", "GCV(g2) for y2","GCV(g2) for y3"), lty = c(1, 2,3),col = c(2, 3,4),lwd = 1)
    
    for (i2 in 1:LO){
      if (GCV_g1h1[i2]==min(GCV_g1h1)){
        bwg1h1 <- bw_seq[i2]
      }
      if (GCV_g1h2[i2]==min(GCV_g1h2)){
        bwg1h2 <- bw_seq[i2]
      }
      if (GCV_g1h3[i2]==min(GCV_g1h3)){
        bwg1h3 <- bw_seq[i2]
      }
      
      if (GCV_g2h1[i2]==min(GCV_g2h1)){
        bwg2h1 <- bw_seq[i2]
      }
      if (GCV_g2h2[i2]==min(GCV_g2h2)){
        bwg2h2 <- bw_seq[i2]
      }
      if (GCV_g2h3[i2]==min(GCV_g2h3)){
        bwg2h3 <- bw_seq[i2]
      }
    }
  }
  #-------------------------------------------------------------------------------
  if(no.y==5){
    for (i in 1:LO){
      for (j in 1:n){
        W1[,j] <- NWW(t1,index1[j],bw=bw_seq[i])
        W2[,j] <- NWW(t2,index2[j],kernel="gaussian",bw=bw_seq[i])
      }
      zeros <- matrix(0,n,1)
      
      S_g1  <- T1%*%solve(t(T1)%*%W1%*%T1)%*%t(T1)%*%W1
      S_g2  <- T2%*%solve(t(T2)%*%W2%*%T2)%*%t(T2)%*%W2
      
      xtil_g1  <- (diag(n)-S_g1)%*%x
      ytil1_g1 <- (diag(n)-S_g1)%*%y1
      ytil2_g1 <- (diag(n)-S_g1)%*%y2
      ytil3_g1 <- (diag(n)-S_g1)%*%y3
      ytil4_g1 <- (diag(n)-S_g1)%*%y4
      ytil5_g1 <- (diag(n)-S_g1)%*%y5
      
      xtil_g2 <- (diag(n)-S_g2)%*%x
      ytil1_g2 <- (diag(n)-S_g2)%*%y1
      ytil2_g2 <- (diag(n)-S_g2)%*%y2
      ytil3_g2 <- (diag(n)-S_g2)%*%y3
      ytil4_g2 <- (diag(n)-S_g2)%*%y4
      ytil5_g2 <- (diag(n)-S_g2)%*%y5
      
      betahat1_g1 <-solve(t(xtil_g1)%*%xtil_g1)%*%t(xtil_g1)%*%ytil1_g1 
      betahat2_g1 <-solve(t(xtil_g1)%*%xtil_g1)%*%t(xtil_g1)%*%ytil2_g1
      betahat3_g1 <-solve(t(xtil_g1)%*%xtil_g1)%*%t(xtil_g1)%*%ytil3_g1
      betahat4_g1 <-solve(t(xtil_g1)%*%xtil_g1)%*%t(xtil_g1)%*%ytil4_g1
      betahat5_g1 <-solve(t(xtil_g1)%*%xtil_g1)%*%t(xtil_g1)%*%ytil5_g1
      
      betahat1_g2 <-solve(t(xtil_g2)%*%xtil_g2)%*%t(xtil_g2)%*%ytil1_g2 
      betahat2_g2 <-solve(t(xtil_g2)%*%xtil_g2)%*%t(xtil_g2)%*%ytil2_g2
      betahat3_g2 <-solve(t(xtil_g2)%*%xtil_g2)%*%t(xtil_g2)%*%ytil3_g2
      betahat4_g2 <-solve(t(xtil_g2)%*%xtil_g2)%*%t(xtil_g2)%*%ytil4_g2
      betahat5_g2 <-solve(t(xtil_g2)%*%xtil_g2)%*%t(xtil_g2)%*%ytil5_g2
      
      beta1hat <- rowMeans(betahat1_g1,betahat1_g2)
      beta2hat <- rowMeans(betahat2_g1,betahat2_g2)
      beta3hat <- rowMeans(betahat3_g1,betahat3_g2)
      beta4hat <- rowMeans(betahat4_g1,betahat4_g2)
      beta5hat <- rowMeans(betahat5_g1,betahat5_g2)
      
      g1hat1 <- W1%*%(y1-x%*%beta1hat)
      g1hat2 <- W1%*%(y2-x%*%beta2hat)
      g1hat3 <- W1%*%(y3-x%*%beta3hat)
      g1hat4 <- W1%*%(y4-x%*%beta4hat)
      g1hat5 <- W1%*%(y5-x%*%beta5hat)
      
      g2hat1 <- W2%*%(y1-x%*%beta1hat)
      g2hat2 <- W2%*%(y2-x%*%beta2hat)
      g2hat3 <- W2%*%(y3-x%*%beta3hat)
      g2hat4 <- W2%*%(y4-x%*%beta4hat)
      g2hat5 <- W2%*%(y5-x%*%beta5hat)
      
      p0 <- tr(S_g1)
      p1 <- tr(S_g2)
      
      GCV_g1h1[i] <-gcvfunc(g1,g1hat1,p0)
      GCV_g1h2[i] <-gcvfunc(g1,g1hat2,p0)
      GCV_g1h3[i] <-gcvfunc(g1,g1hat3,p0)
      GCV_g1h4[i] <-gcvfunc(g1,g1hat4,p0)
      GCV_g1h5[i] <-gcvfunc(g1,g1hat5,p0)
      
      GCV_g2h1[i] <-gcvfunc(g2,g2hat1,p1)
      GCV_g2h2[i] <-gcvfunc(g2,g2hat2,p1)
      GCV_g2h3[i] <-gcvfunc(g2,g2hat3,p1)
      GCV_g2h4[i] <-gcvfunc(g2,g2hat4,p1)
      GCV_g2h5[i] <-gcvfunc(g2,g2hat5,p1)
    }
    
    ming <- min(scale(GCV_g1h1))
    maxg <- max(scale(GCV_g1h1))
    plot(bw_seq,scale(GCV_g1h1),type="l",ylab="GCV score",xlab="bandwidth",ylim=c(ming,maxg),col=2,lty=1,main="Bandwidth Selection for g(t1)")
    par(new=TRUE)
    plot(bw_seq,scale(GCV_g1h2),type="l",ylab="GCV score",xlab="bandwidth",ylim=c(ming,maxg),col=3,lty=2)
    par(new=TRUE)
    plot(bw_seq,scale(GCV_g1h3),type="l",ylab="GCV score",xlab="bandwidth",ylim=c(ming,maxg),col=4,lty=3)
    par(new=TRUE)
    plot(bw_seq,scale(GCV_g1h4),type="l",ylab="GCV score",xlab="bandwidth",ylim=c(ming,maxg),col=5,lty=4)
    par(new=TRUE)
    plot(bw_seq,scale(GCV_g1h5),type="l",ylab="GCV score",xlab="bandwidth",ylim=c(ming,maxg),col=6,lty=5)
    grid()
    legend(x = "topright",legend = c("GCV(g1) for y1", "GCV(g1) for y2","GCV(g1) for y3","GCV(g1) for y4","GCV(g1) for y5"), lty = c(1, 2,3,4,5),col = c(2, 3,4,5,6),lwd = 1)
    
    ming <- min(scale(GCV_g2h1))
    maxg <- max(scale(GCV_g2h1))
    plot(bw_seq,scale(GCV_g2h1),type="l",ylab="GCV score",xlab="bandwidth",ylim=c(ming,maxg),col=2,lty=1,main="Bandwidth Selection for g(t2)")
    par(new=TRUE)
    plot(bw_seq,scale(GCV_g2h2),type="l",ylab="GCV score",xlab="bandwidth",ylim=c(ming,maxg),col=3,lty=2)
    par(new=TRUE)
    plot(bw_seq,scale(GCV_g2h3),type="l",ylab="GCV score",xlab="bandwidth",ylim=c(ming,maxg),col=4,lty=3)
    par(new=TRUE)
    plot(bw_seq,scale(GCV_g2h4),type="l",ylab="GCV score",xlab="bandwidth",ylim=c(ming,maxg),col=5,lty=4)
    par(new=TRUE)
    plot(bw_seq,scale(GCV_g2h5),type="l",ylab="GCV score",xlab="bandwidth",ylim=c(ming,maxg),col=6,lty=5)
    grid()
    legend(x = "topright",legend = c("GCV(g2) for y1", "GCV(g2) for y2","GCV(g2) for y3","GCV(g2) for y4","GCV(g2) for y5"), lty = c(1, 2,3,4,5),col = c(2, 3,4,5,6),lwd = 1)
    
    for (i2 in 1:LO){
      if (GCV_g1h1[i2]==min(GCV_g1h1)){
        bwg1h1 <- bw_seq[i2]
      }
      if (GCV_g1h2[i2]==min(GCV_g1h2)){
        bwg1h2 <- bw_seq[i2]
      }
      if (GCV_g1h3[i2]==min(GCV_g1h3)){
        bwg1h3 <- bw_seq[i2]
      }
      if (GCV_g1h4[i2]==min(GCV_g1h4)){
        bwg1h4 <- bw_seq[i2]
      }
      if (GCV_g1h5[i2]==min(GCV_g1h5)){
        bwg1h5 <- bw_seq[i2]
      }
      
      
      if (GCV_g2h1[i2]==min(GCV_g2h1)){
        bwg2h1 <- bw_seq[i2]
      }
      if (GCV_g2h2[i2]==min(GCV_g2h2)){
        bwg2h2 <- bw_seq[i2]
      }
      if (GCV_g2h3[i2]==min(GCV_g2h3)){
        bwg2h3 <- bw_seq[i2]
      }
      if (GCV_g2h4[i2]==min(GCV_g2h4)){
        bwg2h4 <- bw_seq[i2]
      }
      if (GCV_g2h5[i2]==min(GCV_g2h5)){
        bwg2h5 <- bw_seq[i2]
      }
    }
  }
  #ESTIMATION PART---------------------------------------------------------------
  if(no.y==3){
    ones   <-matrix(1,n,1) 
    W1h1     <-matrix(0,n,n)
    W1h2     <-matrix(0,n,n)
    W1h3     <-matrix(0,n,n)
    
    W2h1     <-matrix(0,n,n)
    W2h2     <-matrix(0,n,n)
    W2h3     <-matrix(0,n,n)
    for (j in 1:n){
      W1h1[,j] <- NWW(t1,index1[j],bw=bwg1h1)
      W1h2[,j] <- NWW(t1,index1[j],bw=bwg1h2)
      W1h3[,j] <- NWW(t1,index1[j],bw=bwg1h3)
      
      W2h1[,j] <- NWW(t2,index2[j],bw=bwg2h1)
      W2h2[,j] <- NWW(t2,index2[j],bw=bwg2h2)
      W2h3[,j] <- NWW(t2,index2[j],bw=bwg2h3)
    }
    
    S_g1h1  <- T1%*%solve(t(T1)%*%W1h1%*%T1)%*%t(T1)%*%W1h1
    S_g1h2  <- T1%*%solve(t(T1)%*%W1h2%*%T1)%*%t(T1)%*%W1h2
    S_g1h3  <- T1%*%solve(t(T1)%*%W1h3%*%T1)%*%t(T1)%*%W1h3
    
    S_g2h1  <- T2%*%solve(t(T2)%*%W2h1%*%T2)%*%t(T2)%*%W2h1
    S_g2h2  <- T2%*%solve(t(T2)%*%W2h2%*%T2)%*%t(T2)%*%W2h2
    S_g2h3  <- T2%*%solve(t(T2)%*%W2h3%*%T2)%*%t(T2)%*%W2h3
    #------------------------------------------------------------------------------
    xtil1_g1  <- (diag(n)-S_g1h1)%*%x
    xtil2_g1  <- (diag(n)-S_g1h2)%*%x
    xtil3_g1  <- (diag(n)-S_g1h3)%*%x
    
    ytil1_g1 <- (diag(n)-S_g1h1)%*%y1
    ytil2_g1 <- (diag(n)-S_g1h2)%*%y2
    ytil3_g1 <- (diag(n)-S_g1h3)%*%y3
    
    xtil1_g2  <- (diag(n)-S_g2h1)%*%x
    xtil2_g2  <- (diag(n)-S_g2h2)%*%x
    xtil3_g2  <- (diag(n)-S_g2h3)%*%x
    
    ytil1_g2 <- (diag(n)-S_g2h1)%*%y1
    ytil2_g2 <- (diag(n)-S_g2h2)%*%y2
    ytil3_g2 <- (diag(n)-S_g2h3)%*%y3
    
    betahat1_g1 <-solve(t(xtil1_g1)%*%xtil1_g1)%*%t(xtil1_g1)%*%ytil1_g1 
    betahat2_g1 <-solve(t(xtil2_g1)%*%xtil2_g1)%*%t(xtil2_g1)%*%ytil2_g1
    betahat3_g1 <-solve(t(xtil3_g1)%*%xtil3_g1)%*%t(xtil3_g1)%*%ytil3_g1
    
    betahat1_g2 <-solve(t(xtil1_g2)%*%xtil1_g2)%*%t(xtil1_g2)%*%ytil1_g2 
    betahat2_g2 <-solve(t(xtil2_g2)%*%xtil2_g2)%*%t(xtil2_g2)%*%ytil2_g2
    betahat3_g2 <-solve(t(xtil3_g2)%*%xtil3_g2)%*%t(xtil3_g2)%*%ytil3_g2
    
    beta1hat <- rowMeans(betahat1_g1,betahat1_g2)
    beta2hat <- rowMeans(betahat2_g1,betahat2_g2)
    beta3hat <- rowMeans(betahat3_g1,betahat3_g2)
    
    g1hat1 <- W1h1%*%(y1-x%*%beta1hat)
    g1hat2 <- W1h2%*%(y2-x%*%beta2hat)
    g1hat3 <- W1h3%*%(y3-x%*%beta3hat)
    
    g2hat1 <- W2h1%*%(y1-x%*%beta1hat)
    g2hat2 <- W2h2%*%(y2-x%*%beta2hat)
    g2hat3 <- W2h3%*%(y3-x%*%beta3hat)
    
    y1hat <- x%*%beta1hat+g1hat1+g2hat1
    y2hat <- x%*%beta2hat+g1hat2+g2hat2
    y3hat <- x%*%beta3hat+g1hat3+g2hat3
    
    
    res <- new.env()
    
    res$beta1hat <- beta1hat
    res$beta2hat <- beta2hat
    res$beta3hat <- beta3hat
    #------------------------------------------------------------------------------    
    res$g1hat1 <- g1hat1
    res$g1hat2 <- g1hat2
    res$g1hat3 <- g1hat3
    
    res$g2hat1 <- g2hat1
    res$g2hat2 <- g2hat2
    res$g2hat3 <- g2hat3
    
    #------------------------------------------------------------------------------
    res$y1hat <- y1hat
    res$y2hat <- y2hat
    res$y3hat <- y3hat
    
  }
  if(no.y==5){
    ones     <-matrix(1,n,1) 
    W1h1     <-matrix(0,n,n)
    W1h2     <-matrix(0,n,n)
    W1h3     <-matrix(0,n,n)
    W1h4     <-matrix(0,n,n)
    W1h5     <-matrix(0,n,n)
    
    W2h1     <-matrix(0,n,n)
    W2h2     <-matrix(0,n,n)
    W2h3     <-matrix(0,n,n)
    W2h4     <-matrix(0,n,n)
    W2h5     <-matrix(0,n,n)
    
    for (j in 1:n){
      W1h1[,j] <- NWW(t1,index1[j],bw=bwg1h1)
      W1h2[,j] <- NWW(t1,index1[j],bw=bwg1h2)
      W1h3[,j] <- NWW(t1,index1[j],bw=bwg1h3)
      W1h4[,j] <- NWW(t1,index1[j],bw=bwg1h4)
      W1h5[,j] <- NWW(t1,index1[j],bw=bwg1h5)
      
      W2h1[,j] <- NWW(t2,index2[j],bw=bwg2h1)
      W2h2[,j] <- NWW(t2,index2[j],bw=bwg2h2)
      W2h3[,j] <- NWW(t2,index2[j],bw=bwg2h3)
      W2h4[,j] <- NWW(t2,index2[j],bw=bwg2h4)
      W2h5[,j] <- NWW(t2,index2[j],bw=bwg2h5)
    }
    
    S_g1h1  <- T1%*%solve(t(T1)%*%W1h1%*%T1)%*%t(T1)%*%W1h1
    S_g1h2  <- T1%*%solve(t(T1)%*%W1h2%*%T1)%*%t(T1)%*%W1h2
    S_g1h3  <- T1%*%solve(t(T1)%*%W1h3%*%T1)%*%t(T1)%*%W1h3
    S_g1h4  <- T1%*%solve(t(T1)%*%W1h4%*%T1)%*%t(T1)%*%W1h4
    S_g1h5  <- T1%*%solve(t(T1)%*%W1h5%*%T1)%*%t(T1)%*%W1h5
    
    S_g2h1  <- T2%*%solve(t(T2)%*%W2h1%*%T2)%*%t(T2)%*%W2h1
    S_g2h2  <- T2%*%solve(t(T2)%*%W2h2%*%T2)%*%t(T2)%*%W2h2
    S_g2h3  <- T2%*%solve(t(T2)%*%W2h3%*%T2)%*%t(T2)%*%W2h3
    S_g2h4  <- T2%*%solve(t(T2)%*%W2h4%*%T2)%*%t(T2)%*%W2h4
    S_g2h5  <- T2%*%solve(t(T2)%*%W2h5%*%T2)%*%t(T2)%*%W2h5
    #------------------------------------------------------------------------------
    xtil1_g1  <- (diag(n)-S_g1h1)%*%x
    xtil2_g1  <- (diag(n)-S_g1h2)%*%x
    xtil3_g1  <- (diag(n)-S_g1h3)%*%x
    xtil4_g1  <- (diag(n)-S_g1h4)%*%x
    xtil5_g1  <- (diag(n)-S_g1h5)%*%x
    
    ytil1_g1 <- (diag(n)-S_g1h1)%*%y1
    ytil2_g1 <- (diag(n)-S_g1h2)%*%y2
    ytil3_g1 <- (diag(n)-S_g1h3)%*%y3
    ytil4_g1 <- (diag(n)-S_g1h4)%*%y4
    ytil5_g1 <- (diag(n)-S_g1h5)%*%y5
    
    xtil1_g2  <- (diag(n)-S_g2h1)%*%x
    xtil2_g2  <- (diag(n)-S_g2h2)%*%x
    xtil3_g2  <- (diag(n)-S_g2h3)%*%x
    xtil4_g2  <- (diag(n)-S_g2h4)%*%x
    xtil5_g2  <- (diag(n)-S_g2h5)%*%x
    
    ytil1_g2 <- (diag(n)-S_g2h1)%*%y1
    ytil2_g2 <- (diag(n)-S_g2h2)%*%y2
    ytil3_g2 <- (diag(n)-S_g2h3)%*%y3
    ytil4_g2 <- (diag(n)-S_g2h4)%*%y4
    ytil5_g2 <- (diag(n)-S_g2h5)%*%y5
    
    betahat1_g1 <-solve(t(xtil1_g1)%*%xtil1_g1)%*%t(xtil1_g1)%*%ytil1_g1 
    betahat2_g1 <-solve(t(xtil2_g1)%*%xtil2_g1)%*%t(xtil2_g1)%*%ytil2_g1
    betahat3_g1 <-solve(t(xtil3_g1)%*%xtil3_g1)%*%t(xtil3_g1)%*%ytil3_g1
    betahat4_g1 <-solve(t(xtil4_g1)%*%xtil4_g1)%*%t(xtil4_g1)%*%ytil4_g1
    betahat5_g1 <-solve(t(xtil5_g1)%*%xtil5_g1)%*%t(xtil5_g1)%*%ytil5_g1
    
    betahat1_g2 <-solve(t(xtil1_g2)%*%xtil1_g2)%*%t(xtil1_g2)%*%ytil1_g2 
    betahat2_g2 <-solve(t(xtil2_g2)%*%xtil2_g2)%*%t(xtil2_g2)%*%ytil2_g2
    betahat3_g2 <-solve(t(xtil3_g2)%*%xtil3_g2)%*%t(xtil3_g2)%*%ytil3_g2
    betahat4_g2 <-solve(t(xtil4_g2)%*%xtil4_g2)%*%t(xtil4_g2)%*%ytil4_g2
    betahat5_g2 <-solve(t(xtil5_g2)%*%xtil5_g2)%*%t(xtil5_g2)%*%ytil5_g2
    
    beta1hat <- rowMeans(betahat1_g1,betahat1_g2)
    beta2hat <- rowMeans(betahat2_g1,betahat2_g2)
    beta3hat <- rowMeans(betahat3_g1,betahat3_g2)
    beta4hat <- rowMeans(betahat4_g1,betahat4_g2)
    beta5hat <- rowMeans(betahat5_g1,betahat5_g2)
    
    g1hat1 <- W1h1%*%(y1-x%*%beta1hat)
    g1hat2 <- W1h2%*%(y2-x%*%beta2hat)
    g1hat3 <- W1h3%*%(y3-x%*%beta3hat)
    g1hat4 <- W1h4%*%(y4-x%*%beta4hat)
    g1hat5 <- W1h5%*%(y5-x%*%beta5hat)
    
    g2hat1 <- W2h1%*%(y1-x%*%beta1hat)
    g2hat2 <- W2h2%*%(y2-x%*%beta2hat)
    g2hat3 <- W2h3%*%(y3-x%*%beta3hat)
    g2hat4 <- W2h4%*%(y4-x%*%beta4hat)
    g2hat5 <- W2h5%*%(y5-x%*%beta5hat)
    
    y1hat <- x%*%beta1hat+g1hat1+g2hat1
    y2hat <- x%*%beta2hat+g1hat2+g2hat2
    y3hat <- x%*%beta3hat+g1hat3+g2hat3
    y4hat <- x%*%beta4hat+g1hat4+g2hat4
    y5hat <- x%*%beta5hat+g1hat5+g2hat5
    
    res <- new.env()
    
    res$beta1hat <- beta1hat
    res$beta2hat <- beta2hat
    res$beta3hat <- beta3hat
    res$beta4hat <- beta4hat
    res$beta5hat <- beta5hat
    #------------------------------------------------------------------------------    
    res$g1hat1 <- g1hat1
    res$g1hat2 <- g1hat2
    res$g1hat3 <- g1hat3
    res$g1hat4 <- g1hat4
    res$g1hat5 <- g1hat5
    res$g2hat1 <- g2hat1
    res$g2hat2 <- g2hat2
    res$g2hat3 <- g2hat3
    res$g2hat4 <- g2hat4
    res$g2hat5 <- g2hat5
    #------------------------------------------------------------------------------
    res$y1hat <- y1hat
    res$y2hat <- y2hat
    res$y3hat <- y3hat
    res$y4hat <- y4hat
    res$y5hat <- y5hat
  }
  #-------------------------------------------------------------------------------
  return(res)
  #-----------------------------------------------------------------------------
  #-----------------------------------------------------------------------------
}
MRLC<-function(x,y,nc,allf){
  library(optR)
  library(condSURV)
  #x: matrix of parametric covariates
  #y: response variable(s) all
  #nc: Matrix of nonparametric covariate
  #Prepare data-------------------------------------------------------------------
  t1 <- nc[,1]
  t2 <- nc[,2]
  g1 <- allf[,1]
  g2 <- allf[,2]
  dimY  <- dim(y)
  no.y <- dimY[2]
  if (no.y==3){
    y1 <- y[,1]
    y2 <- y[,2]
    y3 <- y[,3]  
  }
  if (no.y==5){
    y1 <- y[,1]
    y2 <- y[,2]
    y3 <- y[,3]
    y4 <- y[,4]
    y5 <- y[,5]
  }
  n <- length(y1)
  r <- 2
  p <- 3
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
  ones   <-matrix(1,n,1) 
  W1     <-matrix(0,n,n)
  W2     <-matrix(0,n,n)
  GCV_g1h1 <- 0
  GCV_g1h2 <- 0
  GCV_g1h3 <- 0
  GCV_g1h4 <- 0
  GCV_g1h5 <- 0
  
  GCV_g2h1 <- 0
  GCV_g2h2 <- 0
  GCV_g2h3 <- 0
  GCV_g2h4 <- 0
  GCV_g2h5 <- 0
  
  LO     <-30 
  index1 <-seq(min(t1)-0.1,max(t1)+0.1,length.out=n)
  index2 <-seq(min(t2)-0.5,max(t2)+0.5,length.out=n)
  T1     <-matrix(c(ones,(t1-index1),(t1-index1)^2),n,3)
  T2     <-matrix(c(ones,(t2-index2),(t2-index2)^2),n,3)
  bw_seq <-seq(0.015,0.9,length.out=LO)
  if(no.y==3){
    for (i in 1:LO){
      for (j in 1:n){
        W1[,j] <- NWW(t1,index1[j],bw=bw_seq[i])
        W2[,j] <- NWW(t2,index2[j],kernel="gaussian",bw=bw_seq[i])
      }
      zeros <- matrix(0,n,1)
      
      S_g1  <- W1
      S_g2  <- W2
      
      xtil_g1  <- (diag(n)-S_g1)%*%x
      ytil1_g1 <- (diag(n)-S_g1)%*%y1
      ytil2_g1 <- (diag(n)-S_g1)%*%y2
      ytil3_g1 <- (diag(n)-S_g1)%*%y3
      
      xtil_g2 <- (diag(n)-S_g2)%*%x
      ytil1_g2 <- (diag(n)-S_g2)%*%y1
      ytil2_g2 <- (diag(n)-S_g2)%*%y2
      ytil3_g2 <- (diag(n)-S_g2)%*%y3
      
      betahat1_g1 <-solve(t(xtil_g1)%*%xtil_g1)%*%t(xtil_g1)%*%ytil1_g1 
      betahat2_g1 <-solve(t(xtil_g1)%*%xtil_g1)%*%t(xtil_g1)%*%ytil2_g1
      betahat3_g1 <-solve(t(xtil_g1)%*%xtil_g1)%*%t(xtil_g1)%*%ytil3_g1
      
      betahat1_g2 <-solve(t(xtil_g2)%*%xtil_g2)%*%t(xtil_g2)%*%ytil1_g2 
      betahat2_g2 <-solve(t(xtil_g2)%*%xtil_g2)%*%t(xtil_g2)%*%ytil2_g2
      betahat3_g2 <-solve(t(xtil_g2)%*%xtil_g2)%*%t(xtil_g2)%*%ytil3_g2
      
      beta1hat <- rowMeans(betahat1_g1,betahat1_g2)
      beta2hat <- rowMeans(betahat2_g1,betahat2_g2)
      beta3hat <- rowMeans(betahat3_g1,betahat3_g2)
      
      g1hat1 <- W1%*%(y1-x%*%beta1hat)
      g1hat2 <- W1%*%(y2-x%*%beta2hat)
      g1hat3 <- W1%*%(y3-x%*%beta3hat)
      
      g2hat1 <- W2%*%(y1-x%*%beta1hat)
      g2hat2 <- W2%*%(y2-x%*%beta2hat)
      g2hat3 <- W2%*%(y3-x%*%beta3hat)
      
      p0 <- tr(S_g1)
      p1 <- tr(S_g2)
      
      GCV_g1h1[i] <-gcvfunc(g1,g1hat1,p0)
      GCV_g1h2[i] <-gcvfunc(g1,g1hat2,p0)
      GCV_g1h3[i] <-gcvfunc(g1,g1hat3,p0)
      
      GCV_g2h1[i] <-gcvfunc(g2,g2hat1,p1)
      GCV_g2h2[i] <-gcvfunc(g2,g2hat2,p1)
      GCV_g2h3[i] <-gcvfunc(g2,g2hat3,p1)
    }
    ming <- min(scale(GCV_g1h1))
    maxg <- max(scale(GCV_g1h1))
    plot(bw_seq,scale(GCV_g1h1),type="l",ylab="GCV score",xlab="bandwidth",ylim=c(ming,maxg),col=2,lty=1,main="Bandwidth Selection for g(t1)")
    par(new=TRUE)
    plot(bw_seq,scale(GCV_g1h2),type="l",ylab="GCV score",xlab="bandwidth",ylim=c(ming,maxg),col=3,lty=2)
    par(new=TRUE)
    plot(bw_seq,scale(GCV_g1h3),type="l",ylab="GCV score",xlab="bandwidth",ylim=c(ming,maxg),col=4,lty=3)
    grid()
    legend(x = "topright",legend = c("GCV(g1) for y1", "GCV(g1) for y2","GCV(g1) for y3"), lty = c(1, 2,3),col = c(2, 3,4),lwd = 1)
    
    ming <- min(scale(GCV_g2h1))
    maxg <- max(scale(GCV_g2h1))
    plot(bw_seq,scale(GCV_g2h1),type="l",ylab="GCV score",xlab="bandwidth",ylim=c(ming,maxg),col=2,lty=1,main="Bandwidth Selection for g(t2)")
    par(new=TRUE)
    plot(bw_seq,scale(GCV_g2h2),type="l",ylab="GCV score",xlab="bandwidth",ylim=c(ming,maxg),col=3,lty=2)
    par(new=TRUE)
    plot(bw_seq,scale(GCV_g2h3),type="l",ylab="GCV score",xlab="bandwidth",ylim=c(ming,maxg),col=4,lty=3)
    grid()
    legend(x = "topright",legend = c("GCV(g2) for y1", "GCV(g2) for y2","GCV(g2) for y3"), lty = c(1, 2,3),col = c(2, 3,4),lwd = 1)
    
    for (i2 in 1:LO){
      if (GCV_g1h1[i2]==min(GCV_g1h1)){
        bwg1h1 <- bw_seq[i2]
      }
      if (GCV_g1h2[i2]==min(GCV_g1h2)){
        bwg1h2 <- bw_seq[i2]
      }
      if (GCV_g1h3[i2]==min(GCV_g1h3)){
        bwg1h3 <- bw_seq[i2]
      }
      
      if (GCV_g2h1[i2]==min(GCV_g2h1)){
        bwg2h1 <- bw_seq[i2]
      }
      if (GCV_g2h2[i2]==min(GCV_g2h2)){
        bwg2h2 <- bw_seq[i2]
      }
      if (GCV_g2h3[i2]==min(GCV_g2h3)){
        bwg2h3 <- bw_seq[i2]
      }
    }
  }
  #-------------------------------------------------------------------------------
  if(no.y==5){
    for (i in 1:LO){
      for (j in 1:n){
        W1[,j] <- NWW(t1,index1[j],bw=bw_seq[i])
        W2[,j] <- NWW(t2,index2[j],kernel="gaussian",bw=bw_seq[i])
      }
      zeros <- matrix(0,n,1)
      
      S_g1  <- W1
      S_g2  <- W2
      
      xtil_g1  <- (diag(n)-S_g1)%*%x
      ytil1_g1 <- (diag(n)-S_g1)%*%y1
      ytil2_g1 <- (diag(n)-S_g1)%*%y2
      ytil3_g1 <- (diag(n)-S_g1)%*%y3
      ytil4_g1 <- (diag(n)-S_g1)%*%y4
      ytil5_g1 <- (diag(n)-S_g1)%*%y5
      
      xtil_g2 <- (diag(n)-S_g2)%*%x
      ytil1_g2 <- (diag(n)-S_g2)%*%y1
      ytil2_g2 <- (diag(n)-S_g2)%*%y2
      ytil3_g2 <- (diag(n)-S_g2)%*%y3
      ytil4_g2 <- (diag(n)-S_g2)%*%y4
      ytil5_g2 <- (diag(n)-S_g2)%*%y5
      
      betahat1_g1 <-solve(t(xtil_g1)%*%xtil_g1)%*%t(xtil_g1)%*%ytil1_g1 
      betahat2_g1 <-solve(t(xtil_g1)%*%xtil_g1)%*%t(xtil_g1)%*%ytil2_g1
      betahat3_g1 <-solve(t(xtil_g1)%*%xtil_g1)%*%t(xtil_g1)%*%ytil3_g1
      betahat4_g1 <-solve(t(xtil_g1)%*%xtil_g1)%*%t(xtil_g1)%*%ytil4_g1
      betahat5_g1 <-solve(t(xtil_g1)%*%xtil_g1)%*%t(xtil_g1)%*%ytil5_g1
      
      betahat1_g2 <-solve(t(xtil_g2)%*%xtil_g2)%*%t(xtil_g2)%*%ytil1_g2 
      betahat2_g2 <-solve(t(xtil_g2)%*%xtil_g2)%*%t(xtil_g2)%*%ytil2_g2
      betahat3_g2 <-solve(t(xtil_g2)%*%xtil_g2)%*%t(xtil_g2)%*%ytil3_g2
      betahat4_g2 <-solve(t(xtil_g2)%*%xtil_g2)%*%t(xtil_g2)%*%ytil4_g2
      betahat5_g2 <-solve(t(xtil_g2)%*%xtil_g2)%*%t(xtil_g2)%*%ytil5_g2
      
      beta1hat <- rowMeans(betahat1_g1,betahat1_g2)
      beta2hat <- rowMeans(betahat2_g1,betahat2_g2)
      beta3hat <- rowMeans(betahat3_g1,betahat3_g2)
      beta4hat <- rowMeans(betahat4_g1,betahat4_g2)
      beta5hat <- rowMeans(betahat5_g1,betahat5_g2)
      
      g1hat1 <- W1%*%(y1-x%*%beta1hat)
      g1hat2 <- W1%*%(y2-x%*%beta2hat)
      g1hat3 <- W1%*%(y3-x%*%beta3hat)
      g1hat4 <- W1%*%(y4-x%*%beta4hat)
      g1hat5 <- W1%*%(y5-x%*%beta5hat)
      
      g2hat1 <- W2%*%(y1-x%*%beta1hat)
      g2hat2 <- W2%*%(y2-x%*%beta2hat)
      g2hat3 <- W2%*%(y3-x%*%beta3hat)
      g2hat4 <- W2%*%(y4-x%*%beta4hat)
      g2hat5 <- W2%*%(y5-x%*%beta5hat)
      
      p0 <- tr(S_g1)
      p1 <- tr(S_g2)
      
      GCV_g1h1[i] <-gcvfunc(g1,g1hat1,p0)
      GCV_g1h2[i] <-gcvfunc(g1,g1hat2,p0)
      GCV_g1h3[i] <-gcvfunc(g1,g1hat3,p0)
      GCV_g1h4[i] <-gcvfunc(g1,g1hat4,p0)
      GCV_g1h5[i] <-gcvfunc(g1,g1hat5,p0)
      
      GCV_g2h1[i] <-gcvfunc(g2,g2hat1,p1)
      GCV_g2h2[i] <-gcvfunc(g2,g2hat2,p1)
      GCV_g2h3[i] <-gcvfunc(g2,g2hat3,p1)
      GCV_g2h4[i] <-gcvfunc(g2,g2hat4,p1)
      GCV_g2h5[i] <-gcvfunc(g2,g2hat5,p1)
    }
    
    ming <- min(scale(GCV_g1h1))
    maxg <- max(scale(GCV_g1h1))
    plot(bw_seq,scale(GCV_g1h1),type="l",ylab="GCV score",xlab="bandwidth",ylim=c(ming,maxg),col=2,lty=1,main="Bandwidth Selection for g(t1)")
    par(new=TRUE)
    plot(bw_seq,scale(GCV_g1h2),type="l",ylab="GCV score",xlab="bandwidth",ylim=c(ming,maxg),col=3,lty=2)
    par(new=TRUE)
    plot(bw_seq,scale(GCV_g1h3),type="l",ylab="GCV score",xlab="bandwidth",ylim=c(ming,maxg),col=4,lty=3)
    par(new=TRUE)
    plot(bw_seq,scale(GCV_g1h4),type="l",ylab="GCV score",xlab="bandwidth",ylim=c(ming,maxg),col=5,lty=4)
    par(new=TRUE)
    plot(bw_seq,scale(GCV_g1h5),type="l",ylab="GCV score",xlab="bandwidth",ylim=c(ming,maxg),col=6,lty=5)
    grid()
    legend(x = "topright",legend = c("GCV(g1) for y1", "GCV(g1) for y2","GCV(g1) for y3","GCV(g1) for y4","GCV(g1) for y5"), lty = c(1, 2,3,4,5),col = c(2, 3,4,5,6),lwd = 1)
    
    ming <- min(scale(GCV_g2h1))
    maxg <- max(scale(GCV_g2h1))
    plot(bw_seq,scale(GCV_g2h1),type="l",ylab="GCV score",xlab="bandwidth",ylim=c(ming,maxg),col=2,lty=1,main="Bandwidth Selection for g(t2)")
    par(new=TRUE)
    plot(bw_seq,scale(GCV_g2h2),type="l",ylab="GCV score",xlab="bandwidth",ylim=c(ming,maxg),col=3,lty=2)
    par(new=TRUE)
    plot(bw_seq,scale(GCV_g2h3),type="l",ylab="GCV score",xlab="bandwidth",ylim=c(ming,maxg),col=4,lty=3)
    par(new=TRUE)
    plot(bw_seq,scale(GCV_g2h4),type="l",ylab="GCV score",xlab="bandwidth",ylim=c(ming,maxg),col=5,lty=4)
    par(new=TRUE)
    plot(bw_seq,scale(GCV_g2h5),type="l",ylab="GCV score",xlab="bandwidth",ylim=c(ming,maxg),col=6,lty=5)
    grid()
    legend(x = "topright",legend = c("GCV(g2) for y1", "GCV(g2) for y2","GCV(g2) for y3","GCV(g2) for y4","GCV(g2) for y5"), lty = c(1, 2,3,4,5),col = c(2, 3,4,5,6),lwd = 1)
    
    for (i2 in 1:LO){
      if (GCV_g1h1[i2]==min(GCV_g1h1)){
        bwg1h1 <- bw_seq[i2]
      }
      if (GCV_g1h2[i2]==min(GCV_g1h2)){
        bwg1h2 <- bw_seq[i2]
      }
      if (GCV_g1h3[i2]==min(GCV_g1h3)){
        bwg1h3 <- bw_seq[i2]
      }
      if (GCV_g1h4[i2]==min(GCV_g1h4)){
        bwg1h4 <- bw_seq[i2]
      }
      if (GCV_g1h5[i2]==min(GCV_g1h5)){
        bwg1h5 <- bw_seq[i2]
      }
      
      
      if (GCV_g2h1[i2]==min(GCV_g2h1)){
        bwg2h1 <- bw_seq[i2]
      }
      if (GCV_g2h2[i2]==min(GCV_g2h2)){
        bwg2h2 <- bw_seq[i2]
      }
      if (GCV_g2h3[i2]==min(GCV_g2h3)){
        bwg2h3 <- bw_seq[i2]
      }
      if (GCV_g2h4[i2]==min(GCV_g2h4)){
        bwg2h4 <- bw_seq[i2]
      }
      if (GCV_g2h5[i2]==min(GCV_g2h5)){
        bwg2h5 <- bw_seq[i2]
      }
    }
  }
  #ESTIMATION PART---------------------------------------------------------------
  if(no.y==3){
    ones   <-matrix(1,n,1) 
    W1h1     <-matrix(0,n,n)
    W1h2     <-matrix(0,n,n)
    W1h3     <-matrix(0,n,n)
    
    W2h1     <-matrix(0,n,n)
    W2h2     <-matrix(0,n,n)
    W2h3     <-matrix(0,n,n)
    for (j in 1:n){
      W1h1[,j] <- NWW(t1,index1[j],bw=bwg1h1)
      W1h2[,j] <- NWW(t1,index1[j],bw=bwg1h2)
      W1h3[,j] <- NWW(t1,index1[j],bw=bwg1h3)
      
      W2h1[,j] <- NWW(t2,index2[j],bw=bwg2h1)
      W2h2[,j] <- NWW(t2,index2[j],bw=bwg2h2)
      W2h3[,j] <- NWW(t2,index2[j],bw=bwg2h3)
    }
    
    S_g1h1  <- W1h1
    S_g1h2  <- W1h2
    S_g1h3  <- W1h3
    
    S_g2h1  <- W2h1
    S_g2h2  <- W2h2
    S_g2h3  <- W2h3
    #------------------------------------------------------------------------------
    xtil1_g1  <- (diag(n)-S_g1h1)%*%x
    xtil2_g1  <- (diag(n)-S_g1h2)%*%x
    xtil3_g1  <- (diag(n)-S_g1h3)%*%x
    
    ytil1_g1 <- (diag(n)-S_g1h1)%*%y1
    ytil2_g1 <- (diag(n)-S_g1h2)%*%y2
    ytil3_g1 <- (diag(n)-S_g1h3)%*%y3
    
    xtil1_g2  <- (diag(n)-S_g2h1)%*%x
    xtil2_g2  <- (diag(n)-S_g2h2)%*%x
    xtil3_g2  <- (diag(n)-S_g2h3)%*%x
    
    ytil1_g2 <- (diag(n)-S_g2h1)%*%y1
    ytil2_g2 <- (diag(n)-S_g2h2)%*%y2
    ytil3_g2 <- (diag(n)-S_g2h3)%*%y3
    
    betahat1_g1 <-solve(t(xtil1_g1)%*%xtil1_g1)%*%t(xtil1_g1)%*%ytil1_g1 
    betahat2_g1 <-solve(t(xtil2_g1)%*%xtil2_g1)%*%t(xtil2_g1)%*%ytil2_g1
    betahat3_g1 <-solve(t(xtil3_g1)%*%xtil3_g1)%*%t(xtil3_g1)%*%ytil3_g1
    
    betahat1_g2 <-solve(t(xtil1_g2)%*%xtil1_g2)%*%t(xtil1_g2)%*%ytil1_g2 
    betahat2_g2 <-solve(t(xtil2_g2)%*%xtil2_g2)%*%t(xtil2_g2)%*%ytil2_g2
    betahat3_g2 <-solve(t(xtil3_g2)%*%xtil3_g2)%*%t(xtil3_g2)%*%ytil3_g2
    
    beta1hat <- rowMeans(betahat1_g1,betahat1_g2)
    beta2hat <- rowMeans(betahat2_g1,betahat2_g2)
    beta3hat <- rowMeans(betahat3_g1,betahat3_g2)
    
    g1hat1 <- W1h1%*%(y1-x%*%beta1hat)
    g1hat2 <- W1h2%*%(y2-x%*%beta2hat)
    g1hat3 <- W1h3%*%(y3-x%*%beta3hat)
    
    g2hat1 <- W2h1%*%(y1-x%*%beta1hat)
    g2hat2 <- W2h2%*%(y2-x%*%beta2hat)
    g2hat3 <- W2h3%*%(y3-x%*%beta3hat)
    
    y1hat <- x%*%beta1hat+g1hat1+g2hat1
    y2hat <- x%*%beta2hat+g1hat2+g2hat2
    y3hat <- x%*%beta3hat+g1hat3+g2hat3
    
    #plot(g1hat1,type="l",ylim=c(min(g1),max(g1)))
    #par(new=TRUE)
    #plot(g1hat2,type="l",ylim=c(min(g1),max(g1)))
    #par(new=TRUE)
    #plot(g1hat3,type="l",ylim=c(min(g1),max(g1)))
    #par(new=TRUE)
    #plot(g1,type="l",col="red")
    
    #plot(g2hat1,type="l",ylim=c(min(g2),max(g2)))
    #par(new=TRUE)
    #plot(g2hat2,type="l",ylim=c(min(g2),max(g2)))
    #par(new=TRUE)
    #plot(g2hat3,type="l",ylim=c(min(g2),max(g2)))
    #par(new=TRUE)
    #plot(g2,type="l",col="red")
    
    res <- new.env()
    
    res$beta1hat <- beta1hat
    res$beta2hat <- beta2hat
    res$beta3hat <- beta3hat
    #------------------------------------------------------------------------------    
    res$g1hat1 <- g1hat1
    res$g1hat2 <- g1hat2
    res$g1hat3 <- g1hat3
    
    res$g2hat1 <- g2hat1
    res$g2hat2 <- g2hat2
    res$g2hat3 <- g2hat3
    
    #------------------------------------------------------------------------------
    res$y1hat <- y1hat
    res$y2hat <- y2hat
    res$y3hat <- y3hat
    
  }
  if(no.y==5){
    ones     <-matrix(1,n,1) 
    W1h1     <-matrix(0,n,n)
    W1h2     <-matrix(0,n,n)
    W1h3     <-matrix(0,n,n)
    W1h4     <-matrix(0,n,n)
    W1h5     <-matrix(0,n,n)
    
    W2h1     <-matrix(0,n,n)
    W2h2     <-matrix(0,n,n)
    W2h3     <-matrix(0,n,n)
    W2h4     <-matrix(0,n,n)
    W2h5     <-matrix(0,n,n)
    
    for (j in 1:n){
      W1h1[,j] <- NWW(t1,index1[j],bw=bwg1h1)
      W1h2[,j] <- NWW(t1,index1[j],bw=bwg1h2)
      W1h3[,j] <- NWW(t1,index1[j],bw=bwg1h3)
      W1h4[,j] <- NWW(t1,index1[j],bw=bwg1h4)
      W1h5[,j] <- NWW(t1,index1[j],bw=bwg1h5)
      
      W2h1[,j] <- NWW(t2,index2[j],bw=bwg2h1)
      W2h2[,j] <- NWW(t2,index2[j],bw=bwg2h2)
      W2h3[,j] <- NWW(t2,index2[j],bw=bwg2h3)
      W2h4[,j] <- NWW(t2,index2[j],bw=bwg2h4)
      W2h5[,j] <- NWW(t2,index2[j],bw=bwg2h5)
    }
    
    S_g1h1  <- W1h1
    S_g1h2  <- W1h2
    S_g1h3  <- W1h3
    S_g1h4  <- W1h4
    S_g1h5  <- W1h5
    
    S_g2h1  <- W2h1
    S_g2h2  <- W2h2
    S_g2h3  <- W2h3
    S_g2h4  <- W2h4
    S_g2h5  <- W2h5
    #------------------------------------------------------------------------------
    xtil1_g1  <- (diag(n)-S_g1h1)%*%x
    xtil2_g1  <- (diag(n)-S_g1h2)%*%x
    xtil3_g1  <- (diag(n)-S_g1h3)%*%x
    xtil4_g1  <- (diag(n)-S_g1h4)%*%x
    xtil5_g1  <- (diag(n)-S_g1h5)%*%x
    
    ytil1_g1 <- (diag(n)-S_g1h1)%*%y1
    ytil2_g1 <- (diag(n)-S_g1h2)%*%y2
    ytil3_g1 <- (diag(n)-S_g1h3)%*%y3
    ytil4_g1 <- (diag(n)-S_g1h4)%*%y4
    ytil5_g1 <- (diag(n)-S_g1h5)%*%y5
    
    xtil1_g2  <- (diag(n)-S_g2h1)%*%x
    xtil2_g2  <- (diag(n)-S_g2h2)%*%x
    xtil3_g2  <- (diag(n)-S_g2h3)%*%x
    xtil4_g2  <- (diag(n)-S_g2h4)%*%x
    xtil5_g2  <- (diag(n)-S_g2h5)%*%x
    
    ytil1_g2 <- (diag(n)-S_g2h1)%*%y1
    ytil2_g2 <- (diag(n)-S_g2h2)%*%y2
    ytil3_g2 <- (diag(n)-S_g2h3)%*%y3
    ytil4_g2 <- (diag(n)-S_g2h4)%*%y4
    ytil5_g2 <- (diag(n)-S_g2h5)%*%y5
    
    betahat1_g1 <-solve(t(xtil1_g1)%*%xtil1_g1)%*%t(xtil1_g1)%*%ytil1_g1 
    betahat2_g1 <-solve(t(xtil2_g1)%*%xtil2_g1)%*%t(xtil2_g1)%*%ytil2_g1
    betahat3_g1 <-solve(t(xtil3_g1)%*%xtil3_g1)%*%t(xtil3_g1)%*%ytil3_g1
    betahat4_g1 <-solve(t(xtil4_g1)%*%xtil4_g1)%*%t(xtil4_g1)%*%ytil4_g1
    betahat5_g1 <-solve(t(xtil5_g1)%*%xtil5_g1)%*%t(xtil5_g1)%*%ytil5_g1
    
    betahat1_g2 <-solve(t(xtil1_g2)%*%xtil1_g2)%*%t(xtil1_g2)%*%ytil1_g2 
    betahat2_g2 <-solve(t(xtil2_g2)%*%xtil2_g2)%*%t(xtil2_g2)%*%ytil2_g2
    betahat3_g2 <-solve(t(xtil3_g2)%*%xtil3_g2)%*%t(xtil3_g2)%*%ytil3_g2
    betahat4_g2 <-solve(t(xtil4_g2)%*%xtil4_g2)%*%t(xtil4_g2)%*%ytil4_g2
    betahat5_g2 <-solve(t(xtil5_g2)%*%xtil5_g2)%*%t(xtil5_g2)%*%ytil5_g2
    
    beta1hat <- rowMeans(betahat1_g1,betahat1_g2)
    beta2hat <- rowMeans(betahat2_g1,betahat2_g2)
    beta3hat <- rowMeans(betahat3_g1,betahat3_g2)
    beta4hat <- rowMeans(betahat4_g1,betahat4_g2)
    beta5hat <- rowMeans(betahat5_g1,betahat5_g2)
    
    g1hat1 <- W1h1%*%(y1-x%*%beta1hat)
    g1hat2 <- W1h2%*%(y2-x%*%beta2hat)
    g1hat3 <- W1h3%*%(y3-x%*%beta3hat)
    g1hat4 <- W1h4%*%(y4-x%*%beta4hat)
    g1hat5 <- W1h5%*%(y5-x%*%beta5hat)
    
    g2hat1 <- W2h1%*%(y1-x%*%beta1hat)
    g2hat2 <- W2h2%*%(y2-x%*%beta2hat)
    g2hat3 <- W2h3%*%(y3-x%*%beta3hat)
    g2hat4 <- W2h4%*%(y4-x%*%beta4hat)
    g2hat5 <- W2h5%*%(y5-x%*%beta5hat)
    
    y1hat <- x%*%beta1hat+g1hat1+g2hat1
    y2hat <- x%*%beta2hat+g1hat2+g2hat2
    y3hat <- x%*%beta3hat+g1hat3+g2hat3
    y4hat <- x%*%beta4hat+g1hat4+g2hat4
    y5hat <- x%*%beta5hat+g1hat5+g2hat5

    
    res <- new.env()
    
    res$beta1hat <- beta1hat
    res$beta2hat <- beta2hat
    res$beta3hat <- beta3hat
    res$beta4hat <- beta4hat
    res$beta5hat <- beta5hat
    #------------------------------------------------------------------------------    
    res$g1hat1 <- g1hat1
    res$g1hat2 <- g1hat2
    res$g1hat3 <- g1hat3
    res$g1hat4 <- g1hat4
    res$g1hat5 <- g1hat5
    res$g2hat1 <- g2hat1
    res$g2hat2 <- g2hat2
    res$g2hat3 <- g2hat3
    res$g2hat4 <- g2hat4
    res$g2hat5 <- g2hat5
    #------------------------------------------------------------------------------
    res$y1hat <- y1hat
    res$y2hat <- y2hat
    res$y3hat <- y3hat
    res$y4hat <- y4hat
    res$y5hat <- y5hat
  }
  #-------------------------------------------------------------------------------
  return(res)
  #-----------------------------------------------------------------------------
  #-----------------------------------------------------------------------------
}
#-------------------------------------------------------------------------------
sim  <- 500
n    <- 200 #(50, 100, 200)
CL   <- 0.30#(5%, 15%, 30%)
no.z <- 5 #(3,5)
#ZERO MATRICES------------------------------------------------------------------
LP.beta1hat <- matrix(0,3,sim)
LP.beta2hat <- matrix(0,3,sim)
LP.beta3hat <- matrix(0,3,sim)
LP.beta4hat <- matrix(0,3,sim)
LP.beta5hat <- matrix(0,3,sim)

LL.beta1hat <- matrix(0,3,sim)
LL.beta2hat <- matrix(0,3,sim)
LL.beta3hat <- matrix(0,3,sim)
LL.beta4hat <- matrix(0,3,sim)
LL.beta5hat <- matrix(0,3,sim)

LC.beta1hat <- matrix(0,3,sim)
LC.beta2hat <- matrix(0,3,sim)
LC.beta3hat <- matrix(0,3,sim)
LC.beta4hat <- matrix(0,3,sim)
LC.beta5hat <- matrix(0,3,sim)

LP.g1hat1 <- matrix(0,n,sim)
LP.g1hat2 <- matrix(0,n,sim)
LP.g1hat3 <- matrix(0,n,sim)
LP.g1hat4 <- matrix(0,n,sim)
LP.g1hat5 <- matrix(0,n,sim)

LL.g1hat1 <- matrix(0,n,sim)
LL.g1hat2 <- matrix(0,n,sim)
LL.g1hat3 <- matrix(0,n,sim)
LL.g1hat4 <- matrix(0,n,sim)
LL.g1hat5 <- matrix(0,n,sim)

LC.g1hat1 <- matrix(0,n,sim)
LC.g1hat2 <- matrix(0,n,sim)
LC.g1hat3 <- matrix(0,n,sim)
LC.g1hat4 <- matrix(0,n,sim)
LC.g1hat5 <- matrix(0,n,sim)

LP.g2hat1 <- matrix(0,n,sim)
LP.g2hat2 <- matrix(0,n,sim)
LP.g2hat3 <- matrix(0,n,sim)
LP.g2hat4 <- matrix(0,n,sim)
LP.g2hat5 <- matrix(0,n,sim)

LL.g2hat1 <- matrix(0,n,sim)
LL.g2hat2 <- matrix(0,n,sim)
LL.g2hat3 <- matrix(0,n,sim)
LL.g2hat4 <- matrix(0,n,sim)
LL.g2hat5 <- matrix(0,n,sim)

LC.g2hat1 <- matrix(0,n,sim)
LC.g2hat2 <- matrix(0,n,sim)
LC.g2hat3 <- matrix(0,n,sim)
LC.g2hat4 <- matrix(0,n,sim)
LC.g2hat5 <- matrix(0,n,sim)

LP.y1hat <- matrix(0,n,sim)
LP.y2hat <- matrix(0,n,sim)
LP.y3hat <- matrix(0,n,sim)
LP.y4hat <- matrix(0,n,sim)
LP.y5hat <- matrix(0,n,sim)

LL.y1hat <- matrix(0,n,sim)
LL.y2hat <- matrix(0,n,sim)
LL.y3hat <- matrix(0,n,sim)
LL.y4hat <- matrix(0,n,sim)
LL.y5hat <- matrix(0,n,sim)

LC.y1hat <- matrix(0,n,sim)
LC.y2hat <- matrix(0,n,sim)
LC.y3hat <- matrix(0,n,sim)
LC.y4hat <- matrix(0,n,sim)
LC.y5hat <- matrix(0,n,sim)

LP.bias1 <- matrix(0,3,sim)
LP.bias2 <- matrix(0,3,sim)
LP.bias3 <- matrix(0,3,sim)
LP.bias4 <- matrix(0,3,sim)
LP.bias5 <- matrix(0,3,sim)

LL.bias1 <- matrix(0,3,sim)
LL.bias2 <- matrix(0,3,sim)
LL.bias3 <- matrix(0,3,sim)
LL.bias4 <- matrix(0,3,sim)
LL.bias5 <- matrix(0,3,sim)

LC.bias1 <- matrix(0,3,sim)
LC.bias2 <- matrix(0,3,sim)
LC.bias3 <- matrix(0,3,sim)
LC.bias4 <- matrix(0,3,sim)
LC.bias5 <- matrix(0,3,sim)

#-------------------------------------------------------------------------------
for (s in 1:sim){
  set.seed(184375+s)
  data <- datagen(n,CL,no.z)
  x <- data$x
  y <- data$Y
  nc <- data$nc
  allf <- data$allf
  rbeta <- data$allbeta
  
  res1 <- MRLP(x,y,nc,allf)
  res2 <- MRLL(x,y,nc,allf)
  res3 <- MRLC(x,y,nc,allf)

if (no.z==3){
LP.beta1hat[,s] <- res1$beta1hat
LP.beta2hat[,s] <- res1$beta2hat
LP.beta3hat[,s] <- res1$beta3hat

LL.beta1hat[,s] <- res2$beta1hat
LL.beta2hat[,s] <- res2$beta2hat
LL.beta3hat[,s] <- res2$beta3hat

LC.beta1hat[,s] <- res3$beta1hat
LC.beta2hat[,s] <- res3$beta2hat
LC.beta3hat[,s] <- res3$beta3hat

LP.g1hat1[,s] <- res1$g1hat1
LP.g1hat2[,s] <- res1$g1hat2
LP.g1hat3[,s] <- res1$g1hat3

LL.g1hat1[,s] <- res2$g1hat1
LL.g1hat2[,s] <- res2$g1hat2
LL.g1hat3[,s] <- res2$g1hat3

LC.g1hat1[,s] <- res3$g1hat1
LC.g1hat2[,s] <- res3$g1hat2
LC.g1hat3[,s] <- res3$g1hat3

LP.g2hat1[,s] <- res1$g2hat1
LP.g2hat2[,s] <- res1$g2hat2
LP.g2hat3[,s] <- res1$g2hat3

LL.g2hat1[,s] <- res2$g2hat1
LL.g2hat2[,s] <- res2$g2hat2
LL.g2hat3[,s] <- res2$g2hat3

LC.g2hat1[,s] <- res3$g2hat1
LC.g2hat2[,s] <- res3$g2hat2
LC.g2hat3[,s] <- res3$g2hat3

LP.y1hat[,s] <- res1$y1hat
LP.y2hat[,s] <- res1$y2hat
LP.y3hat[,s] <- res1$y3hat

LL.y1hat[,s] <- res2$y1hat
LL.y2hat[,s] <- res2$y2hat
LL.y3hat[,s] <- res2$y3hat

LC.y1hat[,s] <- res3$y1hat
LC.y2hat[,s] <- res3$y2hat
LC.y3hat[,s] <- res3$y3hat

LP.bias1[,s] <- res1$beta1hat-data$allbeta[,1]
LP.bias2[,s] <- res1$beta2hat-data$allbeta[,2]
LP.bias3[,s] <- res1$beta3hat-data$allbeta[,3]

LL.bias1[,s] <- res2$beta1hat-data$allbeta[,1]
LL.bias2[,s] <- res2$beta2hat-data$allbeta[,2]
LL.bias3[,s] <- res2$beta3hat-data$allbeta[,3]

LC.bias1[,s] <- res3$beta1hat-data$allbeta[,1]
LC.bias2[,s] <- res3$beta2hat-data$allbeta[,2]
LC.bias3[,s] <- res3$beta3hat-data$allbeta[,3]

}
  if (no.z==5){
    LP.beta1hat[,s] <- res1$beta1hat
    LP.beta2hat[,s] <- res1$beta2hat
    LP.beta3hat[,s] <- res1$beta3hat
    LP.beta4hat[,s] <- res1$beta4hat
    LP.beta5hat[,s] <- res1$beta5hat
    
    LL.beta1hat[,s] <- res2$beta1hat
    LL.beta2hat[,s] <- res2$beta2hat
    LL.beta3hat[,s] <- res2$beta3hat
    LL.beta4hat[,s] <- res2$beta4hat
    LL.beta5hat[,s] <- res2$beta5hat
    
    LC.beta1hat[,s] <- res3$beta1hat
    LC.beta2hat[,s] <- res3$beta2hat
    LC.beta3hat[,s] <- res3$beta3hat
    LC.beta4hat[,s] <- res3$beta4hat
    LC.beta5hat[,s] <- res3$beta5hat
    
    LP.g1hat1[,s] <- res1$g1hat1
    LP.g1hat2[,s] <- res1$g1hat2
    LP.g1hat3[,s] <- res1$g1hat3
    LP.g1hat4[,s] <- res1$g1hat4
    LP.g1hat5[,s] <- res1$g1hat5
    
    LL.g1hat1[,s] <- res2$g1hat1
    LL.g1hat2[,s] <- res2$g1hat2
    LL.g1hat3[,s] <- res2$g1hat3
    LL.g1hat4[,s] <- res2$g1hat4
    LL.g1hat5[,s] <- res2$g1hat5
    
    LC.g1hat1[,s] <- res3$g1hat1
    LC.g1hat2[,s] <- res3$g1hat2
    LC.g1hat3[,s] <- res3$g1hat3
    LC.g1hat4[,s] <- res3$g1hat4
    LC.g1hat5[,s] <- res3$g1hat5
    
    LP.g2hat1[,s] <- res1$g2hat1
    LP.g2hat2[,s] <- res1$g2hat2
    LP.g2hat3[,s] <- res1$g2hat3
    LP.g2hat4[,s] <- res1$g2hat4
    LP.g2hat5[,s] <- res1$g2hat5
    
    LL.g2hat1[,s] <- res2$g2hat1
    LL.g2hat2[,s] <- res2$g2hat2
    LL.g2hat3[,s] <- res2$g2hat3
    LL.g2hat4[,s] <- res2$g2hat4
    LL.g2hat5[,s] <- res2$g2hat5
    
    LC.g2hat1[,s] <- res3$g2hat1
    LC.g2hat2[,s] <- res3$g2hat2
    LC.g2hat3[,s] <- res3$g2hat3
    LC.g2hat4[,s] <- res3$g2hat4
    LC.g2hat5[,s] <- res3$g2hat5
    
    LP.y1hat[,s] <- res1$y1hat
    LP.y2hat[,s] <- res1$y2hat
    LP.y3hat[,s] <- res1$y3hat
    LP.y4hat[,s] <- res1$y4hat
    LP.y5hat[,s] <- res1$y5hat
    
    LL.y1hat[,s] <- res2$y1hat
    LL.y2hat[,s] <- res2$y2hat
    LL.y3hat[,s] <- res2$y3hat
    LL.y4hat[,s] <- res2$y4hat
    LL.y5hat[,s] <- res2$y5hat
    
    LC.y1hat[,s] <- res3$y1hat
    LC.y2hat[,s] <- res3$y2hat
    LC.y3hat[,s] <- res3$y3hat
    LC.y4hat[,s] <- res3$y4hat
    LC.y5hat[,s] <- res3$y5hat
    
    LP.bias1[,s] <- res1$beta1hat-data$allbeta[,1]
    LP.bias2[,s] <- res1$beta2hat-data$allbeta[,2]
    LP.bias3[,s] <- res1$beta3hat-data$allbeta[,3]
    LP.bias4[,s] <- res1$beta4hat-data$allbeta[,4]
    LP.bias5[,s] <- res1$beta5hat-data$allbeta[,5]
    
    LL.bias1[,s] <- res2$beta1hat-data$allbeta[,1]
    LL.bias2[,s] <- res2$beta2hat-data$allbeta[,2]
    LL.bias3[,s] <- res2$beta3hat-data$allbeta[,3]
    LL.bias4[,s] <- res2$beta4hat-data$allbeta[,4]
    LL.bias5[,s] <- res2$beta5hat-data$allbeta[,5]
    
    LC.bias1[,s] <- res3$beta1hat-data$allbeta[,1]
    LC.bias2[,s] <- res3$beta2hat-data$allbeta[,2]
    LC.bias3[,s] <- res3$beta3hat-data$allbeta[,3]
    LC.bias4[,s] <- res3$beta4hat-data$allbeta[,4]
    LC.bias5[,s] <- res3$beta5hat-data$allbeta[,5]
  }
message(s,"th simulation ends")  
}
se.beta1 <- matrix(c(sd(LP.beta1hat[1,]),sd(LP.beta1hat[2,]),sd(LP.beta1hat[3,]),sd(LL.beta1hat[1,]),sd(LL.beta1hat[2,]),sd(LL.beta1hat[3,]),sd(LC.beta1hat[1,]),sd(LC.beta1hat[2,]),sd(LC.beta1hat[3,])),3,3)
colnames(se.beta1) <- c("LP","LL","LC")

se.beta2 <- matrix(c(sd(LP.beta2hat[1,]),sd(LP.beta2hat[2,]),sd(LP.beta2hat[3,]),sd(LL.beta2hat[1,]),sd(LL.beta2hat[2,]),sd(LL.beta2hat[3,]),sd(LC.beta2hat[1,]),sd(LC.beta2hat[2,]),sd(LC.beta2hat[3,])),3,3)
colnames(se.beta2) <- c("LP","LL","LC")

se.beta3 <- matrix(c(sd(LP.beta3hat[1,]),sd(LP.beta3hat[2,]),sd(LP.beta3hat[3,]),sd(LL.beta3hat[1,]),sd(LL.beta3hat[2,]),sd(LL.beta3hat[3,]),sd(LC.beta3hat[1,]),sd(LC.beta3hat[2,]),sd(LC.beta3hat[3,])),3,3)
colnames(se.beta3) <- c("LP","LL","LC")
if (no.z==5){
  se.beta4 <- matrix(c(sd(LP.beta4hat[1,]),sd(LP.beta4hat[2,]),sd(LP.beta4hat[3,]),sd(LL.beta4hat[1,]),sd(LL.beta4hat[2,]),sd(LL.beta4hat[3,]),sd(LC.beta4hat[1,]),sd(LC.beta4hat[2,]),sd(LC.beta4hat[3,])),3,3)
  colnames(se.beta4) <- c("LP","LL","LC")
  
  se.beta5 <- matrix(c(sd(LP.beta5hat[1,]),sd(LP.beta5hat[2,]),sd(LP.beta5hat[3,]),sd(LL.beta5hat[1,]),sd(LL.beta5hat[2,]),sd(LL.beta5hat[3,]),sd(LC.beta5hat[1,]),sd(LC.beta5hat[2,]),sd(LC.beta5hat[3,])),3,3)
  colnames(se.beta5) <- c("LP","LL","LC")
}

BIAS_beta1 <- matrix(c(mean(LP.bias1[1,]),mean(LP.bias1[2,]),mean(LP.bias1[3,]),mean(LL.bias1[1,]),mean(LL.bias1[2,]),mean(LL.bias1[3,]),mean(LC.bias1[1,]),mean(LC.bias1[2,]),mean(LC.bias1[3,])),3,3)
BIAS_beta2 <- matrix(c(mean(LP.bias2[1,]),mean(LP.bias2[2,]),mean(LP.bias2[3,]),mean(LL.bias2[1,]),mean(LL.bias2[2,]),mean(LL.bias2[3,]),mean(LC.bias2[1,]),mean(LC.bias2[2,]),mean(LC.bias2[3,])),3,3)
BIAS_beta3 <- matrix(c(mean(LP.bias3[1,]),mean(LP.bias3[2,]),mean(LP.bias3[3,]),mean(LL.bias3[1,]),mean(LL.bias3[2,]),mean(LL.bias3[3,]),mean(LC.bias3[1,]),mean(LC.bias3[2,]),mean(LC.bias3[3,])),3,3)
if (no.z==5){
  BIAS_beta4 <- matrix(c(mean(LP.bias4[1,]),mean(LP.bias4[2,]),mean(LP.bias4[3,]),mean(LL.bias4[1,]),mean(LL.bias4[2,]),mean(LL.bias4[3,]),mean(LC.bias4[1,]),mean(LC.bias4[2,]),mean(LC.bias4[3,])),3,3)
  BIAS_beta5 <- matrix(c(mean(LP.bias5[1,]),mean(LP.bias5[2,]),mean(LP.bias5[3,]),mean(LL.bias5[1,]),mean(LL.bias5[2,]),mean(LL.bias5[3,]),mean(LC.bias5[1,]),mean(LC.bias5[2,]),mean(LC.bias5[3,])),3,3)
}
#TRADOFF PLOTS------------------------------------------------------------------(LP)
library(ggplot2)
library(gridExtra)
ind <- c(1:3)
s.var1  <- scale(se.beta1[,1]^2)
s.bias1 <- scale(BIAS_beta1[,1]) 
tradeoff1 <-data.frame(s.var1,s.bias1) 
tradeoff1[order(tradeoff1$s.var1, decreasing = FALSE), ] 
dfto1 <- data.frame(ind,s.var1)
dfto2 <- data.frame(ind,s.bias1)
p1 <- ggplot()+geom_smooth(data=dfto1,aes(x=ind,y=s.var1,col="variance(y1)"))+geom_smooth(data=dfto2,aes(x=ind,y=s.bias1,col="Bias(y1)"))+xlab("LP Estimates")+ylab("Var.& Bias")+ theme(legend.position="bottom") +theme(legend.title=element_blank())
p1
s.var2  <- scale(se.beta2[,1]^2)
s.bias2 <- scale(BIAS_beta2[,1]) 
tradeoff2 <-data.frame(s.var2,s.bias2) 
tradeoff2[order(tradeoff2$s.var2, decreasing = FALSE), ] 
dfto3 <- data.frame(ind,s.var2)
dfto4 <- data.frame(ind,s.bias2)
p2 <- ggplot()+geom_smooth(data=dfto3,aes(x=ind,y=s.var2,col="variance(y2)"))+geom_smooth(data=dfto4,aes(x=ind,y=s.bias2,col="Bias(y2)"))+xlab("LP Estimates")+ylab("Var.& Bias")+ theme(legend.position="bottom") +theme(legend.title=element_blank())
p2
s.var3  <- scale(se.beta3[,1]^2)
s.bias3 <- scale(BIAS_beta3[,1]) 
tradeoff3 <-data.frame(s.var3,s.bias3) 
tradeoff3[order(tradeoff3$s.var3, decreasing = FALSE), ] 
dfto5 <- data.frame(ind,s.var3)
dfto6 <- data.frame(ind,s.bias3)
p3 <- ggplot()+geom_smooth(data=dfto5,aes(x=ind,y=s.var3,col="variance(y3)"))+geom_smooth(data=dfto6,aes(x=ind,y=s.bias3,col="Bias(y3)"))+xlab("LP Estimates")+ylab("Var.& Bias")+ theme(legend.position="bottom") +theme(legend.title=element_blank())
p3
if (no.z==5){
s.var4  <- scale(se.beta4[,1]^2)
s.bias4 <- scale(BIAS_beta4[,1]) 
tradeoff4 <-data.frame(s.var4,s.bias4) 
tradeoff4[order(tradeoff4$s.var4, decreasing = FALSE), ] 
dfto45 <- data.frame(ind,s.var4)
dfto46 <- data.frame(ind,s.bias4)
p43 <- ggplot()+geom_smooth(data=dfto45,aes(x=ind,y=s.var4,col="variance(y4)"))+geom_smooth(data=dfto46,aes(x=ind,y=s.bias4,col="Bias(y4)"))+xlab("LP Estimates")+ylab("Var.& Bias")+ theme(legend.position="bottom") +theme(legend.title=element_blank())
p43
}
#------------------------------------------------------------------------------(LL)
s2.var1  <- scale(se.beta1[,2]^2)
s2.bias1 <- scale(BIAS_beta1[,2]) 
tradeoff21 <-data.frame(s2.var1,s2.bias1) 
tradeoff21[order(tradeoff21$s2.var1, decreasing = FALSE), ] 
dfto21 <- data.frame(ind,s2.var1)
dfto22 <- data.frame(ind,s2.bias1)
p11 <- ggplot()+geom_smooth(data=dfto21,aes(x=ind,y=s2.var1,col="variance(y1)"))+geom_smooth(data=dfto22,aes(x=ind,y=s2.bias1,col="Bias(y1)"))+xlab("LL Estimates")+ylab("Var.& Bias")+ theme(legend.position="bottom") +theme(legend.title=element_blank())
p11
s2.var2  <- scale(se.beta2[,2]^2)
s2.bias2 <- scale(BIAS_beta2[,2]) 
tradeoff22 <-data.frame(s2.var2,s2.bias2) 
tradeoff22[order(tradeoff22$s2.var2, decreasing = FALSE), ] 
dfto23 <- data.frame(ind,s2.var2)
dfto24 <- data.frame(ind,s2.bias2)
p22 <- ggplot()+geom_smooth(data=dfto23,aes(x=ind,y=s2.var2,col="variance(y2)"))+geom_smooth(data=dfto24,aes(x=ind,y=s2.bias2,col="Bias(y2)"))+xlab("LL Estimates")+ylab("Var.& Bias")+ theme(legend.position="bottom") +theme(legend.title=element_blank())
p22
s2.var3  <- scale(se.beta3[,2]^2)
s2.bias3 <- scale(BIAS_beta3[,2]) 
tradeoff23 <-data.frame(s2.var3,s2.bias3) 
tradeoff23[order(tradeoff23$s2.var3, decreasing = FALSE), ] 
dfto25 <- data.frame(ind,s2.var3)
dfto26 <- data.frame(ind,s2.bias3)
p33 <- ggplot()+geom_smooth(data=dfto25,aes(x=ind,y=s2.var3,col="variance(y3)"))+geom_smooth(data=dfto26,aes(x=ind,y=s2.bias3,col="Bias(y3)"))+xlab("LL Estimates")+ylab("Var.& Bias")+ theme(legend.position="bottom") +theme(legend.title=element_blank())
p33
if (no.z==5){
s2.var5  <- scale(se.beta5[,2]^2)
s2.bias5 <- scale(BIAS_beta5[,2]) 
tradeoff523 <-data.frame(s2.var5,s2.bias5) 
tradeoff523[order(tradeoff523$s2.var5, decreasing = FALSE), ] 
dfto525 <- data.frame(ind,s2.var5)
dfto526 <- data.frame(ind,s2.bias5)
p533 <- ggplot()+geom_smooth(data=dfto525,aes(x=ind,y=s2.var5,col="variance(y5)"))+geom_smooth(data=dfto526,aes(x=ind,y=s2.bias5,col="Bias(y5)"))+xlab("LL Estimates")+ylab("Var.& Bias")+ theme(legend.position="bottom") +theme(legend.title=element_blank())
p533
}
#------------------------------------------------------------------------------(LC)
s3.var1  <- scale(se.beta1[,3]^2)
s3.bias1 <- scale(BIAS_beta1[,3]) 
tradeoff31 <-data.frame(s3.var1,s3.bias1) 
tradeoff31[order(tradeoff31$s3.var1, decreasing = FALSE), ] 
dfto31 <- data.frame(ind,s3.var1)
dfto32 <- data.frame(ind,s3.bias1)
p111 <- ggplot()+geom_smooth(data=dfto31,aes(x=ind,y=s3.var1,col="variance(y1)"))+geom_smooth(data=dfto32,aes(x=ind,y=s3.bias1,col="Bias(y1)"))+xlab("LC Estimates")+ylab("Var.& Bias")+ theme(legend.position="bottom") +theme(legend.title=element_blank())
p111
s3.var2  <- scale(se.beta2[,3]^2)
s3.bias2 <- scale(BIAS_beta2[,3]) 
tradeoff32 <-data.frame(s3.var2,s3.bias2) 
tradeoff32[order(tradeoff32$s3.var2, decreasing = FALSE), ] 
dfto33 <- data.frame(ind,s3.var2)
dfto34 <- data.frame(ind,s3.bias2)
p222 <- ggplot()+geom_smooth(data=dfto33,aes(x=ind,y=s3.var2,col="variance(y2)"))+geom_smooth(data=dfto34,aes(x=ind,y=s3.bias2,col="Bias(y2)"))+xlab("LC Estimates")+ylab("Var.& Bias")+ theme(legend.position="bottom") +theme(legend.title=element_blank())
p222
s3.var3  <- scale(se.beta3[,3]^2)
s3.bias3 <- scale(BIAS_beta3[,3]) 
tradeoff33 <-data.frame(s3.var3,s3.bias3) 
tradeoff33[order(tradeoff33$s3.var3, decreasing = FALSE), ] 
dfto35 <- data.frame(ind,s3.var3)
dfto36 <- data.frame(ind,s3.bias3)
p333 <- ggplot()+geom_smooth(data=dfto35,aes(x=ind,y=s3.var3,col="variance(y3)"))+geom_smooth(data=dfto36,aes(x=ind,y=s3.bias3,col="Bias(y3)"))+xlab("LC Estimates")+ylab("Var.& Bias")+ theme(legend.position="bottom") +theme(legend.title=element_blank())
p333
if (no.z==5){
s3.var4  <- scale(se.beta4[,3]^2)
s3.bias4 <- scale(BIAS_beta4[,3]) 
tradeoff433 <-data.frame(s3.var4,s3.bias4) 
tradeoff433[order(tradeoff433$s3.var4, decreasing = FALSE), ] 
dfto435 <- data.frame(ind,s3.var4)
dfto436 <- data.frame(ind,s3.bias4)
p4333 <- ggplot()+geom_smooth(data=dfto435,aes(x=ind,y=s3.var4,col="variance(y4)"))+geom_smooth(data=dfto436,aes(x=ind,y=s3.bias4,col="Bias(y4)"))+xlab("LC Estimates")+ylab("Var.& Bias")+ theme(legend.position="bottom") +theme(legend.title=element_blank())
p4333
}
#grid.arrange(p1,p2,p3,p11,p22,p33,p111,p222,p333,ncol=2)
#grid.arrange(p1,p22,p333,p2,p11,p111,nrow=3)
grid.arrange(p1,p43, p22,p533,p333,p4333,nrow=3)
#BIAS PLOTS----------------------------------------------------------------------
par(mar=c(4,4,1,1))
Tbias_y1 <- matrix(c(t(LP.bias1),t(LL.bias1),t(LC.bias1)),sim,9)
boxplot(Tbias_y1,names=c("B1","B2","B3","B1","B2","B3","B1","B2","B3"),xlab="Coefficients",ylab="Bias")
abline(a=NULL,b=NULL,h=NULL,v=3.5, lty="dotted")
abline(a=NULL,b=NULL,h=NULL,v=6.5, lty="dotted")
text(2,1.5,"Bias (LP)")
text(5,1.5,"Bias (LL)")
text(8,1.5,"Bias (LC)")
#-------------------------------------------------------------------------------
Tbias_y2 <- matrix(c(t(LP.bias2),t(LL.bias2),t(LC.bias2)),sim,9)
boxplot(Tbias_y2,names=c("B1","B2","B3","B1","B2","B3","B1","B2","B3"),xlab="Coefficients",ylab="Bias")
abline(a=NULL,b=NULL,h=NULL,v=3.5, lty="dotted")
abline(a=NULL,b=NULL,h=NULL,v=6.5, lty="dotted")
text(2,1.5,"Bias (LP)")
text(5,1.5,"Bias (LL)")
text(8,1.5,"Bias (LC)")
#-------------------------------------------------------------------------------
Tbias_y3 <- matrix(c(t(LP.bias3),t(LL.bias3),t(LC.bias3)),sim,9)
boxplot(Tbias_y3,names=c("B1","B2","B3","B1","B2","B3","B1","B2","B3"),xlab="Coefficients",ylab="Bias")
abline(a=NULL,b=NULL,h=NULL,v=3.5, lty="dotted")
abline(a=NULL,b=NULL,h=NULL,v=6.5, lty="dotted")
text(2,1.5,"Bias (LP)")
text(5,1.5,"Bias (LL)")
text(8,1.5,"Bias (LC)")
#-------------------------------------------------------------------------------
if (no.z==5){
  Tbias_y4 <- matrix(c(t(LP.bias4),t(LL.bias4),t(LC.bias4)),sim,9)
  boxplot(Tbias_y4,names=c("B1","B2","B3","B1","B2","B3","B1","B2","B3"),xlab="Coefficients",ylab="Bias")
  abline(a=NULL,b=NULL,h=NULL,v=3.5, lty="dotted")
  abline(a=NULL,b=NULL,h=NULL,v=6.5, lty="dotted")
  text(2,1.3,"Bias (LP)")
  text(5,1.3,"Bias (LL)")
  text(8,1.3,"Bias (LC)")
  #-------------------------------------------------------------------------------
  Tbias_y5 <- matrix(c(t(LP.bias5),t(LL.bias5),t(LC.bias5)),sim,9)
  boxplot(Tbias_y5,names=c("B1","B2","B3","B1","B2","B3","B1","B2","B3"),xlab="Coefficients",ylab="Bias")
  abline(a=NULL,b=NULL,h=NULL,v=3.5, lty="dotted")
  abline(a=NULL,b=NULL,h=NULL,v=6.5, lty="dotted")
  text(2,1.2,"Bias (LP)")
  text(5,1.2,"Bias (LL)")
  text(8,1.2,"Bias (LC)")
}
LP.G1HAT1 <- rowMeans(LP.g1hat1)
LP.G1HAT2 <- rowMeans(LP.g1hat2)
LP.G1HAT3 <- rowMeans(LP.g1hat3)
if (no.z==3){
ming <- min(data$allf[,1])
maxg <- max(data$allf[,1])
for (j in 1:sim){
  plot(data$nc[,1],LP.g1hat3[,j],type="l",ylab="g(t1)",xlab="t1",ylim=c(ming,maxg),col = rgb(red = 0, green = 0, blue = 0, alpha = 0.04))
  par(new=TRUE)
}
par(new=TRUE)
plot(data$nc[,1],data$allf[,1],type="l",ylab="g(t1)",xlab="t1",ylim=c(ming,maxg),col=1,lty=1,main="LP fits of g(t1)")
par(new=TRUE)
plot(data$nc[,1],LP.G1HAT1,type="l",ylab="g(t1)",xlab="t1",ylim=c(ming,maxg),col=2,lty=2,lwd=2)
par(new=TRUE)
plot(data$nc[,1],LP.G1HAT2,type="l",ylab="g(t1)",xlab="t1",ylim=c(ming,maxg),col=3,lty=3,lwd=2)
par(new=TRUE)
plot(data$nc[,1],LP.G1HAT3,type="l",ylab="g(t1)",xlab="t1",ylim=c(ming,maxg),col=4,lty=4,lwd=2)
grid()
legend(x = "topright",legend = c("Real g(t1)", "g(t1)->y1","g(t1)->y2","g(t1)->y3"), lty = c(1,2,3,4),col = c(1,2,3,4),lwd = 1)


par(new=FALSE)
}
if (no.z==5){
  LP.G1HAT4 <- rowMeans(LP.g1hat4)
  LP.G1HAT5 <- rowMeans(LP.g1hat5)
  
  ming <- min(data$allf[,1])
  maxg <- max(data$allf[,1])
  for (j in 1:sim){
    plot(data$nc[,1],LP.g1hat3[,j],type="l",ylab="g(t1)",xlab="t1",ylim=c(ming,maxg),col = rgb(red = 0, green = 0, blue = 0, alpha = 0.04))
    par(new=TRUE)
  }
  par(new=TRUE)
  plot(data$nc[,1],data$allf[,1],type="l",ylab="g(t1)",xlab="t1",ylim=c(ming,maxg),col=1,lty=1,main="LP fits of g(t1)")
  par(new=TRUE)
  plot(data$nc[,1],LP.G1HAT1,type="l",ylab="g(t1)",xlab="t1",ylim=c(ming,maxg),col=2,lty=2,lwd=2)
  par(new=TRUE)
  plot(data$nc[,1],LP.G1HAT2,type="l",ylab="g(t1)",xlab="t1",ylim=c(ming,maxg),col=3,lty=3,lwd=2)
  par(new=TRUE)
  plot(data$nc[,1],LP.G1HAT3,type="l",ylab="g(t1)",xlab="t1",ylim=c(ming,maxg),col=4,lty=4,lwd=2)
  par(new=TRUE)
  plot(data$nc[,1],LP.G1HAT4,type="l",ylab="g(t1)",xlab="t1",ylim=c(ming,maxg),col=5,lty=5,lwd=2)
  par(new=TRUE)
  plot(data$nc[,1],LP.G1HAT5,type="l",ylab="g(t1)",xlab="t1",ylim=c(ming,maxg),col=6,lty=6,lwd=2)
  grid()
  legend(x = "topright",legend = c("Real g(t1)", "g(t1)->y1","g(t1)->y2","g(t1)->y3","g(t1)->y4","g(t1)->y5"), lty = c(1, 2,3,4,5,6),col = c(1,2,3,4,5,6),lwd = 1)
  
  
  par(new=FALSE)
  }

LP.G2HAT1 <- rowMeans(LP.g2hat1)
LP.G2HAT2 <- rowMeans(LP.g2hat2)
LP.G2HAT3 <- rowMeans(LP.g2hat3)

if (no.z==3){
ming <- min(data$allf[,2])
maxg <- max(data$allf[,2])
for (j in 1:sim){
  plot(data$nc[,2],LP.g2hat3[,j],type="l",ylab="g(t2)",xlab="t2",ylim=c(ming,maxg),col = rgb(red = 0, green = 0, blue = 0, alpha = 0.04))
  par(new=TRUE)
}
par(new=TRUE)
plot(data$nc[,2],data$allf[,2],type="l",ylab="g(t2)",xlab="t2",ylim=c(ming,maxg),col=1,lty=1,main="LP fits of g(t2)")
par(new=TRUE)
plot(data$nc[,2],LP.G2HAT1,type="l",ylab="g(t2)",xlab="t2",ylim=c(ming,maxg),col=2,lty=2,lwd=2)
par(new=TRUE)
plot(data$nc[,2],LP.G2HAT2,type="l",ylab="g(t2)",xlab="t2",ylim=c(ming,maxg),col=3,lty=3,lwd=2)
par(new=TRUE)
plot(data$nc[,2],LP.G2HAT3,type="l",ylab="g(t2)",xlab="t2",ylim=c(ming,maxg),col=4,lty=4,lwd=2)
grid()
legend(x = "topright",legend = c("Real g(t2)", "g(t2)->y1","g(t2)->y2","g(t2)->y3"), lty = c(1, 2,3,4),col = c(1,2,3,4),lwd = 1)


par(new=FALSE)
}
if (no.z==5){
  LP.G2HAT4 <- rowMeans(LP.g2hat4)
  LP.G2HAT5 <- rowMeans(LP.g2hat5)
  
  ming <- min(data$allf[,2])
  maxg <- max(data$allf[,2])
  for (j in 1:sim){
    plot(data$nc[,2],LP.g2hat3[,j],type="l",ylab="g(t2)",xlab="t2",ylim=c(ming,maxg),col = rgb(red = 0, green = 0, blue = 0, alpha = 0.04))
    par(new=TRUE)
  }
  par(new=TRUE)
  plot(data$nc[,2],data$allf[,2],type="l",ylab="g(t2)",xlab="t2",ylim=c(ming,maxg),col=1,lty=1,main="LP fits of g(t2)")
  par(new=TRUE)
  plot(data$nc[,2],LP.G2HAT1,type="l",ylab="g(t2)",xlab="t2",ylim=c(ming,maxg),col=2,lty=2,lwd=2)
  par(new=TRUE)
  plot(data$nc[,2],LP.G2HAT2,type="l",ylab="g(t2)",xlab="t2",ylim=c(ming,maxg),col=3,lty=3,lwd=2)
  par(new=TRUE)
  plot(data$nc[,2],LP.G2HAT3,type="l",ylab="g(t2)",xlab="t2",ylim=c(ming,maxg),col=4,lty=4,lwd=2)
  par(new=TRUE)
  plot(data$nc[,2],LP.G2HAT4,type="l",ylab="g(t2)",xlab="t2",ylim=c(ming,maxg),col=5,lty=5,lwd=2)
  par(new=TRUE)
  plot(data$nc[,2],LP.G2HAT5,type="l",ylab="g(t2)",xlab="t2",ylim=c(ming,maxg),col=6,lty=6,lwd=2)
  grid()
  legend(x = "topright",legend = c("Real g(t2)", "g(t2)->y1","g(t2)->y2","g(t2)->y3","g(t2)->y4","g(t2)->y5"), lty = c(1,2,3,4,5,6),col = c(1,2,3,4,5,6),lwd = 1)
  
  par(new=FALSE)
  }

LL.G1HAT1 <- rowMeans(LL.g1hat1)
LL.G1HAT2 <- rowMeans(LL.g1hat2)
LL.G1HAT3 <- rowMeans(LL.g1hat3)

if (no.z==3){
ming <- min(data$allf[,1])
maxg <- max(data$allf[,1])
for (j in 1:sim){
  plot(data$nc[,1],LL.g1hat3[,j],type="l",ylab="g(t1)",xlab="t1",ylim=c(ming,maxg),col = rgb(red = 0, green = 0, blue = 0, alpha = 0.04))
  par(new=TRUE)
}
par(new=TRUE)
plot(data$nc[,1],data$allf[,1],type="l",ylab="g(t1)",xlab="t1",ylim=c(ming,maxg),col=1,lty=1,main="LL fits of g(t1)")
par(new=TRUE)
plot(data$nc[,1],LL.G1HAT1,type="l",ylab="g(t1)",xlab="t1",ylim=c(ming,maxg),col=2,lty=2,lwd=2)
par(new=TRUE)
plot(data$nc[,1],LL.G1HAT2,type="l",ylab="g(t1)",xlab="t1",ylim=c(ming,maxg),col=3,lty=3,lwd=2)
par(new=TRUE)
plot(data$nc[,1],LL.G1HAT3,type="l",ylab="g(t1)",xlab="t1",ylim=c(ming,maxg),col=4,lty=4,lwd=2)
grid()
legend(x = "topright",legend = c("Real g(t1)", "g(t1)->y1","g(t1)->y2","g(t1)->y3"), lty = c(1, 2,3,4),col = c(1,2,3,4),lwd = 1)


par(new=FALSE)
}
if (no.z==5){
  LL.G1HAT4 <- rowMeans(LL.g1hat4)
  LL.G1HAT5 <- rowMeans(LL.g1hat5)
  
  ming <- min(data$allf[,1])
  maxg <- max(data$allf[,1])
  for (j in 1:sim){
    plot(data$nc[,1],LL.g1hat3[,j],type="l",ylab="g(t1)",xlab="t1",ylim=c(ming,maxg),col = rgb(red = 0, green = 0, blue = 0, alpha = 0.04))
    par(new=TRUE)
  }
  par(new=TRUE)
  plot(data$nc[,1],data$allf[,1],type="l",ylab="g(t1)",xlab="t1",ylim=c(ming,maxg),col=1,lty=1,main="LL fits of g(t1)")
  par(new=TRUE)
  plot(data$nc[,1],LL.G1HAT1,type="l",ylab="g(t1)",xlab="t1",ylim=c(ming,maxg),col=2,lty=2,lwd=2)
  par(new=TRUE)
  plot(data$nc[,1],LL.G1HAT2,type="l",ylab="g(t1)",xlab="t1",ylim=c(ming,maxg),col=3,lty=3,lwd=2)
  par(new=TRUE)
  plot(data$nc[,1],LL.G1HAT3,type="l",ylab="g(t1)",xlab="t1",ylim=c(ming,maxg),col=4,lty=4,lwd=2)
  par(new=TRUE)
  plot(data$nc[,1],LL.G1HAT4,type="l",ylab="g(t1)",xlab="t1",ylim=c(ming,maxg),col=5,lty=5,lwd=2)
  par(new=TRUE)
  plot(data$nc[,1],LL.G1HAT5,type="l",ylab="g(t1)",xlab="t1",ylim=c(ming,maxg),col=6,lty=6,lwd=2)
  grid()
  legend(x = "topright",legend = c("Real g(t1)", "g(t1)->y1","g(t1)->y2","g(t1)->y3","g(t1)->y4","g(t1)->y5"), lty = c(1, 2,3,4,5,6),col = c(1,2,3,4,5,6),lwd = 1)
  
  par(new=FALSE)
  }

LL.G2HAT1 <- rowMeans(LL.g2hat1)
LL.G2HAT2 <- rowMeans(LL.g2hat2)
LL.G2HAT3 <- rowMeans(LL.g2hat3)
if (no.z==3){
ming <- min(data$allf[,2])
maxg <- max(data$allf[,2])
for (j in 1:sim){
  plot(data$nc[,2],LL.g2hat3[,j],type="l",ylab="g(t2)",xlab="t2",ylim=c(ming,maxg),col = rgb(red = 0, green = 0, blue = 0, alpha = 0.04))
  par(new=TRUE)
}
par(new=TRUE)
plot(data$nc[,2],data$allf[,2],type="l",ylab="g(t2)",xlab="t2",ylim=c(ming,maxg),col=1,lty=1,main="LL fits of g(t2)")
par(new=TRUE)
plot(data$nc[,2],LL.G2HAT1,type="l",ylab="g(t2)",xlab="t2",ylim=c(ming,maxg),col=2,lty=2,lwd=2)
par(new=TRUE)
plot(data$nc[,2],LL.G2HAT2,type="l",ylab="g(t2)",xlab="t2",ylim=c(ming,maxg),col=3,lty=3,lwd=2)
par(new=TRUE)
plot(data$nc[,2],LL.G2HAT3,type="l",ylab="g(t2)",xlab="t2",ylim=c(ming,maxg),col=4,lty=4,lwd=2)
grid()
legend(x = "topright",legend = c("Real g(t2)", "g(t2)->y1","g(t2)->y2","g(t2)->y3"), lty = c(1, 2,3,4),col = c(1,2,3,4),lwd = 1)

par(new=FALSE)
}
if (no.z==5){
  LL.G2HAT4 <- rowMeans(LL.g2hat4)
  LL.G2HAT5 <- rowMeans(LL.g2hat5)
  
  ming <- min(data$allf[,2])
  maxg <- max(data$allf[,2])
  
  for (j in 1:sim){
    plot(data$nc[,2],LL.g2hat3[,j],type="l",ylab="g(t2)",xlab="t2",ylim=c(ming,maxg),col = rgb(red = 0, green = 0, blue = 0, alpha = 0.04))
    par(new=TRUE)
  }
  par(new=TRUE)
  plot(data$nc[,2],data$allf[,2],type="l",ylab="g(t2)",xlab="t2",ylim=c(ming,maxg),col=1,lty=1,main="LL fits of g(t2)")
  par(new=TRUE)
  plot(data$nc[,2],LL.G2HAT1,type="l",ylab="g(t2)",xlab="t2",ylim=c(ming,maxg),col=2,lty=2,lwd=2)
  par(new=TRUE)
  plot(data$nc[,2],LL.G2HAT2,type="l",ylab="g(t2)",xlab="t2",ylim=c(ming,maxg),col=3,lty=3,lwd=2)
  par(new=TRUE)
  plot(data$nc[,2],LL.G2HAT3,type="l",ylab="g(t2)",xlab="t2",ylim=c(ming,maxg),col=4,lty=4,lwd=2)
  par(new=TRUE)
  plot(data$nc[,2],LL.G2HAT4,type="l",ylab="g(t2)",xlab="t2",ylim=c(ming,maxg),col=5,lty=5,lwd=2)
  par(new=TRUE)
  plot(data$nc[,2],LL.G2HAT5,type="l",ylab="g(t2)",xlab="t2",ylim=c(ming,maxg),col=6,lty=6,lwd=2)
  grid()
  legend(x = "topright",legend = c("Real g(t2)", "g(t2)->y1","g(t2)->y2","g(t2)->y3","g(t2)->y4","g(t2)->y5"), lty = c(1, 2,3,4,5,6),col = c(1,2,3,4,5,6),lwd = 1)
  
  par(new=FALSE)
  }

LC.G1HAT1 <- rowMeans(LC.g1hat1)
LC.G1HAT2 <- rowMeans(LC.g1hat2)
LC.G1HAT3 <- rowMeans(LC.g1hat3)
if (no.z==3){
ming <- min(data$allf[,1])
maxg <- max(data$allf[,1])
for (j in 1:sim){
  plot(data$nc[,1],LC.g1hat3[,j],type="l",ylab="g(t1)",xlab="t1",ylim=c(ming,maxg),col = rgb(red = 0, green = 0, blue = 0, alpha = 0.04))
  par(new=TRUE)
}
par(new=TRUE)
plot(data$nc[,1],data$allf[,1],type="l",ylab="g(t1)",xlab="t1",ylim=c(ming,maxg),col=1,lty=1,main="LC fits of g(t1)")
par(new=TRUE)
plot(data$nc[,1],LC.G1HAT1,type="l",ylab="g(t1)",xlab="t1",ylim=c(ming,maxg),col=2,lty=2,lwd=2)
par(new=TRUE)
plot(data$nc[,1],LC.G1HAT2,type="l",ylab="g(t1)",xlab="t1",ylim=c(ming,maxg),col=3,lty=3,lwd=2)
par(new=TRUE)
plot(data$nc[,1],LC.G1HAT3,type="l",ylab="g(t1)",xlab="t1",ylim=c(ming,maxg),col=4,lty=4,lwd=2)
grid()
legend(x = "topright",legend = c("Real g(t1)", "g(t1)->y1","g(t1)->y2","g(t1)->y3"), lty = c(1, 2,3,4),col = c(1,2,3,4),lwd = 1)
par(new=FALSE)
}
if (no.z==5){
  LC.G1HAT4 <- rowMeans(LC.g1hat4)
  LC.G1HAT5 <- rowMeans(LC.g1hat5)
  
  ming <- min(data$allf[,1])
  maxg <- max(data$allf[,1])
  for (j in 1:sim){
    plot(data$nc[,1],LC.g1hat3[,j],type="l",ylab="g(t1)",xlab="t1",ylim=c(ming,maxg),col = rgb(red = 0, green = 0, blue = 0, alpha = 0.04))
    par(new=TRUE)
  }
  par(new=TRUE)
  plot(data$nc[,1],data$allf[,1],type="l",ylab="g(t1)",xlab="t1",ylim=c(ming,maxg),col=1,lty=1,main="LC fits of g(t1)")
  par(new=TRUE)
  plot(data$nc[,1],LC.G1HAT1,type="l",ylab="g(t1)",xlab="t1",ylim=c(ming,maxg),col=2,lty=2,lwd=2)
  par(new=TRUE)
  plot(data$nc[,1],LC.G1HAT2,type="l",ylab="g(t1)",xlab="t1",ylim=c(ming,maxg),col=3,lty=3,lwd=2)
  par(new=TRUE)
  plot(data$nc[,1],LC.G1HAT3,type="l",ylab="g(t1)",xlab="t1",ylim=c(ming,maxg),col=4,lty=4,lwd=2)
  par(new=TRUE)
  plot(data$nc[,1],LC.G1HAT4,type="l",ylab="g(t1)",xlab="t1",ylim=c(ming,maxg),col=5,lty=5,lwd=2)
  par(new=TRUE)
  plot(data$nc[,1],LC.G1HAT5,type="l",ylab="g(t1)",xlab="t1",ylim=c(ming,maxg),col=6,lty=6,lwd=2)
  grid()
  legend(x = "topright",legend = c("Real g(t1)", "g(t1)->y1","g(t1)->y2","g(t1)->y3","g(t1)->y4","g(t1)->y5"), lty = c(1, 2,3,4,5,6),col = c(1,2,3,4,5,6),lwd = 1)
  par(new=FALSE)
}

LC.G2HAT1 <- rowMeans(LC.g2hat1)
LC.G2HAT2 <- rowMeans(LC.g2hat2)
LC.G2HAT3 <- rowMeans(LC.g2hat3)
if (no.z==3){
ming <- min(data$allf[,2])
maxg <- max(data$allf[,2])
for (j in 1:sim){
  plot(data$nc[,2],LC.g2hat3[,j],type="l",ylab="g(t2)",xlab="t2",ylim=c(ming,maxg),col = rgb(red = 0, green = 0, blue = 0, alpha = 0.04))
  par(new=TRUE)
}
par(new=TRUE)
plot(data$nc[,2],data$allf[,2],type="l",ylab="g(t2)",xlab="t2",ylim=c(ming,maxg),col=1,lty=1,main="LC fits of g(t2)")
par(new=TRUE)
plot(data$nc[,2],LC.G2HAT1,type="l",ylab="g(t2)",xlab="t2",ylim=c(ming,maxg),col=2,lty=2,lwd=2)
par(new=TRUE)
plot(data$nc[,2],LC.G2HAT2,type="l",ylab="g(t2)",xlab="t2",ylim=c(ming,maxg),col=3,lty=3,lwd=2)
par(new=TRUE)
plot(data$nc[,2],LC.G2HAT3,type="l",ylab="g(t2)",xlab="t2",ylim=c(ming,maxg),col=4,lty=4,lwd=2)
grid()
legend(x = "topright",legend = c("Real g(t2)", "g(t2)->y1","g(t2)->y2","g(t2)->y3"), lty = c(1, 2,3,4),col = c(1,2,3,4),lwd = 1)
par(new=FALSE)
}
if (no.z==5){
  LC.G2HAT4 <- rowMeans(LC.g2hat4)
  LC.G2HAT5 <- rowMeans(LC.g2hat5)
  
  ming <- min(data$allf[,2])
  maxg <- max(data$allf[,2])
  for (j in 1:sim){
    plot(data$nc[,2],LC.g2hat3[,j],type="l",ylab="g(t2)",xlab="t2",ylim=c(ming,maxg),col = rgb(red = 0, green = 0, blue = 0, alpha = 0.04))
    par(new=TRUE)
  }
  par(new=TRUE)
  plot(data$nc[,2],data$allf[,2],type="l",ylab="g(t2)",xlab="t2",ylim=c(ming,maxg),col=1,lty=1,main="LC fits of g(t2)")
  par(new=TRUE)
  plot(data$nc[,2],LC.G2HAT1,type="l",ylab="g(t2)",xlab="t2",ylim=c(ming,maxg),col=2,lty=2,lwd=2)
  par(new=TRUE)
  plot(data$nc[,2],LC.G2HAT2,type="l",ylab="g(t2)",xlab="t2",ylim=c(ming,maxg),col=3,lty=3,lwd=2)
  par(new=TRUE)
  plot(data$nc[,2],LC.G2HAT3,type="l",ylab="g(t2)",xlab="t2",ylim=c(ming,maxg),col=4,lty=4,lwd=2)
  par(new=TRUE)
  plot(data$nc[,2],LC.G2HAT4,type="l",ylab="g(t2)",xlab="t2",ylim=c(ming,maxg),col=5,lty=5,lwd=2)
  par(new=TRUE)
  plot(data$nc[,2],LC.G2HAT5,type="l",ylab="g(t2)",xlab="t2",ylim=c(ming,maxg),col=6,lty=6,lwd=2)
  grid()
  legend(x = "topright",legend = c("Real g(t2)", "g(t2)->y1","g(t2)->y2","g(t2)->y3","g(t2)->y4","g(t2)->y5"), lty = c(1, 2,3,4,5,6),col = c(1,2,3,4,5,6),lwd = 1)
  par(new=FALSE)
  }
#-------------------------------------------------------------------------------

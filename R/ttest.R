# t-test

sc.cf <- function(Y1,Y0,T1,T0,K,lsei_type){

  T01   <- T0 + T1
  r     <- min(floor(T0/K),T1)

  Y1.pre  <- Y1[1:T0]
  Y0.pre  <- Y0[1:T0,]
  Y1.post <- Y1[(T0+1):T01]
  Y0.post <- Y0[(T0+1):T01,]

  tau.mat <- matrix(NA,K,1)
  for (k in 1:K){
    Hk            <- (T0-(r*K))+seq((k-1)*r+1,k*r,1)
    w.Hk          <- sc(Y1.pre[-Hk],Y0.pre[-Hk,],lsei_type=lsei_type)$w.hat
    tau.mat[k,1]  <- mean(Y1.post-Y0.post%*%w.Hk) - mean(Y1.pre[Hk]-Y0.pre[Hk,]%*%w.Hk)
  }

  tau.hat <- mean(tau.mat)
  se.hat  <- sqrt(1+((K*r)/T1))*sd(tau.mat)/sqrt(K)
  t.hat   <- tau.hat/se.hat # this has a t_{K-1} distribution

  return(list(t.hat=t.hat,tau.hat=tau.hat,se.hat=se.hat))
}

did.cf <- function(Y1,Y0,T1,T0,K){

  T01   <- T0 + T1
  r     <- min(floor(T0/K),T1)

  Y1.pre  <- Y1[1:T0]
  Y0.pre  <- Y0[1:T0,]
  Y1.post <- Y1[(T0+1):T01]
  Y0.post <- Y0[(T0+1):T01,]

  tau.mat <- matrix(NA,K,1)
  for (k in 1:K){
    Hk            <- (T0-(r*K))+seq((k-1)*r+1,k*r,1)
    tau.mat[k,1]  <- mean(Y1.post-rowMeans(Y0.post)) - mean(Y1.pre[Hk]-rowMeans(Y0.pre[Hk,]))
  }

  tau.hat <- mean(tau.mat)
  se.hat  <- sqrt(1+((K*r)/T1))*sd(tau.mat)/sqrt(K)
  t.hat   <- tau.hat/se.hat # this has a t_{K-1} distribution

  return(list(t.hat=t.hat,tau.hat=tau.hat,se.hat=se.hat))
}


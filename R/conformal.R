# Moving block permutations
movingblock <- function(Y1,Y0,T1,T0,theta0,estimation_method,lsei_type){
  T01 <- T0+T1

  Y1_0 <- Y1
  Y1_0[(T0+1):T01] <- Y1[(T0+1):T01]-theta0

  if (estimation_method=="classo"){
    u.hat <- classo(Y1_0,Y0)$u.hat
  }
  if (estimation_method=="sc"){
    u.hat <- sc(Y1_0,Y0,lsei_type)$u.hat
  }
  if (estimation_method=="did"){
    u.hat <- did(Y1_0,Y0)$u.hat
  }
  sub.size  <- T1
  u.hat.c   <- c(u.hat,u.hat)
  S.vec     <- matrix(NA,T01,1)
  for (s in 1:(T01)){
    S.vec[s,1]  <- sum(abs(u.hat.c[s:(s+sub.size-1)]))
  }
  p <- mean(S.vec>=S.vec[T0+1])
  return(p)
}

# All/iid permutations (use random subsample of size n_perm all permutations)
iid <- function(Y1,Y0,T1,T0,theta0,estimation_method,n_perm,lsei_type){
  T01 <- T0+T1

  Y1_0 <- Y1
  Y1_0[(T0+1):T01] <- Y1[(T0+1):T01]-theta0

  if (estimation_method=="classo"){
    u.hat <- classo(Y1_0,Y0)$u.hat
  }
  if (estimation_method=="sc"){
    u.hat <- sc(Y1_0,Y0,lsei_type)$u.hat
  }
  if (estimation_method=="did"){
    u.hat <- did(Y1_0,Y0)$u.hat
  }
  post.ind  <- ((T0+1):T01)
  pre.ind   <- (1:T0)
  S.vec     <- matrix(NA,n_perm,1)

  Sq <- sum(abs(u.hat[post.ind]))
  for (r in 1:n_perm){
    u.hat.p     <- u.hat[sample(1:T01,replace=F)]
    S.vec[r,1]  <- sum(abs(u.hat.p[post.ind]))
  }
  p <- 1/(n_perm+1)*(1+sum(S.vec>=Sq))
  return(p)
}

# Confidence interval via test inversion
confidence_interval <- function(Y1,Y0,T1,T0,estimation_method,alpha,ci_grid,lsei_type){

  lb <- ub <- rep(NA,T1)
  for (t in 1:T1){
    indices   <- c(1:T0,T0+t)
    Y1_temp   <- Y1[indices]
    Y0_temp   <- Y0[indices,]

    ps_temp <- rep(NA,length(ci_grid))
    for (ind in 1:length(ci_grid)){
      Y1_0_temp <- Y1_temp
      Y1_0_temp[(T0+1)] <- Y1_temp[(T0+1)] - ci_grid[ind]
      if (estimation_method=="classo"){
        u_hat <- classo(Y1_0_temp,Y0_temp)$u.hat
      }
      if (estimation_method=="sc"){
        u_hat <- sc(Y1_0_temp,Y0_temp,lsei_type)$u.hat
      }
      if (estimation_method=="did"){
        u_hat <- did(Y1_0_temp,Y0_temp)$u.hat
      }
      ps_temp[ind] <- mean(abs(u_hat)>=abs(u_hat[(T0+1)]))
    }
    ci_temp <- ci_grid[ps_temp>alpha]
    lb[t]     <- min(ci_temp,na.rm=TRUE)
    ub[t]     <- max(ci_temp,na.rm=TRUE)
  }

  return(list(lb=lb,ub=ub))
}

# difference-in-differences
did <- function(Y1,Y0){
  u.hat <- Y1 - mean(Y1-rowMeans(Y0)) - rowMeans(Y0)
  return(list(u.hat=u.hat))
}

# synthetic control
sc <- function(Y1,Y0,lsei_type){
  J <- dim(Y0)[2]
  e <- matrix(1,1,J)
  f <- 1
  g <- diag(x=1,J,J)
  h <- matrix(0,J,1)
  w.hat <- limSolve::lsei(A=Y0,B=Y1,E=e,F=f,G=g,H=h,type=lsei_type)$X
  u.hat <- Y1-Y0%*%w.hat
  return(list(u.hat=u.hat,w.hat=w.hat))
}

# constrained lasso
classo <- function(Y1,Y0){
  J         <- dim(Y0)[2]
  w         <- CVXR::Variable((J+1))
  objective <- CVXR::Minimize(mean((Y1-cbind(1,Y0)%*%w)^2))
  prob      <- CVXR::Problem(objective,constraints = list(sum(abs(w[2:(J+1)])) <= 1))
  result    <- CVXR::solve(prob)
  w.hat     <- result$getValue(w)
  u.hat     <- Y1 - cbind(1,Y0)%*%w.hat
  return(list(u.hat=u.hat,w.hat=w.hat))
}


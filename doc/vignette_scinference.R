## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----install------------------------------------------------------------------
#install.packages("devtools")
#devtools::install_github("kwuthrich/scinference")
library(scinference)

## ----help---------------------------------------------------------------------
?scinference

## ----setup--------------------------------------------------------------------
set.seed(12345)

J   <- 50
T0  <- 50 
T1  <- 5

w       <- rep(0,J)
w[1:3]  <- 1/3
Y0      <- matrix(rnorm((T0+T1)*J),(T0+T1),J)
Y1      <- Y0 %*% w + rnorm(T0+T1)

Y1[(T0+1):(T0+T1)] <- Y1[(T0+1):(T0+T1)] + 2

## ----null---------------------------------------------------------------------
scinference(Y1,Y0,T1=T1,T0=T0,theta0=4,estimation_method="sc",permutation_method="mb")$p_val
scinference(Y1,Y0,T1=T1,T0=T0,theta0=4,estimation_method="did",permutation_method="mb")$p_val
scinference(Y1,Y0,T1=T1,T0=T0,theta0=4,estimation_method="classo",permutation_method="mb")$p_val

scinference(Y1,Y0,T1=T1,T0=T0,theta0=4,estimation_method="did",permutation_method="iid")$p_val
scinference(Y1,Y0,T1=T1,T0=T0,theta0=4,estimation_method="sc",permutation_method="iid")$p_val
scinference(Y1,Y0,T1=T1,T0=T0,theta0=4,estimation_method="classo",permutation_method="iid")$p_val

## ----ci-----------------------------------------------------------------------

obj <- scinference(Y1,Y0,T1=T1,T0=T0,estimation_method="sc",ci=TRUE,ci_grid=seq(-2,8,0.1))

plot(1, ylab="", xlab="Time", main="90% pointwise CIs",xlim = c(1,5), ylim=c(-1,7), type="n")
lines(1:5, obj$lb, lwd = 2)
lines(1:5, obj$ub, lwd = 2)
abline(h=0,col="darkgrey", lwd=1)
  


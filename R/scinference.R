#' Inference estimation_methods for synthetic control
#'
#' The function \code{scinference} implements the inference methods for synthetic controls proposed by Chernozhukov et al. (2020a,b).
#' The methods apply to a canonical synthetic control setup with 1 treated unit and J control units.
#' The treated unit is untreated for the first \code{T0} periods and treated for the remaining \code{T1=T-T0} periods.
#'
#' @param Y1 outcome data for treated unit (T x 1 vector)
#' @param Y0 outcome data for control units (T x J matrix)
#' @param T1 number of post-treatment periods
#' @param T0 number of pre-treatment periods
#' @param inference_method inference method; conformal inference ("conformal",default), t-test ("ttest")
#' @param alpha significance level; default \code{alpha} = 0.1
#' @param ci logical, indicating whether pointwise confidence intervals are reported; default \code{ci=FALSE}
#' @param theta0 null hypothesis for treatment effect trajectory (T1 x 1 vector or scalar if constant effects null); default \code{theta0}=0
#' @param estimation_method estimation method; difference-in-differences ("did"), synthetic control ("sc", default), constrained lasso ("classo"). Note that constrained lasso is not implemented for the t-test.
#' @param permutation_method permutation method; moving block permutations ("mb", default), iid permutations ("iid")
#' @param ci_grid grid for the confidence interval
#' @param n_perm number of permutation (relevant for iid permutations); default = 5000
#' @param lsei_type option for lsei (package limSolve) used for sc; default = 1
#' @param K K>1 number of cross-fits for t-test; default = 2
#' @return conformal inference: p-value for testing the null that \code{theta=theta0} and pointwise CI (lower bounds and upper bounds) if \code{ci=TRUE}; t-test: ATT estimate, standard error, and confidence interval

#' @export
scinference <-
  function(Y1,
           Y0,
           T1,
           T0,
           inference_method = "conformal",
           alpha = 0.1,
           ci = FALSE,
           theta0 = 0,
           estimation_method = "sc",
           permutation_method = "mb",
           ci_grid = NULL,
           n_perm = 5000,
           lsei_type = 1,
           K = 2) {

  # preliminaries
  if(length(Y1)!=(T0+T1)) stop("length of Y1 needs to be equal to T")
  if(dim(Y0)[1]!=(T0+T1)) stop("number of rows in Y0 needs to be equal to T")
  if(!(inference_method %in% c("conformal","ttest"))) stop("The selected inference method is not available")

  if(inference_method == "conformal"){

    if(!(estimation_method %in% c("did","sc","classo"))) stop("The selected estimation method is not implemented for conformal inference")
    if(!(permutation_method %in% c("iid","mb"))) stop("The selected class of permutations is not available")
    if(length(theta0)!=1){
      if(length(theta0)!=T1){
        stop("length of theta0 should be T1")
      }
    }

    # p-value for overall null hypothesis

    if(permutation_method == "mb") {
      p_val <- movingblock(Y1=Y1,Y0=Y0,T1=T1,T0=T0,theta0=theta0,estimation_method=estimation_method,lsei_type=lsei_type)
    }
    if(permutation_method == "iid") {
      p_val <- iid(Y1=Y1,Y0=Y0,T1=T1,T0=T0,theta0=theta0,estimation_method=estimation_method,n_perm=n_perm,lsei_type=lsei_type)
    }

    # pointwise confidence intervals

    if(ci==TRUE){
      if(is.null(ci_grid)){
        stop("no grid specified for confidence interval")
      }
      obj <- confidence_interval(Y1=Y1,Y0=Y0,T1=T1,T0=T0,estimation_method=estimation_method,alpha=alpha,ci_grid=ci_grid,lsei_type=lsei_type)
      ub <- obj$ub
      lb <- obj$lb
    } else {
      lb <- NA
      ub <- NA
    }
    return(list(p_val=p_val,lb=lb,ub=ub))

  }

  if(inference_method == "ttest"){

    if(!(estimation_method %in% c("did","sc"))) stop("The selected estimation method is not implemented for the t-test")
    if(K==1) stop("K must be strictly than 1")
    if(estimation_method == "did"){
      obj <- did.cf(Y1,Y0,T1,T0,K)
    }
    if(estimation_method == "sc"){
      obj <- sc.cf(Y1,Y0,T1,T0,K,lsei_type)
    }

    att <- obj$tau.hat
    se  <- obj$se.hat
    lb  <- att - qt(1-alpha/2,df=K-1)*se
    ub  <- att + qt(1-alpha/2,df=K-1)*se

    return(list(att=att,se=se,lb=lb,ub=ub))

  }

}








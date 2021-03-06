#' Inference estimation_methods for synthetic control
#'
#' The function \code{scinference} implements the inference methods for synthetic controls proposed by Chernozhukov et al. (2020). The paper is available here: https://arxiv.org/abs/1712.09089.
#' The methods apply to a canonical synthetic control setup with 1 treated unit and J control units.
#' The treated unit is untreated for the first \code{T0} periods and treated for the remaining \code{T1=T-T0} periods.
#'
#' @param Y1 outcome data for treated unit (T x 1 vector)
#' @param Y0 outcome data for control units (T x J matrix)
#' @param T1 number of post-treatment periods
#' @param T0 number of pre-treatment periods
#' @param ci logical, indicating whether pointwise confidence intervals are reported; default \code{ci=FALSE}
#' @param theta0 null hypothesis for treatment effect trajectory (T1 x 1 vector or scalar if constant effects null); default \code{theta0}=0
#' @param estimation_method estimation method; difference-in-differences ("did"), synthetic control ("sc", default), constrained lasso ("classo")
#' @param permutation_method permutation method; moving block permutations ("mb", default), iid permutations ("iid")
#' @param ci_grid grid for the confidence interval
#' @param alpha significance level; default \code{alpha} = 0.1
#' @param n_perm number of permutation (relevant for iid permutations); default = 5000
#' @param lsei_type option for lsei (package limSolve) used for sc
#' @return p-value for testing the null that \code{theta=theta0} and pointwise CI (lower bounds and upper bounds) if \code{ci=TRUE}

#' @export
scinference <-
  function(Y1,
           Y0,
           T1,
           T0,
           ci = FALSE,
           theta0 = 0,
           estimation_method = "sc",
           permutation_method = "mb",
           ci_grid = NULL,
           alpha = 0.1,
           n_perm = 5000,
           lsei_type = 2) {

  # preliminaries
  if(length(Y1)!=(T0+T1)) stop("length of Y1 needs to be equal to T")
  if(dim(Y0)[1]!=(T0+T1)) stop("number of rows in Y0 needs to be equal to T")
  if(!(estimation_method %in% c("did","sc","classo"))) stop("The selected estimation method is not available")
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








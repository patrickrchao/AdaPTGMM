#' Generate Spline Basis
#'
#' @param x one dimensional covariate
#' @param ndf number of degrees of freedom, 1 corresponds to intercept model
#'
#' @details \code{ndf}=1 corresponds to an intercept only model, meaning we disregard \code{x} as a whole.
#' We utilize a natural cubic spline basis
#' @noRd
generate_spline <- function(x,ndf){
  num_hypo <- length(x)
  # Create natural spline basis
  # Add column of ones
  ifelse(ndf > 1,
    spline <- cbind(rep(1,num_hypo), splines::ns(x, df = ndf-1)),
    spline <- matrix(rep(1,num_hypo))
    )

  return(spline)

}

#' Weighted Mean computation
#'
#' @params values Vector of values
#' @params weights Vector of nonnegative weights
#'
#' @return weighted mean
#' @noRd
weighted_mean <- function(values,weights){
  return(sum(values*weights)/sum(weights))
}

#' Perform checks for valid parameters
#'
#' @noRd
input_checks <- function(x,p_values,z,ndf,nclasses,niter,alpha_m,zeta,lambda,masking_shape,alphas){
  if(is.null(x) | (is.null(p_values) & is.null(z))){
    stop("Invalid inputs for x, p_values, and test_statistics. None were inputted.")
  }

  if(length(x) != max(length(p_values), length(z))){
    stop("Invalid inputs, x and p_values/test statistics have different length.")
  }
  if(min(p_values<0) | max(p_values) > 1){
    stop("Invalid p-values, p-values must be in range [0,1].")
  }
  if(alpha_m<0 | alpha_m>1 | alpha_m>lambda | lambda<0 | lambda>1 | zeta<0 | zeta>1 | lambda+alpha_m/zeta>1){
    stop("Invalid input for alpha_m, zeta, lambda, must all be between 0 and 1 and 0<alpha_m<=lambda<lambda+alpha_m/zeta<=1.")
  }
  if(masking_shape!="tent" & masking_shape != "comb"){
    stop("Invalid masking shape inputted, must be `tent` or `comb`.")
  }
  if(min(alphas) < 0 | max(alphas) > 1){
    stop("Invalid alphas inputted, alphas must be in range [0,1].")
  }
  if(niter < 0){
    stop("Invalid number of iterations.")
  }
  if(min(nclasses) < 2){
    stop("Invalid number of classes, minimum 2 classes.")
  }
  if(min(ndf) < 1){
    stop("Invalid number of degrees of freedom, minimum 1 (intercept only model).")
  }
}

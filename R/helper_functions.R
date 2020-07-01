#' Weighted Mean computation
#'
#' @param values Vector of values
#' @param weights Vector of nonnegative weights
#'
#' @return weighted mean
#' @noRd
.weighted_mean <- function(values,weights){
  return(sum(values*weights)/sum(weights))
}

#' Perform checks for valid parameters
#'
#' @noRd
.input_checks <- function(x,pvals,z,testing,rendpoint,lendpoint,ndf,nclasses,niter_fit,niter_ms,nfit,masking_params,masking_shape,alphas){
  alpha_m <- masking_params$alpha_m
  zeta <- masking_params$zeta
  lambda <- masking_params$lambda
  if(is.null(x) | (is.null(pvals) & is.null(z))){
    stop("Invalid inputs for x, pvals, and test_statistics. None were inputted.")
  }
  if(!is.data.frame(x)){
    stop("Invalid input, x must be a dataframe.")
  }
  if(testing == "interval"){
    if(is.null(z)){
      stop("Invalid input for z. Must include test statistics to perform interval null testing.")
    }else if(!is.null(pvals)){
      warning("Inputted p-values will be ignored in interval null testing, they will be computed automatically.")
    }
    if(is.null(rendpoint)){
      stop("Invalid input for rendpoint. Must include right endpoint for interval null testing.")
    }
  }
  #if(length(x) != max(length(pvals), length(z))){
  #  stop("Invalid inputs, x and pvals/test statistics have different length.")
  #}
  if(!is.null(pvals)){
     if((min(pvals<0) | max(pvals) > 1)){
        stop("Invalid p-values, p-values must be in range [0,1].")
     }
  }

  if(alpha_m<0 | alpha_m>1 | alpha_m>lambda | lambda<0 | lambda>1 | zeta<0  | lambda+alpha_m*zeta>1){
    stop("Invalid input for alpha_m and lambda, must all be between 0 and 1 and 0<alpha_m<=lambda<lambda+alpha_m*zeta<=1.")
  }
  if(masking_shape!="tent" & masking_shape != "comb"){
    stop("Invalid masking shape inputted, must be `tent` or `comb`.")
  }
  if(min(alphas) < 0 | max(alphas) > 1){
    stop("Invalid alphas inputted, alphas must be in range [0,1].")
  }
  if(niter_fit <= 0 | niter_ms <= 0 | nfit <= 0){
    stop("Invalid number of iterations for niter_fit, niter_ms, or nfit, must be an integer greater than 0.")
  }
  if(min(nclasses) < 2){
    stop("Invalid number of classes, minimum 2 classes.")
  }
  if(min(ndf) < 1){
    stop("Invalid number of degrees of freedom, minimum 1 (intercept only model).")
  }
  if(!is.null(rendpoint) & !is.null(lendpoint)){
    if(lendpoint > rendpoint){
      stop("Invalid input for endpoints. Right endpoint must be greater than left endpoint.")
    }
  }
}

#' Helper function to find inverse of a function
#' Used to find test statistics that correspond to mirrored p-values
#' From Mike Axiak, https://stackoverflow.com/questions/10081479/solving-for-the-inverse-of-a-function-in-r
#' @param f function
#'
#' @noRd
.inverse = function (f, lower = 0, upper = 200) {
  return(function (y) uniroot((function (x) f(x) - y), lower = lower, upper = upper,tol=.Machine$double.eps^2)[1])
}

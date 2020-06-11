#' Create args class
#'
#' Class containing all relevant arguments for model
#'
#' @param testing "\code{one_sided}" or "\code{interval}".
#' @param rendpoint Right endpoint of interval null
#' @param lendpoint Left endpoint of interval null
#' @param alpha_m The maximum possible rejected p-value.
#' @param zeta Controls minimum possible number of rejections.
#' @param lambda Controls where p-values are mirrored.
#' @param masking_shape "\code{tent}" or "\code{comb}".
#' @param niter Number of iterations in EM procedure
#' @param ndf Degrees of freedom of spline basis, number of dimensions of \eqn{\beta}.
#' @param nclasses Number of classes in Gaussian Mixture Model, minimum 2.
#' @return args class
#' @noRd
construct_args <- function(testing,rendpoint,lendpoint,alpha_m,zeta,lambda,masking_shape,niter,n,ndf=NULL,nclasses=NULL){

  all_a <- c("s","b")
  if(testing=="one_sided"){
    z_to_p <- function(z) pnorm(z,lower.tail = FALSE)
    p_to_z <- function(p) -qnorm(p)
    jacobian <- prob_jacobian_one_sided
  }else if(testing=="interval"){
    if(is.null(lendpoint)){
      lendpoint <- -1 * rendpoint
    }
    radius <-  (rendpoint-lendpoint)/2
    z_to_p <- function(z) pnorm(abs(z)+radius,lower.tail=FALSE)+pnorm(-abs(z)+radius)
    p_to_z_inv <- inverse(z_to_p,lower = 0)
    p_to_z <- function(z) unlist(mapply(p_to_z_inv,z))
    all_a <- c(all_a,"s_neg","b_neg")

    jacobian <- function(z,mean,var)prob_jacobian_interval(z,mean,var,radius=radius)
  }else{
    stop("Invalid testing type inputted. Valid forms of testing: `one_sided` and `interval`.")
  }

  args <- list(testing = testing,
               rendpoint = rendpoint,
               lendpoint = lendpoint,
               alpha_m = alpha_m,
               zeta = zeta,
               lambda = lambda,
               niter = niter,
               masking_shape = masking_shape,
               p_to_z = p_to_z,
               z_to_p = z_to_p,
               ndf = ndf,
               nclasses = nclasses,
               all_a = all_a,
               n  = n,
               jacobian = jacobian
               )
  class(args) <- "args"
  return(args)
}

#' Create data class
#'
#' Class containing all relevant data inputs from user
#' @param x Vector of covariates
#' @param z Vector of test statistics
#' @param p_values Vector of p-values
#' @param args args class containing masking function arguments
#'
#' @return data class
#' @noRd
construct_data <- function(x,p_values,z,args){
  if(args$testing == "interval"){
    center <- (args$rendpoint + args$lendpoint)/2
    z <- abs(z-center)
  }
  data <- list(x = x,
               p_values = p_values,
               z = z
               )
  class(data) <- "data"

  data <- data_preprocessing(data,args)

  return(data)
}




#' Initialize Parameters for Gaussian Mixture Model
#'
#' @param args args class
#'
#' @details Randomly samples the entries of the beta matrix, and sets the first entry to 2
#' to ensure high probability of class 0. The mu parameter is randomly drawn from -6 to 6 for interval testing,
#' and randomly drawn from 1 to 6 for one sided testing.
#' T
#' @return params class containing beta, mu, tau, var (tau^2+1)
#' @noRd
initialize_params <- function(args){

  nclasses <- args$nclasses

  beta <- NULL
  if(args$testing == "interval"){
    mu <-  c(0, runif(nclasses - 1, min = -6, max = 6))
  }else{
    mu <-  c(0, runif(nclasses - 1, min = 1, max = 6))
  }

  tau <-   c(0, runif(nclasses - 1, min = 0.1, max = 3))
  var <- tau^2 + 1

  params <- list(beta=beta, mu=mu, var=var)
  class(params) <- "params"
  return(params)
}

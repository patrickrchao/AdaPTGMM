#' Perform Expectation Maximization Algorithm on Gaussian Mixture Model
#'
#' @description Will perform \code{niter} iterations of the expectation maximization procedure.
#' This procedure estimates \code{gammas}, the probabilities for each class, then uses these to fit \eqn{\beta}.
#' With the updated beta, data, and parameters, we estimate the probability of a class and true p-value based on
#' the covariate and masked value.
#'
#' @param model Model class with data, args, and initialize parameters
#'
#' @return Model class with updated parameters.
EM <- function(model){
  data <- model$data
  args <- model$args

  for(i in seq(args$niter)){
    w_ika <- e_step_w_ika(model)
    gammas <- e_step_gamma(model,w_ika)
    model <- m_step_beta(model,gammas)
    w_ika <- e_step_w_ika(model)
    model$params <- m_step_mu_tau(model,w_ika)
  }

  return(model)
}

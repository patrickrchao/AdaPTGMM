#' Perform Expectation Maximization Algorithm on Gaussian Mixture Model
#'
#' @description Will perform \code{niter_fit} iterations of the expectation maximization procedure.
#' This procedure estimates \code{gammas}, the probabilities for each class, then uses these to fit \eqn{\beta}.
#' With the updated beta, data, and parameters, we estimate the probability of a class and true p-value based on
#' the covariate and masked value.
#'
#' @param model Model class with data, args, and initialize parameters
#'
#' @return Model class with updated parameters.
EM <- function(model, w_ika=NULL,preset_iter=NULL,return_w_ika = FALSE){


  if(is.null(preset_iter)){
    niter  <- model$args$niter_fit
  }else{
    niter <- preset_iter
  }

  w_ika <- NULL

  for(i in seq(niter)){
    w_ika <- e_step_w_ika(model, w_ika)
    gammas <- e_step_gamma(w_ika)
    # Do not update beta for the first iteration if in model selection
    # This is to initialize the model with an intercept only model
    if(niter == 1 | i > 1 | !is.null(model$params$beta)){
      model <- m_step_beta(model,gammas)
    }
    model$params  <- m_step_mu_tau(model,w_ika)
  }

  if(return_w_ika){
    output <- list(model=model,w_ika=w_ika)
  }else{
    output <- model
  }


  return(output)
}


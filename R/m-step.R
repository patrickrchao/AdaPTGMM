#' Perform Maximization step for beta
#'
#' @param model model variable containing data, args, params
#' @param gammas Estimated gammas from e_step_gammas use as weights, dataframe where columns correspond to classes
#' and rows correspond to hypotheses. Each row sums to 1.
#'
#' @details Performs a multinomial logistic regression
#' @return model
#' @noRd
m_step_beta <- function(model,gammas){

  model_type <- model$args$model_type
  nclasses <- model$args$nclasses
  beta_formula <- model$args$beta_formula

  x <- model$data$x
  rownames(x) <- NULL

  out <- m_step_beta_defaults(model_type,beta_formula,x,gammas,model$params$beta)

  model$params$beta <- out$model_weights
  model$params$df <- out$df
  prob <-  out$fitted_prob
  if(nclasses == 2 & ncol(prob) == 1){
    prob <- cbind(1-prob,prob)
  }
  if(model$args$n != nrow(prob)){

    prob <- prob[1:model$args$n,]
  }
  prob <- as.data.frame(pmax(pmin(as.matrix(prob),1 - 1e-12), 1e-12))
  model$data$class_prob <- prob
  return(model)
}

#' Perform Maximization step for mu and tau
#' @noRd
m_step_mu_tau <- function(model,w_ika){
  args <- model$args
  params <- model$params
  data <- model$data
  z <- w_ika$z

  for (k in 1:args$nclasses){
    subset <- w_ika[w_ika$class == k,]
    se <- data$se[subset$i]
    if(length(unique(se)) == 1){
      params$mu[k] <- weighted.mean(subset$z, subset$value)
      params$var[k] <- max(weighted.mean((subset$z - params$mu[k])^2, subset$value) - se[1]^2, 0)
    }else{

      params$mu[k] <- weighted.mean(subset$z,subset$value/(params$var[k]+se^2))

      for(iter in 1:5){
        grad <- - sum(subset$value / (params$var[k] + se^2)) +
          sum(subset$value*(subset$z-params$mu[k])^2 / (params$var[k]+se^2)^2)
        second_derivative <- sum(subset$value / (params$var[k]+se^2)^2) -
          2*sum(subset$value * (subset$z-params$mu[k])^2 / (params$var[k]+se^2)^3)

        params$var[k] <- max(params$var[k] - 0.5*grad/second_derivative, 0)
      }
    }


  }
  if(any(is.na(params$mu)) | any(is.na(params$var))){

    stop("NA value for mu or variance found. Stopping.")
  }
  return(params)
}


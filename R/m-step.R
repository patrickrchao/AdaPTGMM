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

  old_params <- params

  for (k in 1:args$nclasses){
    subset <- w_ika[w_ika$class == k,]
    se <- data$se[subset$i]

    # Optim for mu and tau^2

    computed_optim <- optim(c(params$mu[k],params$var[k]),function(x){
      mu <- x[1]
      var <- x[2]

      mean(subset$value*(log(se^2+var) + (subset$z-mu)^2/(var+se^2)))
    },gr <- function(x){
      mu <- x[1]
      var <- x[2]
      c(-mean(subset$value*(subset$z-mu)/(var+se^2)),
        mean(subset$value / (var + se^2)) -
          mean(subset$value*(subset$z-mu)^2 / (var+se^2)^2)
      )
    },method="L-BFGS-B",control=list(maxit=10,factr=1e-2),lower=c(-Inf,0),upper=c(Inf,Inf))
    params$mu[k] <- computed_optim$par[1]
    params$var[k] <- max(computed_optim$par[2],0)

  }



  if(any(is.na(params$mu)) | any(is.na(params$var))){

    stop("NA value for mu or variance found. Stopping.")
  }

  return(params)
}


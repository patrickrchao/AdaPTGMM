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
  if(nclasses == 2){
    prob <- cbind(1-prob,prob)
  }
  if(model$args$n != nrow(prob)){

    prob <- prob[1:model$args$n,]
  }
  prob <- pmax(pmin(prob,1 - 1e-12), 1e-12)
  model$data$class_prob <- prob
  return(model)
  # multinom_data <- data.frame(x,gammas)
  # multinom_data$class <- multinom_data$class - 1
  #
  # if(model_type == "mgcv"){
  #   multinom_data <- multinom_data[multinom_data$class == 1,]
  #   multinom_data$class <- multinom_data$value
  # }
  # #if(model_type == "glm"){
  # #}else if(model_type == "gam"){
  #  #formulas <- lapply(c(paste("class",beta_formula),rep(beta_formula,nclasses-2)),as.formula)
  # #}
  # # if the multinom beta model exists, use the previous weights as the starting point for faster convergence
  # if(!is.null(model$params$beta)){
  #   if(model_type == "glm"){
  #     est_beta <- nnet::multinom(beta_formula, multinom_data, weights = value, trace = F,maxit=10,reltol=1e-7,Wts=model$params$beta)
  #   }else if(model_type == "gam"){
  #     est_beta <- VGAM::vgam(beta_formula,multinomial,multinom_data,weights = value,coefstart=model$params$beta,control=vgam.control(maxit=5,bf.maxit = 5,trace=F))
  #   }else if(model_type == "mgcv"){
  #    # est_beta <- mgcv::gam(beta_formula,quasibinomial,multinom_data,weights = value,control=vgam.control(maxit=5,bf.maxit = 5,trace=F))
  #     est_beta <- mgcv::gam(beta_formula, family=quasibinomial, data=multinom_data,
  #                     control=gam.control(maxit=1),outer=gam.outer(start=model$params$beta))
  #     }
  # }else{
  #   if(model_type == "glm"){
  #     est_beta <- nnet::multinom(beta_formula, multinom_data, weights = value, trace = F,maxit=100,reltol=1e-8)
  #   }else if(model_type == "gam"){
  #     est_beta <- vgam(beta_formula,multinomial,multinom_data,weights = value,control=vgam.control(maxit=15,bf.maxit = 15,trace=F))
  #   }else if(model_type == "mgcv"){
  #     est_beta <- mgcv::gam(beta_formula, family=quasibinomial, data=multinom_data)
  #   }
  # }
  #
  # model$data$class_prob <- class_prob(est_beta,nclasses,model$args$n,model_type)
  #
  # if(model_type == "glm"){
  #   model$params$beta <- est_beta$wts
  #   model$params$df <- sum(est_beta$edf)
  # }else if(model_type == "gam"){
  #   model$params$beta <- coef(est_beta)
  #   model$params$df <- nobs(est_beta,type="vlm") -df.residual(est_beta)
  # }else if(model_type == "mgcv"){
  #   model$params$df <- sum(est_beta$edf)
  #   model$params$beta <- est_beta$coefficients
  # }


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
  #for (k in 1:(args$nclasses-1)){
    subset <- w_ika[w_ika$class == k,]
    params$mu[k] <- .weighted_mean(subset$z,subset$value)
    # Minimum variance for convolved Gaussian is 1
    params$var[k] <- max(.weighted_mean((subset$z-params$mu[k])^2,subset$value), 1)
  }
  return(params)
}


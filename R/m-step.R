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
  multinom_data <- data.frame(x,gammas)
  multinom_data$class <- multinom_data$class - 1

  if(model_type == "glm"){
    formula <- as.formula(paste("class",beta_formula))
  }else if(model_type == "gam"){
   formulas <- lapply(c(paste("class",beta_formula),rep(beta_formula,nclasses-2)),as.formula)
  }

  # if the multinom beta model exists, use the previous weights as the starting point for faster convergence
  if(!is.null(model$params$beta)){
    if(model_type == "glm"){
      est_beta <- nnet::multinom(formula, multinom_data, weights = value, trace = FALSE,maxit=5,Wts=model$params$beta,optimizer = c("outer","bfgs"))
    }else if(model_type == "gam"){

     # if(nclasses ==2 ){
    #    est_beta <- bam(formulas, family=quasibinomial(), data=multinom_data,
   #                    weights = value,control=gam.control(maxit=1),coef=model$params$beta)
     # }else{
        est_beta <- gam(formulas, family=mgcv::multinom(K=nclasses-1), data=multinom_data,optimizer = c("efs"),
                        weights = value,control=gam.control(maxit=1),outer=gam.outer(start=model$params$beta))
      #}
    }
  }else{
    if(model_type == "glm"){
      est_beta <- nnet::multinom(formula, multinom_data, weights = value, trace = FALSE,maxit=100)
    }else if(model_type == "gam"){
      #if(nclasses ==2 ){
      #  est_beta <- bam(formulas, family=quasibinomial(), data=multinom_data,
      #                  weights = value,control=gam.control(maxit=3,trace=T))
      #}else{

        est_beta <- gam(formulas, family=mgcv::multinom(K=nclasses-1), data=multinom_data,weights = value,
                        optimizer = c("efs"),control=gam.control(maxit=3))
      #}

    }
  }
  model$params$df <- sum(est_beta$edf)
  model$data$class_prob <- class_prob(est_beta,nclasses,model$args$n,model_type)
  #model$params$beta <- est_beta$wts
  if(model_type == "glm"){
    model$params$beta <- est_beta$wts
  }else if(model_type == "gam"){
    model$params$beta <- est_beta$coefficients
  }


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


m_step_beta_defaults <- function(model_type,formula, x,gammas, model_weights,m_step_custom){
  data <- data.frame(x,gammas)
  colnames(data) <- c(colnames(x),"class","weights")
  data$class <- data$class - 1
  data$class <- as.factor(data$class)

  if(model_type == "nnet"){
    if (!requireNamespace("nnet", quietly = TRUE)){
      stop("package \'nnet\' not found. Please install.")
    }
    out <- m_step_nnet(formula,data,model_weights)
  }else if(model_type == "rrvglm"){
    if (!requireNamespace("VGAM", quietly = TRUE)){
      stop("package \'VGAM\' not found. Please install.")
    }
    out <- m_step_rrvglm(formula,data,model_weights)
  }else if(model_type == "mgcv"){
    if (!requireNamespace("mgcv", quietly = TRUE)){
      stop("package \'mgcv\' not found. Please install.")
    }
    out <- m_step_mgcv(formula,data,model_weights)
  }else if(model_type == "glmnet"){
    if (!requireNamespace("glmnet", quietly = TRUE)){
      stop("package \'glmnet\' not found. Please install.")
    }
    out <- m_step_glmnet(formula,data,model_weights)
  }else if(model_type == "neural"){
    if (!requireNamespace("nnet", quietly = TRUE)){
      stop("package \'nnet\' not found. Please install.")
    }
    out <- m_step_neural(formula,data,model_weights)
  }else if(model_type == "custom"){
    out <- m_step_custom(formula,data,model_weights)
  }else{
    stop("Invalid beta fitting method found.")
  }
  return(out)
}

m_step_nnet <- function(formula, data, model_weights){

  if(is.null(model_weights)){
    est_beta <- nnet::multinom(formula=formula, data=data, weights = weights, trace = F)
  }else{
    est_beta <- nnet::multinom(formula, data, weights = weights, trace = F,Wts=model_weights)
  }
  fitted_prob <- fitted(est_beta)
  new_model_weights <- est_beta$wts
  df <- sum(est_beta$edf)
  return(list(fitted_prob=fitted_prob,
              model_weights = new_model_weights,
              df = df))

}



m_step_neural <- function(formula, data, model_weights){

  if(is.null(model_weights)){
    est_beta <- nnet::nnet(formula=formula, data=data, weights = weights, size=3,trace=F,reltol=1.0e-4)
  }else{
    est_beta <- nnet::nnet(formula=formula, data=data, weights = weights, Wts = model_weights, size=3,trace=F,reltol=1.0e-4)
  }
  fitted_prob <- fitted(est_beta)
  new_model_weights <- est_beta$wts
  df <- length(new_model_weights)
  return(list(fitted_prob=fitted_prob,
              model_weights = new_model_weights,
              df = df))

}


m_step_glmnet <- function(formula, data,model_weights){

  # Does not initialize model_weights
  est_beta <- glmnetUtils::cv.glmnet(formula=formula,data=data,weights=data$weights, family="multinomial", nfolds=3,maxit=1e5,nlambda=5)

  # Code from https://github.com/lihualei71/adaptMT/blob/master/R/safe-model.R

  fitted_prob <- drop(predict(est_beta, formula=formula,newdata=data, s="lambda.min",type="response"))
  new_coef <- coef(est_beta,s="lambda.min")
  vi <- lapply(new_coef,function(x){as.numeric(x!= 0)})
  df <- sum(sapply(vi,sum)) + 1
  return(list(fitted_prob=fitted_prob,
              model_weights = new_coef,
              df = df))

}

m_step_rrvglm <- function(formula, data,model_weights){

  data$weights  <- pmax(pmin(data$weights,1-1e-12), 1e-12)
  est_beta <- suppressWarnings(
    VGAM::rrvglm(formula,multinomial,data,weights = weights,Rank=1,control=rrvglm.control(algorithm="derivative",trace=F))
    )

  fitted_prob <- predict(est_beta,type="response")
  new_model_weights <- coef(est_beta)

  df <- nobs(est_beta,type="vlm") - df.residual(est_beta)
  return(list(fitted_prob=fitted_prob,
              model_weights = new_model_weights,
              df = df))

}

m_step_mgcv <- function(formula, data,model_weights){

  data <- data[data$class == 1,]
  data$class <- pmax(pmin(data$weights,1-1e-12), 1e-12)

  if(is.null(model_weights)){
    est_beta <- mgcv::gam(formula, family=quasibinomial, data=data)
  }else{
    est_beta <- mgcv::gam(formula, family=quasibinomial, data=data,
                          outer=gam.outer(start=model_weights))
  }
  fitted_prob <- data.frame(fitted(est_beta))
  new_model_weights <-  est_beta$coefficients
  df <- sum(est_beta$edf)
  return(list(fitted_prob=fitted_prob,
              model_weights = new_model_weights,
              df = df))
}


custom_beta <- function(formula, data, initialization){

  ########################################################################
  # EDIT THIS SECTION to fit any desired model
  # Formula, data, and weights for each observation in data$weights
  # Compute:
  # 1. The fitted probabilities for the data
  # 2. New model parameters for initialization (or leave as null if undesired)
  # 3. Degrees of freedom for model (for model selection)

  model <- fit_model(formula,data,params=initialization,weights=data$weights)

  fitted_probabilities <- fitted(model)
  new_weights <- model$weights
  df <- model$df

  ########################################################################
  output <- list()
  output$fitted_prob <- fitted_probabilities
  output$model_weights <- new_weights
  output$df <- df

  return(output)

}

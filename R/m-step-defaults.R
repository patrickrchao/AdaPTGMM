m_step_beta_defaults <- function(model_type,formula, x,gammas, model_weights){

  data <- data.frame(x,gammas)
  colnames(data) <- c(colnames(x),"class","weights")
  data$class <- data$class - 1

  if(model_type == "glm"){
    out <- m_step_glm(formula,data,model_weights)
  }else if(model_type == "gam"){
    out <- m_step_gam(formula,data,model_weights)
  }else if(model_type == "mgcv"){
    out <- m_step_mgcv(formula,data,model_weights)
  }else if(model_type == "glmnet"){
    out <- m_step_glmnet(formula,data,model_weights)
  }
  return(out)
}

m_step_glm <- function(formula, data, model_weights){
  if(is.null(model_weights)){
    est_beta <- nnet::multinom(formula, data, weights = weights, trace = F,maxit=100,reltol=1e-7)
  }else{
    est_beta <- nnet::multinom(formula, data, weights = weights, trace = F,maxit=10,reltol=1e-7,Wts=model_weights)

  }

  fitted_prob <- fitted(est_beta)
  new_model_weights <- est_beta$wts
  df <- sum(est_beta$edf)
  return(list(fitted_prob=fitted_prob,
              model_weights = new_model_weights,
              df = df))

}

m_step_glmnet <- function(formula, data,model_weights){

 # if(is.null(model_weights)){
    est_beta <- glmnetUtils::cv.glmnet(formula=formula,data=data,weights=data$weights, family="multinomial", nfolds=3,maxit=1e3,nlambda=5)
    #est_beta <- glm(formula=formula, data=data, weights = weights,family=multinomial())
  #}else{
   # est_beta <- glmnetUtils::cv.glmnet(formula=formula,data=data,weights=data$weights, family="multinomial")
    #est_beta <- glm(formula, data=data, weights = weights,family="multinomial",start=model_weights)
  #}

  # Code from https://github.com/lihualei71/adaptMT/blob/master/R/safe-model.R

  fitted_prob <- drop(predict(est_beta, formula=formula,newdata=data, s="lambda.min",type="response"))
  new_coef <- coef(est_beta,s="lambda.min")
  vi <- lapply(new_coef,function(x){as.numeric(x!= 0)})
  df <- sum(sapply(vi,sum)) + 1
  return(list(fitted_prob=fitted_prob,
              model_weights = new_coef,
              df = df))

}

m_step_gam <- function(formula, data,model_weights){
  if(is.null(model_weights)){
    est_beta <- vgam(formula,multinomial,data,weights = weights,control=vgam.control(maxit=15,bf.maxit = 15,trace=F))
  }else{
    est_beta <- VGAM::vgam(formula,multinomial,data,weights = weights,coefstart=model_weights,control=vgam.control(maxit=5,bf.maxit = 5,trace=F))
  }

  fitted_prob <- predict(est_beta,type="response")
  new_model_weights <- coef(est_beta)
  df <- nobs(est_beta,type="vlm") -df.residual(est_beta)
  return(list(fitted_prob=fitted_prob,
              model_weights = new_model_weights,
              df = df))

}

m_step_mgcv <- function(formula, data,model_weights){
  return(0)

}



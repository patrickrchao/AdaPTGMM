m_step_beta_defaults <- function(model_type,formula, x,gammas, model_weights){
  if(model_type == "dirich"){
    return(m_step_dirich(beta_formula = formula,x = x,gammas=gammas,model_weights = model_weights))
  }
  data <- data.frame(x,gammas)
  colnames(data) <- c(colnames(x),"class","weights")
  data$class <- data$class - 1

  if(model_type == "nnet"){
    out <- m_step_nnet(formula,data,model_weights)
  }else if(model_type == "gam"){
    out <- m_step_gam(formula,data,model_weights)
  }else if(model_type == "mgcv"){
    out <- m_step_mgcv(formula,data,model_weights)
  }else if(model_type == "glmnet"){
    out <- m_step_glmnet(formula,data,model_weights)
  }else{
    warning("Invalid beta fitting method found.")
  }
  return(out)
}

m_step_nnet <- function(formula, data, model_weights){

  if(is.null(model_weights)){
    est_beta <- nnet::multinom(formula=formula, data=data, weights = weights, trace = F)
  }else{
    est_beta <- nnet::multinom(formula, data, weights = weights, trace = F,Wts=model_weights)

  }
  # plot(test_data$x,filter(test_data$weights,rep(1/1000,1000),sides=2),ylim=c(0,1))
  # abline(0,1,lwd=3,col="black")
  # plot(seq(0,1,length.out = 100),predict(est_beta,data.frame(x=seq(0,1,length.out = 100)),type="probs"))
  #temp <- data[data$x>0.95,]
 # print(mean(temp[temp$class==1,]$weights))

  #print(predict(est_beta,data.frame(x=seq(0.9,1,length.out = 10)),type="probs"))
  fitted_prob <- fitted(est_beta)
  new_model_weights <- est_beta$wts
  df <- sum(est_beta$edf)
  return(list(fitted_prob=fitted_prob,
              model_weights = new_model_weights,
              df = df))

}

m_step_glmnet <- function(formula, data,model_weights){


 # if(is.null(model_weights)){
    est_beta <- glmnetUtils::cv.glmnet(formula=formula,data=data,weights=data$weights, family="multinomial", nfolds=3,maxit=1e5,nlambda=5)
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

  data <- data[data$class == 1,]
  data$class <- pmin(pmax(data$weights,1-1e-12),1e-12)

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


m_step_dirich <- function(beta_formula, x,gammas,model_weights){

  nclasses <- round(nrow(gammas)/nrow(x))
  data <- data.frame(matrix(gammas$value,ncol = nclasses))
  data$x <- x$x
  dr_data <- DR_data(data[,1:nclasses])
  testing_formula <- as.Formula("dr_data ~ splines::ns(x,df=3)")
  est_beta <- DirichReg(dr_data ~ splines::ns(x,df=3),data=data,control=list(tol1=1e-3,tol2=1e-5,iterlim=100))

  fitted_prob <- predict(est_beta,type="response")
  new_model_weights <- est_beta$coefficients
  df <- est_beta$npar
  return(list(fitted_prob=fitted_prob,
              model_weights = new_model_weights,
              df = df))

}



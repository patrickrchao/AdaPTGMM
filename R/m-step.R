library(glmnet)
# Fits the logistic glm beta based on the gammas as labels and x as inputs
fit_beta <- function(x,gammas,est_params){
  #if(sum(gammas>1)>0 || sum(gammas<0)>0||  sum(is.na(gammas))>0){
  #browser()
  #}
  logistic_glm_data <- data.frame(x,gammas)
  beta_names <- coefficient_names("Beta",ncol(x))
  colnames(logistic_glm_data) <- c(beta_names,"gammas")
  # -1 to ignore intercept

  formula <- paste("gammas ~ ",paste0(colnames(logistic_glm_data)[c(1:ncol(x))],collapse=" + ")," -1",sep="")
  #glm_beta <- suppressWarnings(glm(formula,data=logistic_glm_data,family=binomial()))

  glm_beta <- suppressWarnings(glmnet(x=data.matrix(x[,-1]),y=matrix(c(1-gammas,gammas),ncol=2),family="binomial",lambda=0.00001,intercept=TRUE))
  #beta <- glm_beta$coefficients

  beta <- c(glm_beta$a0,as.vector(glm_beta$beta))
  return(beta)
}

update_parameters <- function(data,est_params,params){
  gammas <- expectation_gamma(data,est_params,params)
  curr_beta <- fit_beta(data$full_x,gammas,est_params)
  w_ia_full <- calculate_w(data,est_params,params)
  w_ia <- w_ia_full %>% filter(k==1)
  new_mu <- weighted_mean(w_ia$w_ika,w_ia$z)
  est_params$mu[2] <- new_mu
 # new_var <- weighted_mean(w_ia$w_ika,(w_ia$z-est_params$mu)^2)
  #print(paste(round(est_params$var[2],5),round(new_var,5)))
  #est_params$var[2] <- new_var

  if(params$testing_interval){

    for(iter in seq(3)){

      #w_ia_full <- calculate_w(data,est_params,params)
      #w_ia <- w_ia_full %>% filter(k==1)
      gradients <- calculate_gradients(data, est_params,params,w_ia)
      #print(paste0("Iter: ",iter," Update: ",round(gradients$var,5)))
      est_params$var[2] <- max(est_params$var - 0.1/sqrt(iter)* pmin(pmax(gradients$var,-2),2),0.1)
      #print(paste(est_params$mu,est_params$var))

    }
  }else{
    est_params$var[2] <- weighted_mean(w_ia$w_ika,(w_ia$z-est_params$mu)^2)

  }
  #print(paste(round(est_params$var[2],5),round(new_var,5)))
  est_params$beta <- curr_beta
  #likelihood(data,est_params,params,w_ika=w_ia_full)
  return(est_params)
}






calculate_gradients <- function(data,est_params,params,w_ia){
  gradients <- list()
  mu <- est_params$mu[2]
  var <- est_params$var[2]
  z <- w_ia$z
  #browser()
  gradients$var <- sum(w_ia$w_ika*(1/var - (z-mu)^2/(var^2)))
  gradients$var_second <- sum(w_ia$w_ika*(-1/(var^2) +2*(z-mu)^2/(var^3)))
  #print(paste("gradient",gradients$var,"second",gradients$var_second))
  return(gradients)
}


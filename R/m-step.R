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
  glm_beta <- suppressWarnings(glm(formula,data=logistic_glm_data,family=binomial()))
  beta <- glm_beta$coefficients
  return(beta)
}

update_parameters <- function(data,est_params,params){
  gammas <- expectation_gamma(data,est_params,params)
  curr_beta <- fit_beta(data$full_x,gammas,est_params)
  w_ia <- calculate_w(data,est_params,params)
  est_params$mu[2] <- weighted_mean(w_ia$w_ia,w_ia$z)
  if(params$testing_interval){

    for(iter in seq(5)){


      gradients <- calculate_gradients(data, est_params,params,w_ia)
      #print(paste0("Iter: ",iter," Gradient: ",round(gradients$var/gradients$var_second,5)))
      est_params$var[2] <- max(est_params$var - 0.1 * gradients$var/gradients$var_second,0.1)
      #print(paste(est_params$mu,est_params$var))

    }
  }else{
    est_params$var[2] <- weighted_mean(w_ia$w_ia,(w_ia$z-est_params$mu)^2)

  }
  est_params$beta <- curr_beta
  return(est_params)
}






calculate_gradients <- function(data,est_params,params,w_ia){
  gradients <- list()
  mu <- est_params$mu[2]
  var <- est_params$var[2]
  z <- w_ia$z
  gradients$var <- mean(w_ia$w_ia*(1/var - (z-mu)^2/var^2))
  gradients$var_second <- mean(w_ia$w_ia*(-1/var^2 +2*(z-mu)^2/var^3))
  #gradients$var <- mean(w_ia$w_ia/(2*var^2)*
   #                       (mu^2-2*mu*w_ia$w_ia*coth_term-var+z^2))
  return(gradients)
}


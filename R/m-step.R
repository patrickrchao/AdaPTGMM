library(glmnet)
library(nnet)
# Fits the logistic glm beta based on the gammas as labels and x as inputs
fit_beta <- function(x,gammas,params){

  if(sum(gammas>1)>0 || sum(gammas<0)>0||  sum(is.na(gammas))>0 || sum(is.na(x))>0){
  browser()
  }
  #print(colSums(gammas))

  num_classes <- ncol(gammas)
  colnames(gammas) <- 0:(ncol(gammas)-1)
  spread_gammas <- gather(gammas,class,weight)
  gamma_class <- spread_gammas$class
  gamma_weight <- spread_gammas$weight
  logistic_glm_data <- data.frame(x,gamma_class,gamma_weight)
  beta_names <- coefficient_names("Beta",ncol(x))
  colnames(logistic_glm_data) <- c(beta_names,"gammas","weight")
  formula <- paste("gammas ~ ",paste0(colnames(logistic_glm_data)[c(1:ncol(x))],collapse=" + ")," -1",sep="")
  beta_model <- multinom(formula,logistic_glm_data,weight = weight,trace=FALSE)
  beta <- as.matrix(coef(beta_model))

  if (num_classes > 2){
    beta <- t(beta)
  }

    #print(dim(params$beta))
  #print(paste("Correct",num_df,num_classes-1))
  #}

 # browser()
 # glmnet.control(pmin=0)
 #
 #
 #  if(ncol(x)==1){
 #    beta <- matrix(log(colMeans(gammas)),ncol=ncol(gammas))
 #    colnames(beta, do.NULL = FALSE)
 #    colnames(beta) <- 0:(ncol(gammas)-1)
 #
 #  }else{
 #    # empty_cols <- colSums(gammas) == 0
 #    # gammas <- gammas[,!empty_cols]
 #    # glm_beta <- suppressWarnings(glmnet(x=scale(data.matrix(x),T,F),
 #    #                                     y=data.matrix(gammas),family="multinomial",lambda=0,alpha=0,intercept=FALSE,
 #    #                                     penalty.factor = c(0, rep(1, ncol(x)-1))
 #    #                                     ))
 #    #
 #    # #beta <- as.matrix(as.data.frame(lapply(coef(glm_beta),function(x) as.matrix(x)[,"s1"])))#as.matrix(as.data.frame(lapply(coef(glm_beta),as.matrix)))
 #    # beta <- as.matrix(as.data.frame(lapply(coef(glm_beta),as.matrix)))
 #    # row.names(beta)[1] <- "Beta0"
 #    # #colnames(beta) <- names(empty_cols)[!empty_cols]
 #    # zero_mat <- matrix(0, nrow=nrow(beta), ncol=sum(empty_cols))
 #    # colnames(zero_mat, do.NULL = FALSE)
 #    # colnames(zero_mat) <- names(empty_cols)[empty_cols]
 #    # beta <- cbind(beta,zero_mat)
 #    # beta <- beta[,order(colnames(beta))]
 #    # colnames(beta) <- 0:(ncol(beta)-1)
 #
 #   #  empty_cols <- colSums(gammas) == 0
 #   #  gammas <- gammas[,!empty_cols]
 #   # glm_beta <- suppressWarnings(glmnet(scale(x=data.matrix(x[,-1]),T,F),y=data.matrix(gammas),family="multinomial",lambda=0,alpha=0,intercept=TRUE))
 #   #
 #   # #temp_beta <- as.matrix(as.data.frame(lapply(coef(glm_beta),function(x) as.matrix(x)[,"s1"])))#as.matrix(as.data.frame(lapply(coef(glm_beta),as.matrix)))
 #   # temp_beta <- as.matrix(as.data.frame(lapply(coef(glm_beta),as.matrix)))
 #   # row.names(temp_beta)[1] <- "Beta0"
 #   # #colnames(beta) <- names(empty_cols)[!empty_cols]
 #   # zero_mat <- matrix(0, nrow=nrow(temp_beta), ncol=sum(empty_cols))
 #   # colnames(zero_mat, do.NULL = FALSE)
 #   # colnames(zero_mat) <- names(empty_cols)[empty_cols]
 #   # temp_beta <- cbind(temp_beta,zero_mat)
 #   # temp_beta <- temp_beta[,order(colnames(temp_beta))]
 #   # colnames(temp_beta) <- 0:(ncol(temp_beta)-1)
 #   #
 #   # browser()
 #     empty_cols <- colSums(gammas) == 0
 #     gammas <- gammas[,!empty_cols]
 #    glm_beta <- suppressWarnings(glmnet(x=data.matrix(x[,-1]),y=data.matrix(gammas),family="multinomial",lambda=0.01,alpha=0,intercept=TRUE))
 #
 #    #beta <- as.matrix(as.data.frame(lapply(coef(glm_beta),function(x) as.matrix(x)[,"s1"])))#as.matrix(as.data.frame(lapply(coef(glm_beta),as.matrix)))
 #    beta <- as.matrix(as.data.frame(lapply(coef(glm_beta),as.matrix)))
 #    row.names(beta)[1] <- "Beta0"
 #    #colnames(beta) <- names(empty_cols)[!empty_cols]
 #    zero_mat <- matrix(0, nrow=nrow(beta), ncol=sum(empty_cols))
 #    colnames(zero_mat, do.NULL = FALSE)
 #    colnames(zero_mat) <- names(empty_cols)[empty_cols]
 #    beta <- cbind(beta,zero_mat)
 #    beta <- beta[,order(colnames(beta))]
 #    colnames(beta) <- 0:(ncol(beta)-1)
 #  }
 # #glm_beta <- suppressWarnings(glmnet(x=data.matrix(x[,-1]),y=data.matrix(gammas),family="binomial",lambda=0.00001,intercept=TRUE))
 # beta <- glm_beta$coefficients
 #
 # # beta <- c(glm_beta$a0,as.vector(glm_beta$beta))
 #
 #
 # #beta <- as.data.frame(lapply(coef(glm_beta),as.matrix))
 #
 # browser()
  # beta <- data.frame(t(beta))
  # beta$class <- rownames(beta)

  return(beta)
}

update_parameters <- function(data,params,args){

  gammas <- expectation_gamma(data,params,args)
  curr_beta <- fit_beta(data$full_x,gammas,params)
  output <- calculate_w(data,params,args)

  for (class in 1:(args$num_classes-1)){
    w_ia <- output$w_ika %>% filter(k==class)

    params$mu[class+1] <- weighted_mean(w_ia$w_ika,w_ia$z,w_ia)
    new_var <- max(weighted_mean(w_ia$w_ika,(w_ia$z-params$mu[class+1])^2,w_ia),0.0)
    params$var[class+1] <- new_var
  }

  #params$var[2] <- weighted_mean(w_ia$w_ika,(w_ia$z-params$mu)^2,w_ia)
  #gradients <- calculate_gradients(data, params,args,w_ia)

  # new_mu <- weighted_mean(w_ia$w_ika,w_ia$z,w_ia)
  #
  # params$var[2] <- weighted_mean(w_ia$w_ika,(w_ia$z-new_mu)^2,w_ia)
  # #params$var[2] <- weighted_mean(w_ia$w_ika,(w_ia$z-params$mu)^2,w_ia)
  # params$mu[2] <- new_mu
  # gradients <- calculate_gradients(data, params,args,w_ia)
  #print(paste("Mu grad",gradients$mu))
  #print(paste("Var grad",gradients$var))
  #print(paste(round(params$var[2],5),round(new_var,5)))
  #params$var[2] <- new_var

  # if(args$testing_interval){
  #
  #   for(iter in seq(3)){
  #
  #     #w_ia_full <- calculate_w(data,params,args)
  #     #w_ia <- w_ia_full %>% filter(k==1)
  #     gradients <- calculate_gradients(data, params,args,w_ia)
  #     #print(paste0("Iter: ",iter," Update: ",round(gradients$var,5)))
  #     params$var[2] <- max(params$var - 0.1/sqrt(iter)* pmin(pmax(gradients$var,-2),2),0.1)
  #     #print(paste(params$mu,params$var))
  #
  #   }
  # }else{
  #   params$var[2] <- weighted_mean(w_ia$w_ika,(w_ia$z-params$mu)^2)
  #
  # }
  #print(paste(round(params$var[2],5),round(new_var,5)))

  params$beta <- curr_beta
  #params <- sort_args(params)

  return(params)
}




#
#
# calculate_gradients <- function(data,params,args,w_ia){
#   gradients <- list()
#   mu <- params$mu[2]
#   var <- params$var[2]
#   z <- w_ia$z
#   #browser()
#   gradients$mu  <- sum(w_ia$w_ika * (z-mu)/var)
#   gradients$var <- sum(w_ia$w_ika * (1/var - (z-mu)^2/(var^2)))
#   gradients$var_second <- sum(w_ia$w_ika*(-1/(var^2) +2*(z-mu)^2/(var^3)))
#   #print(paste("gradient",gradients$var,"second",gradients$var_second))
#   return(gradients)
# }


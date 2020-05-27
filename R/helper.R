library(splines)
library(stringr)

expit <- function(x){
  return(exp(x)/(1+exp(x)))
}

weighted_mean <- function(weights,values,w_ia){
  total_weight <- sum(weights)
  if(sum(is.na(weights))>0){
    browser()
  }else if(total_weight==0 ){
    return(0)
  }
  if(sum(weights<0)>0){
    browser()
    stop("Warning! Invalid weights found")
    return(0)
  }
  return(sum(values*weights)/total_weight)
}

gaussian_pdf <- function(x,mean,var){
  pdf <- dnorm(x,mean,sqrt(var))
  #pdf[is.na(pdf)] <- 10^(-8)
  return(pdf)
}

generate_spline <- function(x,num_df){
  if(num_df > 1){
    spline_func <- ns(x,df=num_df-1)
    spline_x <- cbind(rep(1,length(x)),spline_func)
  }else{
    spline_x <- matrix(rep(1,length(x)))
    spline_func <- function(x){matrix(rep(1,length(x)))}
  }

  beta_names <- coefficient_names("Beta",ncol(spline_x))
  colnames(spline_x) <- beta_names
  return(list(spline_x=spline_x,spline_func= spline_func))
}

coefficient_names <- function(var_name,num_dim,start=0){
  x_names <- paste(var_name,(start+0):(start+num_dim-1),sep="")
  return(x_names)
}

inverse = function (f, lower = 0, upper = 200) {
  function (y) uniroot((function (x) f(x) - y), lower = lower, upper = upper,tol=.Machine$double.eps^0.5)[1]
}

change_of_density <- function(z, radius=1, mean, var) {
  return(dnorm(z,mean,sqrt(var))/
           (dnorm(-abs(z)+radius,0,1)-dnorm(abs(z)+radius,0,1)))
}

change_of_density_one_sided <- function(z, mean, var) {
  return(exp(z^2/2-(z-mean)^2/(2*var))/sqrt(var))
  #return(dnorm(z,mean,sqrt(var))/
   #        (dnorm(-abs(z)+radius,0,1)-dnorm(abs(z)+radius,0,1)))
}

probability_from_spline <- function(x, betas,num_classes){

  if(!(is.numeric(betas))| !(is.numeric(x))){
    browser()
  }
  if(ncol(betas) == num_classes ){

    # x is n by num_df
    # betas is num_df by (num_classes-1)
    class_prob <- exp(x %*% betas)
    class_prob <- class_prob/rowSums(class_prob)

  }else if(ncol(betas) == num_classes - 1){
    #Assuming using last column as reference column
    class_prob <- data.frame(matrix(ncol = num_classes, nrow = nrow(x)))
    exp_value <- exp(x %*% betas)
    row_sum <- rowSums(exp_value)

    class_prob[,1] <- 1/(1+row_sum)
    class_prob[,2:(num_classes)] <- exp_value/(1+row_sum)

  }else{
    browser()
    print("INVALID SOMETHING WENT WRONG")
  }

  return(class_prob)
}


# This calculates the probability of a_i and masked p_i given class
calculate_conditional_probabilities <- function(data,params,args){

  mask <- data$mask
  num_hypo <- length(data$small_z)
  mu <- params$mu
  var <- params$var
  small_z = data$small_z
  big_z = data$big_z
  z <- data$z
  radius <- args$interval_radius

  beta <- params$beta

  #TODO CHANGE THIS TO PREDICT from glmnet

  class_prob <- probability_from_spline( data$full_x, beta, args$num_classes)
  colnames(class_prob) <- paste0("class_prob_",0:(args$num_classes-1))

  prob <- data.frame(matrix(NA, nrow=length(small_z), ncol=1))
  if(args$testing_interval){
    # Iterate over s,b,-s,-b
    for(class in 0:(args$num_classes-1)){
      for(a in args$all_a){

        name <- paste0(a,"_",class)
        sign <- 1
        if(var[class+1] <- 0){
          prob[,name] <- 0
        }else{
          if(str_detect(a,"neg")){
            sign <- -1
          }
          if(str_detect(a,"s")){
            prob[mask,name] <- change_of_density(sign*small_z[mask], radius, mu[class+1], var[class+1])
            if(!str_detect(a,"neg")){

              prob[!mask,name] <- change_of_density(z[!mask], radius, mu[class+1], var[class+1])
            }else{
              prob[!mask,name] <- 0
            }
          }else{
            prob[mask,name] <- change_of_density(sign*big_z[mask], radius, mu[class+1], var[class+1])/args$zeta
            prob[!mask,name] <- 0
          }
        }
      }
    }

    prob[!mask ,"neg_b_0"] <- 0
  } else {
    prob[,"s_0"] <- 1
    prob[,"b_0"] <- mask * 1 / args$zeta
    for(class in 1:(args$num_classes-1)){
      for(a in args$all_a){
        name <- paste0(a,"_",class)
        if(var[class+1]==0){
          prob[,name] <- 0
        }else{

          if(str_detect(a,"s")){
            prob[mask,name] <- change_of_density_one_sided(small_z[mask], mu[class+1], var[class+1])
            prob[!mask,name] <- change_of_density_one_sided(z[!mask], mu[class+1], var[class+1])
            #prob[mask,name] <- gaussian_pdf(small_z[mask], mu[class+1], var[class+1]) / gaussian_pdf(small_z[mask], 0, 1)
            #prob[!mask,name] <- gaussian_pdf(z[!mask], mu[class+1], var[class+1]) / gaussian_pdf(z[!mask], 0, 1)
          }else{
            prob[mask,name] <- change_of_density_one_sided(big_z[mask], mu[class+1], var[class+1])
            #prob[mask,name] <- gaussian_pdf(big_z[mask], mu[class+1], var[class+1]) / gaussian_pdf(big_z[mask], 0, 1)/args$zeta
            prob[!mask,name] <- 0
          }
        }

      }
    }

  }

  prob <- prob[-c(1)]

  # for(class in 0:(args$num_classes-1)){
  #   class_columns <- str_detect(colnames(prob),toString(class))
  #   prob[,class_columns] <- prob[,class_columns] / rowSums(prob[,class_columns])
  # }
  prob <- cbind(prob,class_prob)


  return(prob)
}

log_likelihood <- function(model,optimal_param=FALSE,w_ika=FALSE,verbose=TRUE){
  data <- model$data
  params <- model$params
  args <- model$args
  full_x <-
  if(!is.numeric(w_ika)){
    w_ika <- calculate_w(data,params,args)

  }
  denominator <- w_ika$denominator
  log_likelihood <- sum(log(denominator))
  if(verbose){
    if(optimal_param){
      print(paste0("Log Likelihood with True Parameters: ",log_likelihood))
    }else{
      print(paste0("Log Calculated Likelihood: ",log_likelihood))
    }
  }
  return(log_likelihood)

}

initialize_params <- function(num_classes,num_df,interval=FALSE){
  #if(num_classes > 2 ){
  beta <-  matrix(sample(-4:4,(num_df*num_classes),replace=TRUE),ncol=num_classes)
  beta[1,] <- sample(seq(-1,1,0.5),num_classes,replace=TRUE)
  #beta[,1] <- 0#*(num_classes)
  beta[1,1] <- 2
  # }else{
  #   beta <-  matrix(sample(-2:2,num_df,replace=TRUE),ncol=1)
  #   beta[1] <- -2
  # }

  #if (interval){
    mu <-  c(0,sample(-6:6,size=num_classes-1,replace=TRUE))
  #}else{
  #  mu <-  c(0,sample(2:6,size=num_classes-1,replace=TRUE))
  #}

  tau <-  c(0,sample(seq(1,3,0.2),size=num_classes-1,replace=TRUE))
  var <- tau^2+1
  params <- list(beta=beta,mu=mu,var=var,tau=tau)
 # params <- sort_args(params)
  return(params)
}

sort_args <- function(args){
  if(length(args$mu) > 2){
    mu_order <- order(args$mu,decreasing = FALSE)
    args$mu <- args$mu[mu_order]
    args$var <- args$var[mu_order]
    args$beta <- args$beta[,mu_order]
    if("tau"%in% args){
      args$tau <- args$tau[mu_order]
    }
  }
  return(args)
}

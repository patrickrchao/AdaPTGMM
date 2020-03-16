library(splines)
library(stringr)

expit <- function(x){
  return(exp(x)/(1+exp(x)))
}

weighted_mean <- function(weights,values){
  total_weight <- sum(weights)
  if(sum(is.na(weights))>0){
    browser()
  }
  if(total_weight==0 || sum(weights<0)>0){
    browser()
    stop("Warning! Invalid weights found")
    return(0)
  }
  return(sum(values*weights)/total_weight)
}

gaussian_pdf <- function(x,mean,var){
  pdf <- dnorm(x,mean,sqrt(var))
  pdf[is.na(pdf)] <- 10^(-8)
  return(pdf)
}

generate_spline <- function(x,num_df){
  spline_x <- ns(x,df=num_df-1)
  spline_x <- cbind(rep(1,length(x)),spline_x)
  beta_names <- coefficient_names("Beta",ncol(spline_x))
  colnames(spline_x) <- beta_names
  return(spline_x)
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

change_of_density_df <- function(df,radius=1) {
  return(dnorm(df["z"],df["mean"],sqrt(df[["var"]]))/
           (dnorm(-abs(df["z"])+radius,0,1)-dnorm(abs(df["z"])+radius,0,1)))
}

calculate_conditional_probabilities <- function(data,est_params,params){

  num_hypo <- length(data$small_z)
  mu <- est_params$mu
  var <- est_params$var
  small_z = data$small_z
  big_z = data$big_z
  radius <- params$interval_radius

  beta <- est_params$beta
  expit_prob <- expit(data$full_x %*% beta)
  prob <- data.frame(matrix(NA, nrow=length(small_z), ncol=1))
  if(params$testing_interval){
    # Iterate over s,b,-s,-b
    for(a in params$all_a){
      for(class in 0:(params$num_classes-1)){
        name <- paste0(a,"_",class)
        sign <- 1

        if(str_detect(a,"neg")){
          sign <- -1
        }

        if(str_detect(a,"s")){
          prob[,name] <- change_of_density(sign*small_z, radius, mu[class+1], var[class+1])
        }else{
          prob[,name] <- change_of_density(sign*big_z,   radius, mu[class+1], var[class+1])/params$zeta
        }

      }
    }
  } else {
    prob[,"s_0"] <- rep(1,length(small_z))
    prob[,"b_0"] <- rep(1,length(small_z)) / params$zeta
    for(a in params$all_a){
      for(class in 1:(params$num_classes-1)){
        name <- paste0(a,"_",class)
        if(str_detect(a,"s")){
          prob[,name] <- gaussian_pdf(small_z, mu[class+1], var[class+1]) / gaussian_pdf(small_z, 0, 1)
        }else{
          prob[,name] <- gaussian_pdf(big_z, mu[class+1], var[class+1]) / gaussian_pdf(big_z, 0, 1)/params$zeta
        }

      }

    }
  }


  prob <- prob[-c(1)]
  # prob$small_prob_null <- small_prob_null
  # prob$big_prob_null <- big_prob_null
  # prob$small_prob_alt <- small_prob_alt
  # prob$big_prob_alt <- big_prob_alt
  prob$expit_prob <- expit_prob
  return(prob)
}

likelihood <- function(data,est_params,params,optimal_param=FALSE,w_ika=FALSE){

  if(!is.numeric(w_ika)){
    w_ika <- calculate_w(data,est_params,params)
  }
  temp <- w_ika %>% select("i","numerator") %>% group_by(i)%>%summarize(mean = mean(numerator))
  likelihood <- mean((log(temp))$mean)


  # w_ika$mean <- est_params$mu[(w_ika$k)+1]
  # w_ika$var  <- est_params$var[(w_ika$k)+1]
  # #w_ika$class_prob <- prob_classes[(w_ika$k)+1]
  # w_ika$density <- apply(w_ika[c("z","mean","var")], 1, change_of_density_df)
  # likelihood <- mean(w_ika$w_ia*w_ika$density)
  if(optimal_param){
    print(paste0("Likelihood with True Parameters: ",likelihood))
  }else{
    print(paste0("Calculated Likelihood: ",likelihood))
  }
  return(likelihood)
  # #,z="z", mean = "mean", var="var")
  # #apply(w_ika$z, 1, change_of_density,radius = params$interval_radius, mean = w_ika$mean, var=w_ika$var)
  #
  # w_ia <- filter(w_ika,a=="s",k==1)%>% select("w_ia")
}

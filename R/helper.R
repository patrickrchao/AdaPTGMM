library(splines)
library(stringr)
expit <- function(x){
  return(exp(x)/(1+exp(x)))
}

weighted_mean <- function(weights,values){
  total_weight <- sum(weights)
  if(total_weight==0 || sum(weights<0)>0){
    browser()
    print("Warning! Invalid weights found")
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
  x_names <- paste("X",0:(ncol(spline_x)-1),sep="")
  colnames(spline_x) <- x_names
  return(spline_x)
}

x_col_names <- function(x_num_dim){
  x_names <- paste("X",0:(x_num_dim-1),sep="")
  return(x_names)
}

inverse = function (f, lower = 0, upper = 100) {
  function (y) uniroot((function (x) f(x) - y), lower = lower, upper = upper)[1]
}

change_of_density <- function(z, radius, mean, var) {
  return(gaussian_pdf(z, mean, var)  /
           (gaussian_pdf(-abs(z) +radius,0,1) - gaussian_pdf(abs(z) + radius,0,1)))
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
      for(class in c(0,1)){
        name <- paste0(a,"_",class)
        sign <- 1

        if(str_detect(a,"neg")){
          sign <- -1
        }

        if(str_detect(a,"s")){
          prob[,name] <- change_of_density(sign*small_z, radius, mu*class, var*class+(1-class))
        }else{
          prob[,name] <- change_of_density(sign*big_z,   radius, mu*class, var*class+(1-class))/params$zeta
        }

      }
    }
  } else {
    prob[,"s_0"] <- rep(1,length(small_z))
    prob[,"b_0"] <- rep(1,length(small_z)) / params$zeta
    for(a in params$all_a){
      for(class in c(1)){
        name <- paste0(a,"_",class)
        if(str_detect(a,"s")){
          prob[,name] <- gaussian_pdf(small_z, mu, var) / gaussian_pdf(small_z, 0, 1)
        }else{
          prob[,name] <- gaussian_pdf(big_z, mu, var) / gaussian_pdf(big_z, 0, 1)/params$zeta
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
# z_to_p_values <- function(z,testing_interval,intervals){
#   if(!testing_interval){
#     p_values <- -qnorm(z)
#   }else{
#     interval_center = (intervals[,"left_end"]+intervals[,"right_end"])/2
#     radius = intervals[,"right_end"]-interval_center
#     centered_z = abs(z-interval_center)
#
#     p_values = pnorm(centered_z+radius,lower.tail=FALSE)+pnorm(-centered_z+radius)
#
#   }
#   return(p_values)
# }
#
# p_values_to_z <- function(z,testing_interval,intervals){
#   if(!testing_interval){
#     z <- -qnorm(z)
#   }else{
#     interval_center = (intervals[,"left_end"]+intervals[,"right_end"])/2
#     radius = intervals[,"right_end"]-interval_center
#     centered_z = abs(z-interval_center)
#
#     p_values = pnorm(centered_z+radius,lower.tail=FALSE)+pnorm(-centered_z+radius)
#
#   }
#   return(p_values)
# }

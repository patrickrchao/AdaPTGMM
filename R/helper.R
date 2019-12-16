library(splines)

expit <- function(x){
  return(exp(x)/(1+exp(x)))
}

weighted_mean <- function(weights,values){
  total_weight <- sum(weights)
  if(total_weight==0 || sum(weights<0)>0){
    print("Warning! Invalid weights found")
    return(0)
  }
  return(sum(values*weights)/total_weight)
}

gaussian_pdf <- function(x,mean,var){
  return(dnorm(x,mean,sqrt(var)))
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
  return((gaussian_pdf(z, mean, var) + gaussian_pdf(-z, mean, var)) /
           (gaussian_pdf(-z +radius,0,1) - gaussian_pdf(z + radius,0,1)))
}

calculate_probabilities <- function(data,est_params,params){

  num_hypo <- length(data$small_z)
  mu <- est_params$mu
  var <- est_params$var
  small_z = data$small_z
  big_z = data$big_z
  radius <- params$interval_radius

  beta <- est_params$beta
  expit_prob <- expit(data$full_x %*% beta)

  small_prob_null <- ifelse(params$testing_interval,
                            change_of_density(small_z, radius, 0, 1),
                            1)
  big_prob_null <- ifelse(params$testing_interval,
                          change_of_density(big_z, radius, 0, 1),
                          1) / params$zeta

  small_prob_alt <- ifelse(rep(params$testing_interval,num_hypo),
                           change_of_density(small_z, radius, mu, var),
                           gaussian_pdf(small_z, mu, var) /gaussian_pdf(small_z, 0, 1))

  big_prob_alt <- ifelse(rep(params$testing_interval,num_hypo),
                         change_of_density(big_z, radius, mu, var),
                         gaussian_pdf(big_z, mu, var) / gaussian_pdf(big_z, 0, 1)) / params$zeta
  prob <- list()
  prob$small_prob_null <- small_prob_null
  prob$big_prob_null <- big_prob_null
  prob$small_prob_alt <- small_prob_alt
  prob$big_prob_alt <- big_prob_alt
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

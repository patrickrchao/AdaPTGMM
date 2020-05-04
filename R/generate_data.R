library(Hmisc)
# Currently only works with num_dim = 1
generate_data <- function(num_samples=1000,num_df= 20,beta = FALSE,mu = FALSE,tau=FALSE,num_classes = 2,interval=FALSE,intervals=c(-1,1)){
  #browser()
  stopifnot(num_classes >= 2)

  params <- initialize_params(num_classes,num_df)
  if(as.logical(beta)){
    true_beta <- beta
  }else{
    true_beta <- params$beta
  }

  if(as.logical(mu)){
    true_mu <- mu
  }else{
    true_mu <-  params$mu
  }

  if(as.logical(tau)){
    true_tau <- tau
  }else{
    true_tau <-  params$tau
  }

  x <- sort(runif(num_samples)-0.5)
  spline_x <- generate_spline(x,num_df)
  colnames(spline_x) <- paste("X",0:(ncol(spline_x)-1),sep="")

 # if(num_classes == 2){
  #  gamma <- rbinom(num_samples,1,expit(spline_x%*%true_beta))
  #   browser()
  # }else{

  class_probs <- probability_from_spline(spline_x, true_beta, num_classes)


  # Subtract 1 since classes are zero-indexed
  gamma <- rMultinom(class_probs,1) - 1
 # }

  # Theta is zero if gamma is zero
  # Otherwise it is drawn from a Gaussian
  theta <- (gamma > 0) * rnorm(n=num_samples,true_mu,true_tau)
  z <- rnorm(n = num_samples, mean = theta)
  if(interval){
    interval_center = (intervals[1]+intervals[2])/2
    radius = intervals[2]-interval_center
    centered_z = abs(z-interval_center)
    z_to_p <- function(centered_z,radius) pnorm(centered_z+radius,lower.tail=FALSE)+pnorm(-centered_z+radius)
    z_to_p_rad <- function(centered_z) z_to_p(centered_z,radius)
    p_values <- z_to_p_rad(centered_z)
  }else{

    p_values <- pnorm(z,lower.tail=FALSE)
  }
  #all_data <- data.frame(x,gamma,theta,z,p_values)


  known <- list(p_values = p_values, z = z, x = x, spline_x = spline_x)
  unknown <- list(beta = true_beta, mu = true_mu, var = true_tau^2+1, gamma = gamma, theta = theta)

  return(list(known = known, unknown = unknown))
}

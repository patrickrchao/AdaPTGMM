library(Hmisc)
# Currently only works with num_dim = 1
generate_data <- function(num_samples=1000,num_df= 20,beta = FALSE,mu = FALSE,tau=FALSE,num_classes = 2){
  #browser()
  stopifnot(num_classes >= 2)

  estimates <- initialize_estimates(num_classes,num_df)
  if(as.logical(beta)){
    true_beta <- beta
  }else{
    true_beta <- estimates$beta
  }

  if(as.logical(mu)){
    true_mu <- mu
  }else{
    true_mu <-  estimates$mu
  }

  if(as.logical(tau)){
    true_tau <- tau
  }else{
    true_tau <-  estimates$tau
  }

  x <- sort(runif(num_samples)-0.5)
  spline_x <- generate_spline(x,num_df)
  if(num_classes == 2){
    gamma <- rbinom(num_samples,1,expit(spline_x%*%true_beta))
  }else{

    class_probs <- exp(spline_x%*%true_beta)
    class_probs <- class_probs / rowSums(class_probs)
    gamma <- rMultinom(class_probs,1)-1
  }

  theta <- (gamma>0)*rnorm(n=num_samples,true_mu,true_tau)
  z <- rnorm(num_samples,theta)
  p_values <- 1-pnorm(z)
  p_values <- pmax(p_values,1e-12)
  p_values <- pmin(p_values,1-1e-12)
  all_data <- data.frame(x,gamma,theta,z,p_values)
  x_names <- paste("X",0:(ncol(spline_x)-1),sep="")
  colnames(spline_x) <- x_names
  #plot(x,expit(spline_x%*% true_beta),type="l",lty=2,lwd=3,ylab="Probability of Class 1",main="True Spline")
  known <- list(p_values=p_values, z = z, x=x, spline_x=spline_x)
  unknown <- list(beta=true_beta,mu=true_mu,var=true_tau^2+1,gamma=gamma,theta=theta)
  return(list(known=known,unknown=unknown))
}


# Currently only works with num_dim = 1
generate_data <- function(num_samples=1000,num_dim=1,num_df= 20,beta = FALSE,mu = FALSE,tau=FALSE){
  #browser()
  stopifnot(num_dim == 1)


  if(as.logical(beta)){
    true_beta <- beta
  }else{

    true_beta <-  matrix(sample(-2:2,num_dim*(num_df),replace=TRUE),ncol=1)
    true_beta[1] <- -2
  }

  if(as.logical(mu)){
    true_mu <- mu
  }else{

    true_mu <-  sample(2:5,size=1)
  }

  if(as.logical(tau)){
    true_tau <- tau
  }else{
    true_tau <-  sample(1:3,size=1)
  }


  x <- sort(runif(num_dim*num_samples)-0.5)
  spline_x <- generate_spline(x,num_df)
  gamma <- rbinom(num_samples,1,expit(spline_x%*%true_beta))
  theta <- gamma*rnorm(n=num_samples,true_mu,true_tau)
  z <- rnorm(num_samples,theta)
  p_values <- 1-pnorm(z)
  p_values <- pmax(p_values,1e-10)
  p_values <- pmin(p_values,1-1e-10)
  all_data <- data.frame(x,gamma,theta,z,p_values)
  x_names <- paste("X",0:(ncol(spline_x)-1),sep="")
  colnames(spline_x) <- x_names
  plot(x,expit(spline_x%*% true_beta),type="l",lty=2,lwd=3,ylab="Probability of Class 1",main="True Spline")
  known <- list(p_values=p_values, z = z, x=x, spline_x=spline_x)
  unknown <- list(beta=true_beta,mu=true_mu,var=true_tau^2+1,gamma=gamma,theta=theta)
  return(list(known=known,unknown=unknown))
}

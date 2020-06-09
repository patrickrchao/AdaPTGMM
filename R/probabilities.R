#' Compute Class Probabilities
#'
#' @param beta fitted multinom
#' @param nclasses number of Gaussian mixture classes
#'
#' @return Dataframe of class probabilities, columns represent classes and rows represent each hypothesis
#'
class_prob <- function(beta,nclasses){

  prob <- fitted(beta)
  # Divide by the number of classes since the input data uses n*nclasses data points
  # fitted(beta) repeats predictions for various classes

  prob <- prob[1:(nrow(prob)/nclasses),]

  if(ncol(prob)==1){
    prob <- cbind(1-prob,prob)
  }

  return(prob)
}

#' Compute the probability of big/small and masked p-value conditioned on class
#'
#' @param model
#'


#' Helper function to compute probability of big/small and masked p-value conditioned on class
#' for specific hypothesis
#'
#' @param z Test statistic value
#' @param mean mean of Gaussian mixture model class k
#' @param var variance of Gaussian mixture model class k
#'
#' @return P[a_i=a,\tilde p_i | \gamma=k]
#'
#' @details Computation is equivalent to phi(z,mean,var)/phi(z,0,1) where phi(z,mu,tau) is the density of a Gaussian
#' random variable with mean mu and variance tau at z.
prob_jacobian <- function(z, mean, var) {
  return(exp(z^2/2-(z-mean)^2/(2*var))/sqrt(var))
}


#' Marginalize over variable
#'
#'

#' Compute Class Probabilities
#'
#' @param beta fitted multinom
#' @param nclasses number of Gaussian mixture classes
#' @param beta_formula beta formula for model
#' @param class_probabilities
#' @return Dataframe of class probabilities, columns represent classes and rows represent each hypothesis
#' @noRd
class_prob <- function(beta_model,nclasses,n,model_type,x=NULL){
  #
  if(model_type == "glm"){
    if(is.null(x)){
      prob <- fitted(beta_model)
    }else{

      prob <- predict(beta_model,type="probs",newdata=x)
    }
  }else if(model_type == "gam"){
    prob <- predict(beta_model,type="response")
  }else if(model_type == "mgcv"){
    prob <- predict(beta_model,type="response")
  }
    # Divide by the number of classes since the input data uses n*nclasses data points
    # fitted(beta) repeats predictions for various classes
  if(nclasses == 2){
    prob <- cbind(1-prob,prob)
  }
  if(n != nrow(prob)){
    prob <- prob[1:(nrow(prob)/nclasses),]
  }
  prob <- pmax(pmin(prob,1 - 1e-12), 1e-12)
  return(prob)
}

#' Helper function to compute probability of big/small and masked p-value conditioned on class in one sided case
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
#' @noRd
prob_jacobian_one_sided <- function(z, mean, var) {
  return(exp(z^2/2-(z-mean)^2/(2*var))/sqrt(var))
}


#' Helper function to compute probability of big/small and masked p-value conditioned on class in interval case
#' for specific hypothesis
#'
#' @param z Test statistic value
#' @param mean mean of Gaussian mixture model class k
#' @param var variance of Gaussian mixture model class k
#'
#' @return P[a_i=a,\tilde p_i | \gamma=k]
#'
#' @details Computation is  phi(z,mean,var)/[phi(-abs(z)+r,0,1)-phi(abs(z)+r,0,1)] where phi(z,mu,tau) is
#' the density of a Gaussian random variable with mean mu and variance tau at z.
#' @noRd
prob_jacobian_interval <- function(z, mean, var, radius) {
    return(dnorm(z,mean,sqrt(var))/
            (dnorm(-abs(z) + radius,0,1) -
             dnorm(abs(z) + radius,0,1))
          )
}


#' Marginalize over variable
#'
#' @param w_ika w_ika dataframe
#' @param margin_vars variables to marginalize over
#'
#' @return w_ika data frame marginalized over variables
#' @noRd
marginalize <- function(w_ika, margin_vars){

  all_vars <- c("i","class","a")
  group_by_vars <- setdiff(all_vars,margin_vars)
  new_w_ika <- w_ika[, .(value=sum(value)), by = group_by_vars]
  #groups <-   dplyr::group_by(w_ika,.dots=group_by_vars)
  #new_w_ika <- dplyr::ungroup(dplyr::summarise(groups,value=sum(value)))

  # Avoid piping for speed
  # Equivalent to
  # new_w_ika  <- w_ika %>%
  #  dplyr::group_by_(.dots=group_by_vars) %>% #groupby variables
  #  dplyr::summarise(value=sum(value))%>% # sum across groups
  #  dplyr::ungroup()
  return(new_w_ika)
}

#' log likelihood function
#'
#' @param model model class
#'
#' @return Returns log likelihood of data given model
#' @noRd
log_likelihood <- function(model){

  w_ika_sums <- e_step_w_ika(model, include_z = FALSE, agg_over_hypotheses = TRUE)
  # w_ika_sums corresponds to sum_{k,a} {P[a_ia=a,\tilde p_i|\gamma_i=k] P[\gamma_i=k|x_i]*\zeta^1{a_ia=b}}
  log_like <- sum(log(w_ika_sums$value))
  return(log_like)
}

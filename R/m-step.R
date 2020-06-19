#' Perform Maximization step for beta
#'
#' @param model model variable containing data, args, params
#' @param gammas Estimated gammas from e_step_gammas, dataframe where columns correspond to classes
#' and rows correspond to hypotheses. Each row sums to 1.
#'
#' @details Performs a multinomial logistic regression from nnet class
#' @return model
#' @noRd
m_step_beta <- function(model,gammas){

  ndf <- model$args$ndf
  nclasses <- model$args$nclasses
  if(ndf == 1){
    # If there is only a single degree of freedom (intercept only model),
    # we may estimate the class probabilities by considering the mean for each class.
    grouped <- dplyr::group_by(gammas,class)
    beta <- dplyr::ungroup(dplyr::summarise(grouped, value=mean(value)))$value
    beta <- beta/sum(beta)
    if(sum(is.na(beta))>0){
      browser()
    }
    if(sum(beta==0)>0){
      browser()
    }
  }else{

    x <- model$data$full_x
    multinom_data <- data.frame(x,gammas)
    x_colnames <- colnames(multinom_data)[seq(ndf)]
    formula <- paste0("class ~ ",paste0(x_colnames,collapse = " + ")," -1")
    # if the multinom beta model exists, use the previous weights as the starting point for faster convergence
    if(!is.null(model$params$beta)){
      beta <- nnet::multinom(formula, multinom_data, weights = value, trace = FALSE,maxit=5,Wts=model$params$beta$wts)
    }else{
      beta <- nnet::multinom(formula, multinom_data, weights = value, trace = FALSE,maxit=100)
    }

  }

  # Update class probabilities
  model$data$class_prob <- class_prob(beta,nclasses)
  # Update beta
  model$params$beta <- beta

  return(model)
}

#' Perform Maximization step for mu and tau
#' @noRd
m_step_mu_tau <- function(model,w_ika){
  args <- model$args
  params <- model$params
  data <- model$data

  z <- w_ika$z
  for (k in 1:(args$nclasses-1)){
    subset <- w_ika[w_ika$class == k,]
    params$mu[k+1] <- weighted_mean(subset$z,subset$value)
    # Minimum variance for convolved Gaussian is 1
    params$var[k+1] <- max(weighted_mean((subset$z-params$mu[k+1])^2,subset$value), 1)
  }
  if(sum(is.na(params$mu))>0 | sum(is.na(params$var))>0){
    browser()
  }
  return(params)
}


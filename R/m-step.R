#' Perform Maximization step for beta
#'
#' @param model model variable containing data, args, params
#' @param gammas Estimated gammas from e_step_gammas use as weights, dataframe where columns correspond to classes
#' and rows correspond to hypotheses. Each row sums to 1.
#'
#' @details Performs a multinomial logistic regression
#' @return model
#' @noRd
m_step_beta <- function(model,gammas){

  model_type <- model$args$model_type
  nclasses <- model$args$nclasses
  beta_formula <- model$args$beta_formula

  x <- model$data$x
  rownames(x) <- NULL

  out <- m_step_beta_defaults(model_type,beta_formula,x,gammas,model$params$beta)

  model$params$beta <- out$model_weights
  model$params$df <- out$df
  prob <-  out$fitted_prob
  if(nclasses == 2 & ncol(prob) == 1){
    prob <- cbind(1-prob,prob)
  }
  if(model$args$n != nrow(prob)){

    prob <- prob[1:model$args$n,]
  }
  prob <- as.data.frame(pmax(pmin(as.matrix(prob), 1 - 1e-12), 1e-12))
  model$data$class_prob <- prob

  return(model)
}

#' Perform Maximization step for mu and tau
#' @noRd
m_step_mu_tau <- function(model,w_ika){
  args <- model$args
  params <- model$params
  data <- model$data
  z <- w_ika$z

  for (k in 1:args$nclasses){
    subset <- w_ika[w_ika$class == k,]
    se <- data$se[subset$i]

    # Optim for mu and tau^2
    if(args$symmetric_modeling){
      sol <- optim_symmetric_solver(params$mu[k],params$var[k],subset$z,subset$value,se)
    }else{
      sol <- optim_solver(params$mu[k],params$var[k],subset$z,subset$value,se)
    }

    params$mu[k] <- sol$mu
    params$var[k] <- sol$var

  }


  if(any(is.na(params$mu)) | any(is.na(params$var))){

    stop("NA value for mu or variance found. Stopping.")
  }

  return(params)
}

optim_solver <- function(init_mu,init_var,z,w_ika,se){
  sol <- optim(c(init_mu,init_var),function(x){
    mu <- x[1]
    var <- x[2]

    mean(w_ika*(log(se^2+var) + (z-mu)^2/(var+se^2)))
  },gr <- function(x){
    mu <- x[1]
    var <- x[2]
    c(-mean(w_ika * 2*(z-mu) / (var+se^2)),
      mean(w_ika / (var + se^2)) -
        mean(w_ika * (z-mu)^2 / (var+se^2)^2)
    )
  },method="L-BFGS-B",control=list(maxit=10,factr=1e-2),lower=c(-Inf,0),upper=c(10,30))
  return(list(mu=sol$par[1],var=sol$par[2]))
}


optim_symmetric_solver <- function(init_mu,init_var,z,w_ika,se){

  sol <- optim(c(init_mu,init_var),function(x){
    mu <- x[1]
    var <- x[2]

    -mean(w_ika*( log(dnorm(z,mu,sqrt(var+se^2)) + dnorm(z,-mu,sqrt(var+se^2)))))
  },gr <- function(x){
    mu <- x[1]
    var <- x[2]
    c(-mean(w_ika * ( z * tanh(mu*z/(se^2+var))-mu)/(se^2+var)),
      -mean(w_ika * (-1/2 / (se^2+var) +
                       (mu^2 - 2 * mu * z * tanh(mu*z/(se^2+var)) + z^2)/(2*(se^2+var)^2)) )
    )
  },method="L-BFGS-B",control=list(maxit=10,factr=1e-2),lower=c(0,0),upper=c(10,30))
  return(list(mu=sol$par[1],var=sol$par[2]))
}



fit_parameters <- function(model,logging=FALSE){
  if(logging){
    return(fit_parameters_logging(model))
  }else{
    data <- model$data
    params <- model$params
    args <- model$args

    for(i in seq(args$iterations)){
      params <- update_parameters(data,params,args)
    }

    model$params <- params
    out <- list(model=model)

    return(out)
  }
}

fit_parameters_logging <- function(model,beta_seq=data.frame(),mu_seq=data.frame(),var_seq=data.frame()){
  data <- model$data
  params <- model$params
  args <- model$args

  if("Trial" %in% colnames(beta_seq)){
    trial = max(beta_seq$Trial)+1
  }else{
    trial = 0
  }
  beta_names <- c(coefficient_names("Beta",ncol(data$full_x)),"Class")
  mu_names <- coefficient_names("Mu",args$num_classes)
  var_names <- coefficient_names("Var",args$num_classes)
  for (curr_class in 0:(args$num_classes-1)){
    beta_class <- params$beta[,curr_class+1]
    beta_seq <- rbind(beta_seq, c(beta_class, curr_class,0, trial))
  }
  #beta_seq <- rbind(beta_seq,c(params$beta,0,trial))
  mu_seq <-   rbind(mu_seq,  c(params$mu,  0,trial))
  var_seq <-  rbind(var_seq, c(params$var, 0,trial))
  # can speed this up
  for(i in seq(args$iterations)){

    params <- update_parameters(data,params,args)
    for (curr_class in 0:(args$num_classes-1)){

      beta_class <- params$beta[,curr_class+1]
      beta_seq <- rbind(beta_seq, c(beta_class, curr_class,i, trial))
    }
    mu_seq <-   rbind(mu_seq,   c(params$mu,   i, trial))
    var_seq <-  rbind(var_seq,  c(params$var,  i, trial))

  }

  colnames(beta_seq) <- c(beta_names,"Iteration","Trial")
  colnames(mu_seq) <- c(mu_names,"Iteration","Trial")
  colnames(var_seq) <- c(var_names,"Iteration","Trial")
  out <- list()
  model$params <- params

  out$model <- model
  out$best_beta = params$beta
  out$beta_seq = beta_seq

  out$best_mu = params$mu
  out$mu_seq = mu_seq

  out$best_var = params$var
  out$var_seq = var_seq
  return(out)
}



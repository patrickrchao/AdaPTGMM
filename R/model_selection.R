#' Model selection for degrees of freedom and number of classes
#' @param data data of model
#' @param args args of model
#' @param beta_formulas List of formulas for beta expansion
#' @param nclasses_list Vector for possible number of classes
#' @param selection Criteria for model selection
#' @param intercept_model Boolean whether to include intercept model
#' @param initialization Initialization procedure, either kmeans or random
#' @param training_proportion Proportion of data used for training, in the \code{selection}='\code{cross_validation}' case.
#'
#' @return Initialized and pretrained model with best performance based on selection criterion
#' @keywords Internal
model_selection <- function(data,args,beta_formulas,nclasses_list,selection,intercept_model,initialization,training_proportion=0.6){

  n <- args$n
  #cat("Model selection starting. Shrink the set of candidate models if it is too time-consuming.\n")
  if(intercept_model){
    beta_formulas <- c("intercept",beta_formulas)
  }

  beta_formulas <- .check_formulas(data$x,beta_formulas)
  # Construct train/validation splits for cross_validation
  if(selection == "cross_validation"){
    datasets <- .split_data(data,n,training_proportion)
    train <- datasets$train
    valid <- datasets$valid
    new_args <- args
    new_args$n <- as.integer(n * training_proportion)
    beta_formulas <- .check_formulas(train$x,beta_formulas)
    beta_formulas <- .check_formulas(valid$x,beta_formulas)
  }else{
    train <- data
    new_args <- args
    valid <- NULL
  }


  # Construct grid of all parameter combinations
  # corresponds to indicies in beta_formulas and nclasses_list
  param_grid <-  expand.grid(1:length(beta_formulas),1:length(nclasses_list))
  colnames(param_grid) <- c("formula","nclasses")
  n_permutations <- nrow(param_grid)
  model_list <- vector("list",n_permutations)

  init_params <- lapply(nclasses_list,function(x)initialize_params(train,x,initialization))

  for(row_index in 1:n_permutations){
    row <- param_grid[row_index,]
    new_args$beta_formula <- beta_formulas[row$formula]
    new_args$nclasses <- nclasses_list[row$nclasses]

    model <- create_model(train, new_args,init_params[[row$nclasses]])
    model <- try(EM(model, preset_iter = args$niter_ms),silent=TRUE)
    if (class(model)[1] == "try-error"){
      param_grid[row_index, "log_like"] <- -Inf
      param_grid[row_index, "penalty"] <- 0
    }else{
      out <- .selection_helper(selection,model,valid,n)

      model <- out$model
      penalty <- out$penalty

      # Store all trained models
      model_list[[row_index]] <- model

      #Update grid loglikelihood and penalty
      param_grid[row_index, "log_like"] <- log_likelihood(model)
      param_grid[row_index, "penalty"] <- penalty
    }
  }

  param_grid$value <- param_grid$log_like - param_grid$penalty
  # Extract best set of parameters from grid
  max_index <- which(param_grid$value == max(param_grid$value), arr.ind = TRUE)
  beta_formula_ind <- param_grid[max_index,"formula"]
  nclasses_ind <- param_grid[max_index,"nclasses"]
  args$beta_formula <- beta_formulas[beta_formula_ind]
  args$nclasses <- nclasses_list[nclasses_ind]

  #cat("Model selection completed.\n")
  # If using cross_validation, model needs to be reinitialized with full data and pretrained
  if(selection == "cross_validation"){
    params <- model_list[[max_index]]$params
    model <- create_model(data,args,params)
    chosen_model <- EM(model,preset_iter = args$niter_ms)
  }else{
    chosen_model <- model_list[[max_index]]
  }
  return(chosen_model)
}

#' Split data for cross validation
#'
#' @param data data class
#' @param n number of total data points
#' @param training_proportion proportion of data for training
#'
#' @noRd
.split_data <- function(data,n,training_proportion){
  random_order <- sample(1:n,size = n)
  split_point <- as.integer(n*training_proportion)

  split_func <- function(y,indices){
    if(is.vector(y)){
      return(y[indices])
    }else{
      out <- data.frame(y[indices,])
      colnames(out) <- colnames(y)
      rownames(out) <- NULL
      return(out)
    }
  }

  train <- lapply(data,function(x)split_func(x,random_order[1:split_point]) )
  valid <- lapply(data,function(x)split_func(x,random_order[(split_point+1):n]) )

  datasets <- list("train"= train, "valid" = valid)
  return(datasets)
}


#' Helper function for model selection
#' For cross validation, returns model for validation set
#' Updates class probability correctly
#'
#' For AIC/BIC, computes penalty
#'
#' Returns model and penalty
#'
#' @param selection AIC/BIC/cross_validation
#' @param model trained model
#' @param valid Validation set (null for AIC/BIC)
#' @param total_n total number of hypotheses including training and validation
#'
#' @return list of model and penalty
#' @noRd
.selection_helper <- function(selection, model, valid, total_n){
  if(selection == "cross_validation"){
    trained_params <- model$params
    valid_args <- model$args
    valid_args$n <- as.integer(total_n - model$args$n)
    model <- create_model(valid, valid_args,trained_params)

    prob <- class_prob(model$params$beta,model$args$nclasses,model$data$x)
    model$data$class_prob <- prob
    penalty <- 0
  } else {
    # df degrees of freedom in beta
    # (class-1)*2 degrees of freedom in mu and tau^2+1
    if(model$args$beta_formula == "intercept"){
      df <- 1
    }else{
      df <- model$params$beta$edf
    }

    d <-  (model$args$nclasses - 1) * (df + 2)
    # BIC is divided by -2 so that large values are desirable
    if(selection == "BIC"){
      penalty <- log(total_n)*d/2
      # AIC is divided by -N/2 so that large values are desirable
    } else if(selection == "AIC"){
      penalty <- d
    } else{
      stop("Invalid model selection method.")
    }
  }
  return(list(model=model,penalty=penalty))
}

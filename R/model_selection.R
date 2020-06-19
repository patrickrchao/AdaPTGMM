#' Model selection for degrees of freedom and number of classes
#' @param data data of model
#' @param args args of model
#' @param ndf_list Vector for possible degrees of freedom
#' @param nclasses_list Vector for possible number of classes
#' @param selection Criteria for model selection
#' @param training_proportion Proportion of data used for training, in the \code{selection}='\code{cross_validation}' case.
#'
#' @return Initialized and pretrained model with best performance based on selection criterion
#' @keywords Internal
model_selection <- function(data,args,ndf_list,nclasses_list,selection,training_proportion=0.6){

  n <- args$n
  print("Model selection starting. Shrink the set of candidate models if it is too time-consuming.")

  # Construct grid of all parameter combinations
  param_grid <-  expand.grid(ndf_list,nclasses_list)
  colnames(param_grid) <- c("ndf","nclasses")
  n_permutations <- nrow(param_grid)
  model_list <- vector("list",n_permutations)

  # Construct train/validation splits for cross_validation
  if(selection == "cross_validation"){
    datasets <- split_data(data,n,training_proportion)
    train <- datasets$train
    valid <- datasets$valid
    new_args <- args
    new_args$n <- as.integer(n * training_proportion)
  }else{
    train <- data
    new_args <- args
  }

  for(row_index in 1:n_permutations){
    row <- param_grid[row_index,]
    new_args$ndf <- row$ndf
    new_args$nclasses <- row$nclasses
    #print(paste(row$ndf,row$nclasses))
    model <- create_model(train, new_args)
    model <- EM(model, preset_iter = args$niter_ms)

    # Store all trained models
    model_list[[row_index]] <- model

    if(selection == "cross_validation"){
      trained_params <- model$params
      valid_args <- new_args
      valid_args$n <- as.integer(n - new_args$n)
      model <- create_model(valid, valid_args)
      model$params <- trained_params

      if(model$args$ndf == 1){
        prob <- data.frame(t(model$params$beta))
      }else{
        prob <- predict(model$params$beta,type="probs",newdata=data.frame(model$data$full_x))
        if(model$args$nclasses == 2){
          prob <- cbind(1-prob,prob)
        }
      }

      model$data$class_prob <- prob

      penalty <- 0
    } else {
      # df degrees of freedom in beta
      # (class-1)*2 degrees of freedom in mu and tau^2+1
      d <- model$args$ndf + (model$args$nclasses - 1)*2
      # BIC is divided by -2 so that large values are desirable
      if(selection == "BIC"){
        penalty <- log(n)*d/2
      # AIC is divided by -N/2 so that large values are desirable
      } else if(selection == "AIC"){
        penalty <- d
      } else{
        stop("Invalid model selection method.")
      }
    }

    #Update grid loglikelihood and penalty

    param_grid[row_index, "log_like"] <- log_likelihood(model)
    param_grid[row_index, "penalty"] <- penalty
  }

  param_grid$value <- param_grid$log_like - param_grid$penalty
  # Extract best set of parameters from grid
  max_index <- which(param_grid$value == max(param_grid$value), arr.ind = TRUE)
  ndf <- param_grid[max_index,"ndf"]
  nclasses <- param_grid[max_index,"nclasses"]
  args$ndf <- ndf
  args$nclasses <- nclasses

  print(paste("Model selection completed,",ndf,"degree(s) of freedom and",nclasses,"classes."))
  # If using cross_validation, model needs to be reinitialized with full data and pretrained
  if(selection == "cross_validation"){
    model <- create_model(data,args)
    model$params <- model_list[[max_index]]$params
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
#' @noRd
split_data <- function(data,n,training_proportion){
  random_order <- sample(1:n,size = n)
  split_point <- as.integer(n*training_proportion)
  train <- lapply(data,function(x) x[random_order[1:split_point]] )
  valid <- lapply(data,function(x) x[random_order[(split_point+1):n]] )
  datasets <- list("train"= train, "valid" = valid)
  return(datasets)
}

#' Model selection for degrees of freedom and number of classes
#' @param data data of model
#' @param args args of model
#' @param beta_formulas List of formulas for beta expansion
#' @param nclasses_list Vector for possible number of classes
#' @param cr Criterion for model selection
#' @param initialization Initialization procedure, either kmeans or random
#'
#' @return Initialized and pretrained model with best performance based on selection criterion
#' @keywords Internal
model_selection <- function(data,args,beta_formulas,nclasses_list,cr,initialization,param_grid=NULL){
  n <- args$n
  niter_ms <- args$niter_ms
  cat("Model selection starting. Shrink the set of candidate models if it is too time-consuming.\n")

  # Construct grid of all parameter combinations
  # corresponds to indices in beta_formulas and nclasses_list
  if(is.null(param_grid)){
    param_grid <-  expand.grid(1:length(beta_formulas),1:length(nclasses_list))
    colnames(param_grid) <- c("formula","nclasses")
  }

  n_permutations <- nrow(param_grid)

  init_params <- lapply(nclasses_list,function(x)initialize_params(data,x,initialization))
  best_value <- Inf
  best_index <- 0

  pb <- txtProgressBar(min = 0, max = n_permutations, style = 3, width = 50)
  for(row_index in 1:n_permutations){
    row <- param_grid[row_index,]
    new_args <- args
    new_args$beta_formula <- beta_formulas[[row$formula]]
    new_args$nclasses <- nclasses_list[row$nclasses]

    model <- create_model(data, new_args,init_params[[row$nclasses]])
    model <- try(EM(model, preset_iter = niter_ms),silent=TRUE)

    if (class(model)[1] == "try-error"){
      #If this is the first time this formula has been encountered
      if(row$nclasses == 1){
        warning(paste0("Invalid beta formula found: ",paste(format(new_args$beta_formula)),". Ignoring formula."))
      }
      param_grid[row_index, "log_like"] <- -Inf
      param_grid[row_index, "penalty"] <- 0
    }else{
      out <- .selection_helper(cr,model)

      #Update grid loglikelihood and penalty
      log_like <- out$log_like
      value <- out$value
      df <- out$df
      param_grid[row_index, "log_like"] <- log_like
      param_grid[row_index, "value"] <- value
      param_grid[row_index, "df"] <- df
      if( value < best_value){
        best_model <- model
        best_value <- value
        best_df <- df
      }

    }
    setTxtProgressBar(pb, row_index)
  }

  if(best_value == Inf){
    stop("All beta formula models are invalid.")
  }
  best_model$params$value <- best_value
  best_model$params$full_df <- best_df
  cat("\n")
  return(best_model)
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
.selection_helper <- function(cr, model){
  nclasses <- model$args$nclasses
  d <- NULL
  # df degrees of freedom in beta
  if(cr == "spread"){
    probs <- big_over_small_prob(model)
    probs <- probs/(probs+1)
    probs <- probs[model$data$mask]
    hist(probs,xlim = c(0,1))
    # Multiply by negative 1 since smaller values are more desirable
    value <- -1 * mean((probs-0.5)^2)
    log_like <- NA
  }else if(cr == "c entropy"){
    probs <- big_over_small_prob(model)
    probs <- probs/(probs+1)
    probs <- probs[model$data$mask]
    value <- -1* mean(round(probs)*log(probs))#mean((probs-0.5)^2)
    log_like <- NA
  }else if(cr == "cheating"){
    probs <- big_over_small_prob(model)
    probs <- probs/(probs+1)

    probs <- probs[model$data$mask]
    labels <-(model$data$a == "b")[model$data$mask]
    hist(probs,xlim = c(0,1))
    value <- -1 * mean(labels*log(probs))#mean((probs-0.5)^2)
    log_like <- NA
  }else{
    log_like <- log_likelihood(model)
    df <- model$params$df
    # nclasses*2 degrees of freedom in mu and tau^2+1
    d <-  nclasses * 2 + df
    n <- model$args$n
    value <- info_cr(log_like,cr,d,n)
  }


  return(list(log_like=log_like,value=value,df=d))
}


# From Lihua Lei adaptMT
# https://github.com/lihualei71/adaptMT/blob/34d2f183e6cbb9842ce329da250a1de1d586b648/R/EM-mix-ms.R
info_cr <- function(log_like, cr, df, n){
  switch(cr,
         "AIC" = 2 * df - 2 * log_like,
         "AICC" = 2 * df * n / (n - df - 1),
         "BIC" = log(n) * df - 2 * log_like,
         "HIC" = 2 * log(log(n)) * df - 2 * log_like,
         NULL)
}



#
#
# model_selection <- function(data,args,beta_formulas,nclasses_list,selection,intercept_model,initialization,training_proportion=0.6){
#
#   n <- args$n
#   cat("Model selection starting. Shrink the set of candidate models if it is too time-consuming.\n")
#   if(intercept_model){
#     beta_formulas <- c("intercept",beta_formulas)
#   }
#   #beta_formulas <- .check_formulas(data$x,beta_formulas)
#   # Construct train/validation splits for cross_validation
#   if(selection == "cross_validation"){
#     datasets <- .split_data(data,n,training_proportion)
#     train <- datasets$train
#     valid <- datasets$valid
#     new_args <- args
#     new_args$n <- as.integer(n * training_proportion)
#     #beta_formulas <- .check_formulas(train$x,beta_formulas)
#     #beta_formulas <- .check_formulas(valid$x,beta_formulas)
#   }else{
#     train <- data
#     new_args <- args
#     valid <- NULL
#   }
#
#
#   # Construct grid of all parameter combinations
#   # corresponds to indicies in beta_formulas and nclasses_list
#   param_grid <-  expand.grid(1:length(beta_formulas),1:length(nclasses_list))
#   colnames(param_grid) <- c("formula","nclasses")
#   n_permutations <- nrow(param_grid)
#   model_list <- vector("list",n_permutations)
#
#   init_params <- lapply(nclasses_list,function(x)initialize_params(train,x,initialization))
#   pb <- txtProgressBar(min = 1, max = row_index, style = 3, width = 50)
#   for(row_index in 1:n_permutations){
#     row <- param_grid[row_index,]
#     new_args$beta_formula <- beta_formulas[[row$formula]]
#     new_args$nclasses <- nclasses_list[row$nclasses]
#
#     model <- create_model(train, new_args,init_params[[row$nclasses]])
#
#     model <- try(EM(model, preset_iter = args$niter_ms,save_model=(selection=="cross_validation")),silent=TRUE)
#
#     if (class(model)[1] == "try-error"){
#       #If this is the first time this formula has been encountered
#       if(row$nclasses == 1){
#         warning(paste0("Invalid beta formula found: ",paste(format(beta_formulas[[3]])),". Ignoring formula."))
#       }
#       param_grid[row_index, "log_like"] <- -Inf
#       param_grid[row_index, "penalty"] <- 0
#     }else{
#       out <- .selection_helper(selection,model,valid,n)
#
#       model <- out$model
#       penalty <- out$penalty
#
#       # Store all trained models
#       model_list[[row_index]] <- model
#
#       #Update grid loglikelihood and penalty
#       param_grid[row_index, "log_like"] <- log_likelihood(model)
#       param_grid[row_index, "penalty"] <- penalty
#     }
#   }
#   param_grid$value <- param_grid$log_like - param_grid$penalty
#   # Extract best set of parameters from grid
#   max_index <- which(param_grid$value == max(param_grid$value), arr.ind = TRUE)
#   beta_formula_ind <- param_grid[max_index,"formula"]
#   nclasses_ind <- param_grid[max_index,"nclasses"]
#   args$beta_formula <- beta_formulas[[beta_formula_ind]]
#   args$nclasses <- nclasses_list[nclasses_ind]
#   cat("Model selection completed.\n")
#   print(args$beta_formula)
#   # If using cross_validation, model needs to be reinitialized with full data and pretrained
#   if(selection == "cross_validation"){
#     params <- model_list[[max_index]]$params
#     model <- create_model(data,args,params)
#     chosen_model <- EM(model,preset_iter = args$niter_ms)
#   }else{
#     chosen_model <- model_list[[max_index]]
#   }
#   return(chosen_model)
# }


#' #' Split data for cross validation
#' #'
#' #' @param data data class
#' #' @param n number of total data points
#' #' @param training_proportion proportion of data for training
#' #'
#' #' @noRd
#' .split_data <- function(data,n,training_proportion){
#'   random_order <- sample(1:n,size = n)
#'   split_point <- as.integer(n*training_proportion)
#'
#'   split_func <- function(y,indices){
#'     if(is.vector(y)){
#'       return(y[indices])
#'     }else{
#'       out <- data.frame(y[indices,])
#'       colnames(out) <- colnames(y)
#'       rownames(out) <- NULL
#'       return(out)
#'     }
#'   }
#'
#'   train <- lapply(data,function(x)split_func(x,random_order[1:split_point]) )
#'   valid <- lapply(data,function(x)split_func(x,random_order[(split_point+1):n]) )
#'
#'   datasets <- list("train"= train, "valid" = valid)
#'   return(datasets)
#' }

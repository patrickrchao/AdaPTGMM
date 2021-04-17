#' Model selection for degrees of freedom and number of classes
#' @param data data of model
#' @param args args of model
#' @param beta_formulas List of formulas for beta expansion
#' @param nclasses_list Vector for possible number of classes
#' @param cr Criterion for model selection
#'
#' @return Initialized and pretrained model with best performance based on selection criterion
#' @keywords Internal
model_selection <- function(data,args,beta_formulas,nclasses_list,cr,param_grid=NULL){
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

  init_params <- lapply(nclasses_list,function(x)initialize_params(data,x,args$all_a,args$symmetric_modeling))
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
    #TODO: Add shortcircuit if only one model selected
    if (class(model)[1] == "try-error"){
      #If this is the first time this formula has been encountered
      if(row$nclasses == 1){
        warning(paste0("Invalid beta formula found: ",paste(format(new_args$beta_formula)),". Ignoring formula."))
      }
      param_grid[row_index, "log_like"] <- -Inf
      param_grid[row_index, "penalty"] <- 0
    }else{
      out <- .selection_helper(cr,model)

      #Update grid log likelihood and penalty
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
  log_like <- log_likelihood(model)
  df <- model$params$df
  # nclasses*2 degrees of freedom in mu and tau^2
  d <-  nclasses * 2 + df
  n <- model$args$n
  value <- info_cr(log_like,cr,d,n)

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



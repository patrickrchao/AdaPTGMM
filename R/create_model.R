#' Construct Model Class
#'
#' @param data data class
#' @param args args class
#' @param params args class, default NULL, only include if precomputed
#'
#' @details
#' Function initializes parameters and computes the basis expanded covariates. Combines data, args, and params
#' into model class
#'
#' @return model class with data, args, and params
#' @noRd
create_model <- function(data,args,params=NULL){
  nclasses <- args$nclasses
  if(is.null(params)){
    params <- initialize_params(data,nclasses,args$all_a,args$symmetric_modeling)
  }
  base_prob <- rep(1/(nclasses),nclasses)
  data$class_prob <- t(replicate(n=args$n,base_prob))
  model <- list(data=data, args=args, params=params)

  testing <- args$testing

  class(model) <- paste0("adaptgmm_model_",testing)

  return(model)
}

#' Construct Model Class
#'
#' @param data data class
#' @param args args class
#'
#' @details
#' Function initializes parameters and computes the basis expanded covariates. Combines data, args, and params
#' into model class
#'
#' @return model class with data, args, and params
#' @noRd
create_model <- function(data,args){
  params <- initialize_params(args)

  nclasses <- args$nclasses
  base_prob <- c(0.9,rep(0.1/(nclasses-1),nclasses-1))
  data$class_prob <- t(replicate(n=args$n,base_prob))

  model <- list(data=data, args=args, params=params)

  testing <- args$testing
  if(testing == "one_sided" | testing == "interval"){
    class(model) <- paste0("adaptgmm_model_",testing)
  }else{
    stop("Invalid testing type inputted. Valid forms of testing: `one_sided` and `interval`.")
  }


  return(model)
}

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

  data$full_x <- generate_spline(data$x,args$ndf)

  data$class_prob <- matrix(c(0.9,rep(0.1/(args$nclasses-1),args$nclasses-1)),ncol=args$nclasses,nrow=args$n,byrow=FALSE)
  model <- list(data=data, args=args, params=params)

  testing <- args$testing
  if(testing == "one_sided" | testing == "interval"){
    class(model) <- paste0("adaptgmm_model_",testing)
  }else{
    stop("Invalid testing type inputted. Valid forms of testing: `one_sided` and `interval`.")
  }


  return(model)
}

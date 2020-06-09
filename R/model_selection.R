#' Model selection for degrees of freedom and number of classes
#' @param data data of model
#' @param args args of model
#' @param ndf_list Vector for possible degrees of freedom
#' @param nclasses_list Vector for possible number of classes
#' @param selection Criteria for model selection
#' @keywords Internal
model_selection <- function(data,args,ndf_list,nclasses_list,selection){
  # Iterate over possible degrees of freedom
  # Iterate over number of classes
  args$ndf <- 3
  args$nclasses <- 3
  model <- create_model(data,args)
  return(model)
}



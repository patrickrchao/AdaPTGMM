#' Create args class
#'
#' Class containing all relevant arguments for model
#'
#' @param testing "\code{one_sided}" or "\code{interval}".
#' @param rendpoint Right endpoint of interval null
#' @param lendpoint Left endpoint of interval null
#' @param masking_params List of masking params,
#' alpha_m The maximum possible rejected p-value.
#' zeta Controls minimum possible number of rejections.
#' lambda Controls where p-values are mirrored.
#' @param masking_shape "\code{tent}" or "\code{comb}".
#' @param niter_fit Number of iterations in EM procedure for model update
#' @param niter_ms Number of iterations in EM procedure for model selection
#' @param nfits Number of model updates in AdaPT procedure
#' @param n Number of hypotheses
#' @param initialization Initialization procedure, kmeans or random
#' @param beta_formula Beta formula for model
#' @param nclasses Number of classes in Gaussian Mixture Model, minimum 2.
#'
#' @return args class
#' @noRd
construct_args <- function(testing,model_type,rendpoint,lendpoint,masking_params,masking_shape,niter_fit,niter_ms,nfits,n,initialization,beta_formula=NULL,nclasses=NULL){

  all_a <- c("s","b")
  if(testing == "one_sided"){
    z_to_p <- function(z) pnorm(z,lower.tail = FALSE)
    p_to_z <- function(p) qnorm(p,lower.tail = FALSE)
    base_prob <- base_prob_one_sided
    #jacobian <- prob_jacobian_one_sided
  }else if(testing == "two_sided"){
    z_to_p <- function(z) 2*pnorm(abs(z),lower.tail = FALSE)
    p_to_z <- function(p) qnorm(p/2,lower.tail = FALSE)
    base_prob <- base_prob_two_sided
    all_a <- c(all_a,"s_neg","b_neg")
  }
  else if(testing == "interval"){
    if(is.null(lendpoint)){
      lendpoint <- -1 * rendpoint
    }
    if(masking_shape == "tent"){
      warning("For interval testing the masking function is automatically set to a comb masking.")
      masking_shape <- "comb"
    }
    radius <-  (rendpoint - lendpoint) / 2
    z_to_p <- function(z) pnorm(abs(z)+radius,lower.tail=FALSE) + pnorm(-abs(z)+radius)
    p_to_z_inv <- .inverse(z_to_p,lower = 0)
    p_to_z <- function(z) unlist(mapply(p_to_z_inv,z))
    all_a <- c(all_a,"s_neg","b_neg")

    base_prob <- function(z,se) base_prob_interval(z,radius,se)
    #jacobian <- function(z,mean,var,se)prob_jacobian_interval(z,mean,var,se,radius=radius)
  }else{
    stop("Invalid testing type inputted. Valid forms of testing: `one_sided`, `interval`, and `two_sided`.")
  }
  alpha_m <- masking_params$alpha_m
  zeta <- masking_params$zeta
  lambda <- masking_params$lambda

  args <- list(testing = testing,
               model_type = model_type,
               rendpoint = rendpoint,
               lendpoint = lendpoint,
               alpha_m = alpha_m,
               zeta = zeta,
               lambda = lambda,
               niter_fit = niter_fit,
               niter_ms = niter_ms,
               nfits = nfits,
               masking_shape = masking_shape,
               p_to_z = p_to_z,
               z_to_p = z_to_p,
               beta_formula = beta_formula,
               nclasses = nclasses,
               all_a = all_a,
               n  = n,
               initialization = initialization,
               base_prob = base_prob
               )
  class(args) <- "args"
  return(args)
}

#' Create data class
#'
#' Class containing all relevant data inputs from user
#' @param x Vector of covariates
#' @param pvals Vector of pvals
#' @param z Vector of test statistics
#' @param se Vector of standard errors
#' @param args args class containing masking function arguments
#'
#' @return data class
#' @noRd
construct_data <- function(x,pvals,z,se,args){
  if(args$testing == "interval"){
    center <- (args$rendpoint + args$lendpoint)/2
    z <- z - center
  }
  # if(args$coerce_categorical){
  #   for(column in colnames(x)){
  #     num_unique <- length(unique(x[[column]]))
  #     if(num_unique/nrow(x)<0.10 & num_unique < 100){
  #       x[[column]] <- as.factor(x[[column]])
  #     }
  #   }
  # }

  x <- .scale_data(x)
  if(is.null(se)){
    se <- rep(1,nrow(x))
  }
  data <- list(x = x,
               pvals = pvals,
               z = z,
               se = se
               )
  class(data) <- "data"

  data <- data_preprocessing(data,args)

  return(data)
}




#' Initialize Parameters for Gaussian Mixture Model
#'
#' @param data data class
#' @param nclasses number of classes in Gaussian mixture model
#' @param initialization initialization type, kmeans, uniform, or random
#' @param testing testing type
#' @details Selects mu and variance by k-means
#'
#' @return params class containing beta, mu, var
#' @noRd
initialize_params <- function(data,nclasses,initialization,testing){

  mask <- data$mask
  a <- data$a

  true_z <- data$z[!mask]
  small_z <- data$small_z[mask]
  big_z <- data$big_z[mask]
  if(initialization == "random"){
    if(testing == "one_sided"){
      mu <- c(0, runif(nclasses-1, 1,4))
    }else{
      mu <- c(0, runif(nclasses-1, -4,4))
    }

    var <- c(0,runif(nclasses-1, 1, 6))
  }else if (initialization == "kmeans"){
    # Use true_z twice to count as double the weight
    all_z <- c(true_z,true_z,small_z,big_z)

    out <- kmeans(all_z, nclasses, nstart=50)

    mu <- as.numeric(out$centers)
    pred <- data.frame(z=all_z,class=out$cluster)
    var <- aggregate(pred$z,list(pred$class),"var") - 1

    colnames(var) <- c("group","value")
    var <- var[order(var$group),]
    var <- pmax(var$value,0)
    # If a class only has one observation, the empirical variance will be zero
    # Set NA values to 0
    var[is.na(var)] <- 0
  }else if(initialization == "uniform"){
    if(testing == "one_sided"){
      all_z <- c(true_z,true_z,small_z,big_z)
    }else{
      all_z <- c(true_z,true_z,small_z,big_z,-1*small_z,-1*big_z)
    }


    mu <- unlist(quantile(all_z, seq(0.2,1,length.out = nclasses+2)[2:(nclasses+1)]))
    mu <- unname(mu)
    var <- runif(nclasses,2,10)

  }else{
    stop("Unknown initialization scheme inputted.")
  }
  beta <- NULL


  params <- list(beta=beta, mu=mu, var=var)
  class(params) <- "params"
  return(params)
}


#' Scale data to [0,1] range
#'
#' @param x data frame
#'
#' @details Scales all features to [0,1] range for `multinom` function
#'
#' @return dataframe of x with scaled data
#' @noRd
.scale_data <- function(x){

  y <- data.frame(apply(x,2,function(z){
    if(is.numeric(z)){
      return((z-min(z))/(max(z)-min(z)))
    }else{
      return(z)
    }
  }) )
  colnames(y) <- colnames(x)
  return(y)
}

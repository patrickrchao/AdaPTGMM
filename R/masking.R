#' Data Preprocessing function
#'
#' Computes p_values and test statistics from dataset
#' @param data data class
#' @param args args class
#'
#' @return data class augmented with big/small test statistics, p_values, a_i, and boolean mask.
#' @noRd
data_preprocessing <- function(data,args){
  # set mask to true
  # preprocess p_values and test_statistics
  # preprocess masking

  # Initialize z and p_values if uninitialized
  if(is.null(data$z)){
    data$p_values <- pmax(pmin(data$p_values, 1 - 1e-20), 1e-20)
    data$z <- args$p_to_z(data$p_values)
  }else if(is.null(data$p_values)){
    data$p_values <- args$z_to_p(data$z)
  }
  # Clamp p-values
  data$p_values <- pmax(pmin(data$p_values, 1 - 1e-20), 1e-20)

  p_values <- data$p_values
  z <- data$z
  num_hypo <- length(p_values)

  lambda <- args$lambda
  alpha_m <- args$alpha_m
  zeta <- args$zeta

  data <- masking(data,args)

  return(data)
}

#' Compute big/small test statistics and p-values based on masking parameters
#'
#' @param data class, contains information for p-values and test statistics
#' @param args class, contains information on masking function and model constants
#'
#' @return data with computed big/small test statistics, p-values, initial mask, and a
#' @noRd
masking <- function(data,args){
  alpha_m <- args$alpha_m
  lambda <- args$lambda
  zeta <- args$zeta

  p_values <- data$p_values
  z <- data$z
  num_hypo <- length(p_values)

  # Determine if hypothesis corresponds to big or small
  a <- rep("NONE",num_hypo)
  a[p_values < alpha_m] <-  "s"
  a[(p_values > lambda) & (p_values < lambda + alpha_m / zeta)] <-  "b"

  mask <- (a != "NONE")

  small_p_values <- rep(NA,length(p_values))
  big_p_values <- rep(NA,length(p_values))

  small_p_values[a=="s"] <- p_values[a=="s"]
  big_p_values[a=="b"] <- p_values[a=="b"]



  if(args$masking_shape == "tent"){
    small_p_values[a=="b"] <- alpha_m + zeta * (lambda - p_values[(a == "b") ])
    big_p_values[a=="s"] <- (alpha_m - p_values[a=="s"]) / zeta + lambda
  }else if(args$masking_shape == "comb"){
    small_p_values[a=="b"] <- zeta * (p_values[a=="b"] - lambda)
    big_p_values[a=="s"] <- p_values[a=="s"] / zeta + lambda
  }else{
    stop("Invalid masking shape.")
  }

  small_z <- rep(NA, length(p_values))
  big_z <- rep(NA, length(p_values))
  small_z[a!="NONE"] <-  args$p_to_z(small_p_values[a!="NONE"])

  if(alpha_m == 0.5 & lambda == 0.5 & zeta == 1 & args$masking_shape == "tent"){
    big_z[a!="NONE"] <-  - small_z[a!="NONE"]
  }else{
    big_z[a!="NONE"] <-  args$p_to_z(big_p_values[a!="NONE"])
  }
  # Add to data class
  data$small_z <- small_z
  data$big_z <- big_z
  data$small_p_values <- small_p_values
  data$big_p_values <- big_p_values
  data$mask <- mask
  data$a <- a

  return(data)
}




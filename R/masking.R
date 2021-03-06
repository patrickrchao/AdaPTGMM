#' Select masking parameters
#'
#' zeta influences the minimum possible number of rejections. The minimum number of possible rejections
#' at FDR level \code{alpha} is \code{1/(zeta alpha)}. Thus for \code{alpha}=0.05, this corresponds to
#' \code{20/zeta}.
select_masking_params <- function(n,alpha_m,zeta,lambda,alpha_level=0.05){
    if(is.null(alpha_m) | is.null(zeta) | is.null(lambda)){
      warning("Masking parameter alpha_m, zeta, or lambda found to be NULL. Automatically selecting masking function. See documentation for details.")
      if(is.null(zeta)){
        zeta <- min(1 / alpha_level, max( 300 / (alpha_level * n), 2) )
      }
      alpha_m <- 0.9 / (zeta + 1)
      lambda <- alpha_m
    }
    masking_params <- list(alpha_m=alpha_m, zeta=zeta, lambda=lambda)
    return(masking_params)
}


#' Data Preprocessing function
#'
#' Computes pvals and test statistics from dataset
#' @param data data class
#' @param args args class
#' @param Boolean whether to resample blue p-values as uniform
#' @return data class augmented with big/small test statistics, pvals, a_i, and boolean mask.
#' @noRd
data_preprocessing <- function(data,args,randomize_pvals){
  # set mask to true
  # preprocess pvals and test_statistics
  # preprocess masking
  se <- data$se
  # Initialize z and pvals if uninitialized
  if(is.null(data$z)){
    data$pvals <- pmax(pmin(data$pvals, 1 - 1e-10), 1e-40)
    if(randomize_pvals){
      lambda <- args$lambda
      alpha_m <- args$alpha_m
      zeta <- args$zeta
      pvals <- data$pvals
      subset <- (pvals >= lambda) & (pvals <= (lambda + alpha_m*zeta))
      num <- sum(subset)
      pvals[subset] <- runif(num,lambda,lambda + alpha_m*zeta)
      data$pvals <- pvals
    }
    data$z <- args$p_to_z(data$pvals) * se
  }else if(is.null(data$pvals)){
    data$pvals <- args$z_to_p(data$z / se)
  }
  # Clamp p-values
  data$pvals <- pmax(pmin(data$pvals, 1 - 1e-10), 1e-40)

  pvals <- data$pvals
  z <- data$z
  num_hypo <- length(pvals)

  lambda <- args$lambda
  alpha_m <- args$alpha_m
  zeta <- args$zeta

  data <- masking(data,args)
  # if(!randomize_pvals){
  #   check_pval_dist(data$pvals,args)
  # }
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
  pvals <- data$pvals
  z <- data$z
  se <- data$se
  num_hypo <- length(pvals)

  # Determine if hypothesis corresponds to big or small
  a <- rep("NONE",num_hypo)
  a[pvals < alpha_m] <-  "s"
  a[(pvals > lambda) & (pvals < lambda + alpha_m * zeta)] <-  "b"

  mask <- (a != "NONE")

  small_pvals <- rep(NA,length(pvals))
  big_pvals <- rep(NA,length(pvals))

  small_pvals[a=="s"] <- pvals[a=="s"]
  big_pvals[a=="b"] <- pvals[a=="b"]



  if(args$masking_shape == "tent"){
    small_pvals[a=="b"] <- alpha_m + (lambda - pvals[(a == "b") ]) / zeta
    big_pvals[a=="s"] <- (alpha_m - pvals[a=="s"]) * zeta + lambda
  }else if(args$masking_shape == "comb"){
    small_pvals[a=="b"] <-  (pvals[a=="b"] - lambda) / zeta
    big_pvals[a=="s"] <- pvals[a=="s"] * zeta + lambda
  }else{
    stop("Invalid masking shape.")
  }
  big_pvals <- pmin(big_pvals,1-1e-12)

  small_z <- rep(NA, length(pvals))
  big_z <- rep(NA, length(pvals))

  small_z[a!="NONE"] <-  args$p_to_z(small_pvals[a!="NONE"]) * se[a!="NONE"]

  if(alpha_m == 0.5 & lambda == 0.5 & zeta == 1 & args$masking_shape == "tent" & args$testing == "one_sided"){
    big_z[a!="NONE"] <-  - small_z[a!="NONE"]
    # Reveal pvals between 0.45 and 0.55 for symmetric masking to mimic AdaPT
    mask[abs(pvals-0.5) < 0.05] <- FALSE
  }else{
    big_z[a!="NONE"] <-  args$p_to_z(big_pvals[a!="NONE"]) * se[a!="NONE"]
  }
  # Add to data class
  data$small_z <- small_z
  data$big_z <- big_z
  data$small_pvals <- small_pvals
  data$big_pvals <- big_pvals
  data$mask <- mask
  data$a <- a
  #data$mask[data$small_pvals > quantile(data$small_pvals,0.9,na.rm = T)] <- FALSE

  return(data)
}




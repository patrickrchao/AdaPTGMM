p_value_preprocessing <- function(data,args){

  lambda <- args$lambda
  alpha_m <- args$alpha_m
  zeta <- args$zeta
  p_values <- data$p_values
  if(!('z' %in% names(data))){
    data$z <- args$p_to_z(p_values)
  }
  # if p is between 0 and alpha_m (0.05) or
  # p is between lambda and lambda+alpha/zeta (0.4 to 0.9)
  # mask is TRUE
  a <- rep("NONE",length(data$p_values))
  a[p_values < alpha_m] <-  "s"
  a[(p_values > lambda) & (p_values < lambda + alpha_m / zeta)] <-  "b"
  mask <- (a != "NONE")
  observed_p_values <- p_values
  big_small <- masking_function(p_values,a,args,data$z)

  small_p_values <- big_small$small_p_values
  big_p_values <- big_small$big_p_values
  small_z <- big_small$small_z
  big_z <- big_small$big_z

  observed_p_values[mask] <- small_p_values[mask]

  data$a <- a
  data$small_p_values <- small_p_values
  data$big_p_values <- big_p_values
  data$observed_p_values <- observed_p_values
  data$mask <- mask
  data$small_z <- small_z
  data$big_z <- big_z

  data <- reveal(data)

  return(data)
}

masking <- function(data,args){

  data$mask <- (data$a != "NONE") & data$mask

  data <- reveal(data,args)
  return(data)
}

# Reveal z and p for unmasked data
reveal <- function(data,args){



  mask <- data$mask

  data$observed_p_values[!mask] <- data$p_values[!mask]
  data$small_p_values[!mask] <- NA
  data$big_p_values[!mask] <- NA
  data$small_z[!mask] <- NA
  data$big_z[!mask] <- NA

  return(data)
}


masking_function <- function(p_values,a,args, true_z){
  alpha_m <- args$alpha_m
  lambda <- args$lambda
  zeta <- args$zeta

  small_p_values <- rep(NA,length(p_values))
  big_p_values <- rep(NA,length(p_values))

  small_p_values[a=="s"] <- p_values[a=="s"]
  big_p_values[a=="b"] <- p_values[a=="b"]



  if(args$tent){
    small_p_values[a=="b"] <- alpha_m + zeta*(lambda - p_values[(a == "b") ])
    big_p_values[a=="s"] <- (alpha_m - p_values[a=="s"])/zeta + lambda

  }
  else{
    small_p_values[a=="b"] <- zeta * (p_values[a=="b"]-lambda)
    big_p_values[a=="s"] <- p_values[a=="s"]/zeta + lambda
  }

  small_z <- rep(NA,length(p_values))
  big_z <- rep(NA,length(p_values))
  small_z[a!="NONE"] <-  args$p_to_z(small_p_values[a!="NONE"])

  if(alpha_m==0.5 && lambda == 0.5 && args$tent){
    big_z[a!="NONE"] <-  - small_z[a!="NONE"]
  }else{
    big_z[a!="NONE"] <-  args$p_to_z(big_p_values[a!="NONE"])
  }
  # If exact values known, replace them
  if (!is.logical(true_z)){

    small_z[a=="s"] <- true_z[a=="s"]
    big_z[a=="b"] <- true_z[a=="b"]
  }
  return(list(small_p_values = small_p_values,big_p_values=big_p_values,small_z = small_z,big_z=big_z))
}

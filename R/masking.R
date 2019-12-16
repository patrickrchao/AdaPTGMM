masking <- function(data,params){#,alpha_m,lambda,zeta){
  mask <- data$mask
  p_i <- data$p_values
  alpha_m <- params$alpha_m
  lambda <- params$lambda
  zeta <- params$zeta
  # if p is between 0 and alpha_m (0.05) or
  # p is between lambda and lambda+alpha/zeta (0.4 to 0.9)
  # mask is TRUE

  mask <- mask &((p_i<alpha_m) | (p_i>lambda & p_i<lambda+alpha_m/zeta))
  if(params$tent){
    masked_p_i <- p_i*((1-mask)| (p_i<alpha_m)) +mask*(p_i>=alpha_m)*(alpha_m+zeta*(lambda-p_i))
  }
  else{
    masked_p_i <- p_i*((1-mask)| (p_i<alpha_m)) +mask*(p_i>=alpha_m)*(zeta*(p_i-lambda))
  }
  data$masked_p_i <- masked_p_i
  data$mask <- mask
  return(data)
}

inverse_masking <- function(data,params){
  mask <- data$mask
  masked_p_i <- data$masked_p_i
  alpha_m <- params$alpha_m
  lambda <- params$lambda
  zeta <- params$zeta
  small <- masked_p_i
  if(params$tent){
    big <- ifelse(mask,(-1/zeta*(masked_p_i-alpha_m)+lambda),(1-mask))
  }else{
    big <- ifelse(mask,(masked_p_i)/zeta+lambda,masked_p_i)
  }
  small_z <-  ifelse(mask,params$p_to_z(small),data$z)
  big_z <-  ifelse(mask,params$p_to_z(big),data$z)
 # small_z <-  ifelse(mask,-qnorm(small),data$z)
#  big_z <-  ifelse(mask,-qnorm(big),data$z)


  data$small <- small
  data$big <- big
  data$small_z <- small_z
  data$big_z <- big_z
  return(data)
}


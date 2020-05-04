
# Returns expectation of gamma given x and using current params of mu, variance
expectation_gamma <- function(data, params, args) {

  names <- paste0("gamma_",0:(args$num_classes-1))

  denominator <-0

  prob <- calculate_conditional_probabilities(data, params, args)
  output <-  data.frame(matrix(ncol = args$num_classes, nrow = length(prob$class_prob_0)))
  colnames(output) <- names
  denominator <- 0

  for (class in 0:(args$num_classes - 1)) {
    numerator <-  0
    for (a in args$all_a) {
      class_prob_str = paste0("class_prob_",class)
      weight_str <- paste0(a,"_",class)

      numerator <- numerator + prob[[class_prob_str]]*prob[[weight_str]]


    }
    output[paste0("gamma_",class)] <- numerator
    denominator <- denominator + numerator
  }

  gammas <- output/denominator

  if (sum(is.na(gammas))>0){
    browser()
  }else if(sum(colSums(gammas) == 0 ) > 0){
    browser()
  }
  # if(args$num_classes == 2){
  #   output <- data.frame(output$gamma_1)
  # }
  # if (args$testing_interval) {
  #   numerator <- prob$class_prob_1 * (prob$s_1 + prob$b_1 + prob$neg_s_1 + prob$neg_b_1)
  #   denominator <- numerator + prob$class_prob_0 * (prob$s_0 + prob$b_0 + prob$neg_s_0 + prob$neg_b_0)
  # } else{
  #   numerator <- prob$class_prob_1 * (prob$s_1 + prob$b_1)
  #   denominator <- numerator + prob$class_prob_0 * (prob$s_0 + prob$b_0)
  # }
  # expected_gamma <- numerator / denominator

  return(gammas)
}



# fit_beta_ridge <- function(x,gammas){
#
#   lambdas <- 10^seq(-8, -9, by = -.7)
#
#   glm_beta <- cv.glmnet(x[,-1], cbind(gammas,1), alpha = 0, family="binomial",lambda = lambdas,intercept=TRUE)
#
#   #glm_beta <- suppressWarnings(glm(formula,data=logistic_glm_data,family=binomial()))
#   beta <- as.matrix(coef(glm_beta))
#   row.names(beta)[1] <- "X0"
#   return(beta)
# }


# this function returns dataframe
# one dimensional case only contains
# w_{i,big}, w_{i,small} (k assumed to be one)
calculate_w <- function(data, params, args) {
  prob <- calculate_conditional_probabilities(data, params, args)

  small_z = data$small_z
  big_z = data$big_z


  output <- data.frame()
  names <- c("w_ika", "z", "k", "a", "i")
  num_data_points <- length(data$small_z)
  all_i = seq(1, num_data_points)

  denominator <-0

  total_rows <- num_data_points * length(args$all_a)*args$num_classes

  output <-  data.frame(matrix(ncol = length(names), nrow = total_rows))
  colnames(output) <- names

  index = 0

  for (a in args$all_a) {
    for (class in 0:(args$num_classes - 1)) {
      class_prob_str = paste0("class_prob_",class)
      weight_str <- paste0(a,"_",class)
      numerator <- prob[[class_prob_str]]*prob[[weight_str]]

      denominator <- denominator + numerator
      #browser()
      sign <- 1
      if(str_detect(a,"neg")){
        sign <- -1
      }
      if(str_detect(a,"s")){
        curr_z <- small_z
      }else{
        curr_z <- big_z
      }
      curr_z[!data$mask] <- data$z[!data$mask]
      output[(index*num_data_points+1):((index+1)*num_data_points),] <-
        c(numerator,sign*curr_z,rep(class,num_data_points),rep(a,num_data_points),all_i)
      index = index + 1
    }
  }

  denominator <- as.numeric(denominator)
  output$w_ika <- as.numeric(output$w_ika)
  output$i <- as.numeric(output$i)
  output$z <- as.numeric(output$z)
  output$k <- as.numeric(output$k)
  output$w_ika <- output$w_ika/as.numeric(denominator)
  return_var <- list(w_ika=output,denominator=denominator)
  return(return_var)
}

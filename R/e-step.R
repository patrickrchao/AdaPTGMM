library(glmnet)
library(pracma)

# Returns expectation of gamma given x and using current estimates of mu, variance
expectation_gamma <- function(data, est_params, params) {
  prob <- calculate_conditional_probabilities(data, est_params, params)
  if (params$testing_interval) {
    numerator <- prob$class_prob_1 * (prob$s_1 + prob$b_1 + prob$neg_s_1 + prob$neg_b_1)
    denominator <- numerator + prob$class_prob_0 * (prob$s_0 + prob$b_0 + prob$neg_s_0 + prob$neg_b_0)
  } else{
    numerator <- prob$class_prob_1 * (prob$s_1 + prob$b_1)
    denominator <- numerator + prob$class_prob_0 * (prob$s_0 + prob$b_0)
  }
  expected_gamma <- numerator / denominator
  return(expected_gamma)
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
calculate_w <- function(data, est_params, params) {
  prob <- calculate_conditional_probabilities(data, est_params, params)

  small_z = data$small_z
  big_z = data$big_z

  #prob_class_1 <- prob$expit_prob
  #prob_class_0 <- 1-prob_class_1
  output <- data.frame()

  names <- c("w_ika", "z", "k", "a", "i")

  num_data_points <- length(data$small_z)
  all_i = seq(1, num_data_points)

  # if (params$testing_interval) {
  #   # for(a in params$all_a){
  #   #   for(class in 0:(params$num_classes-1)){
  #   #
  #   #   }
  #   # }
  #
  #   denominator <-
  #     prob_class_1 * (prob$s_1 + prob$b_1 + prob$neg_s_1 + prob$neg_b_1) +
  #     prob_class_0 * (prob$s_0 + prob$b_0 + prob$neg_s_0 + prob$neg_b_0)
  #
  #   # w_{i,1,s}
  #   numerator <- prob_class_1 * prob$s_1
  #   w_1_s <- data.frame(numerator / denominator, small_z, 1, "s", all_i, numerator)
  #   colnames(w_1_s) <- names
  #
  #   # w_{i,1,b}
  #   numerator <- prob_class_1 * prob$b_1
  #   w_1_b <- data.frame(numerator / denominator, big_z, 1, "b", all_i, numerator)
  #   colnames(w_1_b) <- names
  #
  #   # w_{i,1,-s}
  #   numerator <- prob_class_1 * prob$neg_s_1
  #   w_1_neg_s <- data.frame(numerator / denominator, -1 * small_z, 1, "neg_s",all_i,numerator)
  #   colnames(w_1_neg_s) <- names
  #
  #   # w_{i,1,-b}
  #   numerator <- prob_class_1 * prob$neg_b_1
  #   w_1_neg_b <-
  #     data.frame(numerator / denominator, -1 * big_z, 1, "neg_b", all_i, numerator)
  #   colnames(w_1_neg_b) <- names
  #
  #   # w_{i,0,s}
  #   numerator <- prob_class_0 * prob$s_0
  #   w_0_s <- data.frame(numerator / denominator, small_z, 0, "s", all_i, numerator)
  #   colnames(w_0_s) <- names
  #
  #   # w_{i,b}
  #   numerator <- prob_class_0 * prob$b_0
  #   w_0_b <- data.frame(numerator / denominator, big_z, 0, "b", all_i, numerator)
  #   colnames(w_0_b) <- names
  #
  #   # w_{i,-s}
  #   numerator <- prob_class_0 * prob$neg_s_0
  #   w_0_neg_s <- data.frame(numerator / denominator, -1 * small_z, 0, "neg_s", all_i, numerator)
  #   colnames(w_0_neg_s) <- names
  #
  #   # w_{i,-b}
  #   numerator <- prob_class_0 * prob$neg_b_0
  #   w_0_neg_b <-data.frame(numerator / denominator, -1 * big_z, 0, "neg_b", all_i, numerator)
  #   colnames(w_0_neg_b) <- names
  #
  #   output <-
  #     rbind(w_1_s,
  #           w_1_b,
  #           w_1_neg_s,
  #           w_1_neg_b,
  #           w_0_s,
  #           w_0_b,
  #           w_0_neg_s,
  #           w_0_neg_b)
  # } else{

    denominator <-0
    # for (a in params$all_a) {
    #   for (class in 0:(params$num_classes - 1)) {
    #     class_prob_str = paste0("class_prob_",class)
    #     weight_str <- paste0(a,"_",class)
    #     denominator <- denominator + prob[class_prob_str]*weight_str
    #   }
    # }

    total_rows <- num_data_points * length(params$all_a)*params$num_classes
    # output <- data.frame(w_ika = numeric(),
    #                      z=numeric(),
    #                      k=numeric(),
    #                      a=character(),
    #                      i=numeric(),
    #                      numerator=numeric(),
    #                      stringsAsFactors = FALSE)
    output <-  data.frame(matrix(ncol = length(names), nrow = total_rows))
    colnames(output) <- names


    #output <- rbind(output,blank)



    index = 0

    for (a in params$all_a) {
      for (class in 0:(params$num_classes - 1)) {
        class_prob_str = paste0("class_prob_",class)
        weight_str <- paste0(a,"_",class)
        numerator <- prob[[class_prob_str]]*prob[[weight_str]]
        denominator <- denominator + numerator
        sign <- 1
        if(str_detect(a,"neg")){
          sign <- -1
        }
        if(str_detect(a,"s")){
          curr_z <- small_z
        }else{
          curr_z <- big_z
        }
        output[(index*num_data_points+1):((index+1)*num_data_points),] <-
          c(numerator,sign*curr_z,rep(class,num_data_points),rep(a,num_data_points),all_i)
        index = index + 1
      }
    }
    output$w_ika <- as.numeric(output$w_ika)
    output$i <- as.numeric(output$i)
    output$z <- as.numeric(output$z)
    output$k <- as.numeric(output$k)
    output$w_ika <- output$w_ika/denominator
    return_var <- list(w_ika=output,denominator=denominator)
    return(return_var)


    #denominator <- prob_class_1 * (prob$s_1 + prob$b_1) +
    #  prob_class_0 * (prob$s_0 + prob$b_0)

    # w_{i,1,s}
    # numerator <- prob_class_1 * prob$s_1
    # w_1_s <- data.frame(numerator / denominator, small_z, 1, "s", all_i, numerator)
    #
    # colnames(w_1_s) <- names
    #
    # # w_{i,1,b}
    # numerator <- prob_class_1 * prob$b_1
    # w_1_b <- data.frame(numerator / denominator, big_z, 1, "b", all_i, numerator)
    # colnames(w_1_b) <- names
    #
    # # w_{i,0,s}
    # numerator <- prob_class_0 * prob$s_0
    # w_0_s <- data.frame(numerator / denominator, small_z, 0, "s", all_i, numerator)
    # colnames(w_0_s) <- names
    #
    # # w_{i,0,b}
    # numerator <- prob_class_0 * prob$b_0
    # w_0_b <- data.frame(numerator / denominator, big_z, 0, "b", all_i, numerator)
    # colnames(w_0_b) <- names
    #
    # output <- rbind(w_1_s, w_1_b, w_0_s, w_0_b)
  # }
  return(output)
}

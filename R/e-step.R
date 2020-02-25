library(glmnet)
library(pracma)

# Returns expectation of gamma given x and using current estimates of mu, variance
expectation_gamma<- function(data,est_params,params){
  prob <- calculate_conditional_probabilities(data,est_params,params)
  if(params$testing_interval){
    numerator <- prob$expit_prob * (prob$s_1 + prob$b_1 + prob$neg_s_1 + prob$neg_b_1)
    denominator <- numerator + (1-prob$expit_prob) * (prob$s_0 + prob$b_0 + prob$neg_s_0 + prob$neg_b_0)
  }else{
    numerator <- prob$expit_prob * (prob$s_1 + prob$b_1)
    denominator <- numerator + (1-prob$expit_prob) * (prob$s_0 + prob$b_0)
  }
  expected_gamma <- numerator/denominator
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
calculate_w <- function(data,est_params,params){

  prob <- calculate_conditional_probabilities(data,est_params,params)

  small_z = data$small_z
  big_z = data$big_z

  prob_class_1 <- prob$expit_prob
  prob_class_0 <- 1-prob_class_1
  output <- data.frame()

  names <- c("w_ia","z","a")

  if(params$testing_interval){

    denominator <- prob_class_1 * (prob$s_1 + prob$b_1 + prob$neg_s_1 + prob$neg_b_1) +
      prob_class_0 * (prob$s_0 + prob$b_0 + prob$neg_s_0 + prob$neg_b_0)

    # w_{i,s}
    numerator <- prob_class_1 * prob$s_1
    w_s <- data.frame(numerator/denominator,small_z,"s")
    colnames(w_s) <- names

    # w_{i,b}
    numerator <- prob_class_1 * prob$b_1
    w_b <- data.frame(numerator/denominator,big_z,"b")
    colnames(w_b) <- names

    # w_{i,-s}
    numerator <- prob_class_1 * prob$neg_s_1
    w_neg_s <- data.frame(numerator/denominator,-1*small_z,"neg_s")
    colnames(w_neg_s) <- names

    # w_{i,-b}
    numerator <- prob_class_1 * prob$neg_b_1
    w_neg_b <- data.frame(numerator/denominator,-1*big_z,"neg_b")
    colnames(w_neg_b) <- names


    output <- rbind(w_s,w_b,w_neg_s,w_neg_b)
  }else{

    denominator <- prob_class_1 * (prob$s_1 + prob$b_1) +
      prob_class_0 * (prob$s_0 + prob$b_0)

    # w_{i,1,s}
    numerator <- prob_class_1 * prob$s_1
    w_s <- data.frame(numerator/denominator,small_z,"s")
    colnames(w_s) <- names

    # w_{i,1,b}
    numerator <- prob_class_1 * prob$b_1
    w_b <- data.frame(numerator/denominator,big_z,"b")
    colnames(w_b) <- names

    output <- rbind(w_s,w_b)
  }

  return(output)
}



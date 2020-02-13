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

  mu <- est_params$mu
  var <- est_params$var
  beta <- est_params$beta
  small_z = data$small_z
  big_z = data$big_z

  # P[p_small|gamma=1]
  #prob_small_1 <- gaussian_pdf(small_z,mu,var)/gaussian_pdf(small_z,0,1)
  # P[p_big|gamma=1]
  #prob_big_1 <- gaussian_pdf(big_z,mu,var)/gaussian_pdf(big_z,0,1)

  # P[z_small|gamma = 0]
  #prob_small_0 <- 1#gaussian_pdf(small_z,0,1)
  # P[z_big|gamma = 0]
  #prob_big_0 <- 1#gaussian_pdf(big_z,0,1)

  prob_class_1 <- prob$expit_prob
  prob_class_0 <- 1-prob_class_1
  output <- data.frame()

  names <- c("w_ia","z","a")

  if(params$testing_interval){

    denominator <- prob_class_1 * (prob$s_1 + prob$b_1 + prob$neg_s_1 + prob$neg_b_1) +
      prob_class_0 * (prob$s_0 + prob$b_0 + prob$neg_s_0 + prob$neg_b_0)

    # w_{i,s}
    #numerator <- prob_class_1*prob_small_1
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

    # w_{i,s}
    #numerator <- prob_class_1*prob_small_1
    numerator <- prob_class_1 * prob$s_1
    w_s <- data.frame(numerator/denominator,small_z,"s")
    colnames(w_s) <- names

    # w_{i,b}
    numerator <- prob_class_1 * prob$b_1
    w_b <- data.frame(numerator/denominator,big_z,"b")
    colnames(w_b) <- names

    # # w_{i,small}
    # #numerator <- prob_class_1*prob_small_1
    #
    # numerator <- prob_class_1 * prob$small_prob_alt
    #
    # #denominator <- prob_class_1*(prob_small_1+prob_big_1/params$zeta)+prob_class_0*(prob_small_0+prob_big_0/params$zeta)
    # denominator <- prob_class_1 * (prob$small_prob_alt + prob$big_prob_alt) + prob_class_0 * (prob$small_prob_null + prob$big_prob_null)
    # w_small <- data.frame(numerator/denominator,small_z,"small")#rep("small",nrow(data$full_x)))
    # colnames(w_small) <- names
    #
    # # w_{i,big}
    # numerator <- prob_class_1*prob$big_prob_alt#/params$zeta
    # w_big <- data.frame(numerator/denominator,big_z,"big")#,nrow(data$full_x)))
    #
    # colnames(w_big) <- names
    output <- rbind(w_s,w_b)
  }

  return(output)
}



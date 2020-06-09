#' Perform Maximization step for beta
#'
#' @param model model variable containing data, args, params
#' @param gammas Estimated gammas from e_step_gammas, dataframe where columns correspond to classes
#' and rows correspond to hypotheses. Each row sums to 1.
#'
#' @details Performs a multinomial logistic regression from nnet class
#' @return model
#' @noRd
m_step_beta <- function(model,gammas){

  x <- model$data$full_x
  ndf <- model$args$ndf
  nclasses <- model$args$nclasses

  colnames(gammas) <- 0:(nclasses-1)
  spread_gammas <- tidyr::gather(gammas,class,weight)
  gamma_class <- spread_gammas$class
  gamma_weight <- spread_gammas$weight
  multinom_data <- data.frame(x,gamma_class,gamma_weight)
  x_colnames <- colnames(multinom_data)[seq(ndf)]

  formula <- paste0("gamma_class ~ ",paste0(x_colnames,collapse = " + ")," -1")
  beta <- nnet::multinom(formula, multinom_data, weight = gamma_weight, trace = FALSE)

  # Update class probabilities
  model$data$class_prob <- class_prob(beta,nclasses)
  # Update beta
  model$params$beta <- beta


  return(model)
}



#' Perform Maximization step for mu and tau
#' @importFrom magrittr "%>%"
#' @noRd
m_step_mu_tau <- function(model,w_ika){
  args <- model$args
  params <- model$params
  data <- model$data

  z <- w_ika$z
  for (k in 1:(args$nclasses-1)){
    subset <- w_ika[w_ika$class == k,]
    params$mu[k+1] <- weighted_mean(z,subset$value)
    params$var[k+1] <- weighted_mean((z-params$mu[k+1])^2,subset$value)
    if(sum(is.na(c(params$mu,params$var)))>0){
      browser()
    }
  }


  return(params)
}



# levels(prob$a) <- 1:2
# prob$a <- as.numeric(prob$a)

# total_weight <- prob %>% subset(select=-c(a,class)) %>% rowSums()
# # for (class in 1:(args$num_classes-1)){
# #   w_ia <- output$w_ika %>% filter(k==class)
# #
# #   params$mu[class+1] <- weighted_mean(w_ia$w_ika,w_ia$z,w_ia)
# #   new_var <- max(weighted_mean(w_ia$w_ika,(w_ia$z-params$mu[class+1])^2,w_ia),0.0)
# #   params$var[class+1] <- new_var
# # }
# #prob %>% dplyr::group_by(class) %>% subset(select=-c(class,a))%>% dplyr::group_map(function(x) as.vector(t(x)))
#
# helper_mu <- function(row,data){
#
#   a <- row["a"]
#   values <- row[3:length(row)]
#   weighted_sum <- ifelse(a=="1",values*data$small_z,values*data$big_z)
#   return(weighted_sum)
#
# }
#
# helper_var <- function(row,data,mu){
#   a <- row["a"]
#   class <- as.numeric(row["class"])
#   values <- row[3:length(row)]
#   weighted_sum <- ifelse(a=="1",values*(data$small_z-mu[class+1])^2,values*(data$big_z-mu[class+1])^2)
#   return(weighted_sum)
# }
# # df_to_vec <- function(x) as.vector(t(x))
# # stacked_df <- do.call(rbind,
# #                 prob %>%
# #                   dplyr::group_by(class) %>%
# #                   subset(select=-c(a)) %>%
# #                   dplyr::group_map(~df_to_vec(.x)))
# # full_z <- c(data$small_z, data$big_z)
# # total_weights <- rowSums(stacked_df)
# # params$mu <- as.vector((stacked_df %*% full_z) / total_weights)
# # params$var <- as.vector((stacked_df %*% (full_z-params$mu)^2) / total_weights)
# prob_sum <- cbind(prob$class,total_weight,apply(prob,1,FUN=helper_mu,data))
# for(class in 1:)
# browser()


# library(dplyr)
# library(ggplot2)
# library(tidyr)
# num_df = 6
# data <- generate_data(1000,1,num_df)
# known <- data$known
# unknown <- data$unknown
# x <- known$x
# #spline_x <- known$spline_x
# p_values <- known$p_values
# #z <- known$z
#
# model <- create_model(x,p_values,num_df,iterations=75)
# data <- model$data
# params <- model$params
#
# data <- masking(data,params)
# data <- inverse_masking(data,params)
# plot(data$p_values,data$masked_p_i)
#
# temp_df <- data.frame(data[-3])
# test = temp_df %>% filter(p_values<params$alpha_m)%>% select(small,p_values) %>% mutate(diff = small-p_values)
# stopifnot(sum(abs(test$diff))<1e-8)
# #big small should be equal if mask is false
# test = temp_df %>% filter(mask==FALSE)%>% select(big,small,p_values) %>% mutate(diff = (big+small)/2-p_values)
# stopifnot(sum(abs(test$diff))<1e-8)
#
# # masked_p_values between lambda and alpha_m/zeta + lambda
# # should be equal to big
# test = temp_df %>% filter(p_values>params$lambda & p_values<params$alpha_m/params$zeta+params$lambda)%>% select(big,p_values) %>% mutate(diff = big-p_values)
# stopifnot(sum(abs(test$diff))<1e-8)
#
# #
# test = temp_df %>% mutate(diff1 = -qnorm(big)-big_z,diff2 = -qnorm(small)-small_z)
# stopifnot(sum(abs(test$diff1),abs(test$diff2))<1e-8)
#
#
# # Estimate using the correct parameters
# estimated_gamma <- expectation_gamma(data,unknown,params)
# ggplot(data.frame(estimated_gamma,unknown$gamma),aes(x=estimated_gamma,y=unknown.gamma))+geom_point(alpha=0.5)
# naive_estimate = 0
# print("Gamma Estimation Masked")
# print(paste("Correct estimated:",sum(round(estimated_gamma)==unknown$gamma)))
# print(paste("Naive Estimate:",sum(naive_estimate == unknown$gamma)))
#
# print(paste("RMSE with correct params:",round(sqrt(mean((estimated_gamma-unknown$gamma)^2)),3)))
# print(paste("Naive RMSE:",round(sqrt(mean((naive_estimate -unknown$gamma)^2)),3)))
#
# loess_fit <- loess(unknown$gamma ~ data$masked_p_i)
# plot(data$masked_p_i,estimated_gamma,xlim=c(0,0.05))
#
# res <- data.frame(data$masked_p_i,predict(loess_fit))
# colnames(res) <- c("masked","loess")
# res <- res[order(res$masked), ]
# lines(res$masked, res$loess, col = "blue")
# #plot(data$masked_p_i,unknown$gamma,xlim=c(0,0.1))
# beta_guess <- rep(0,num_df)
# mu_guess <- 2
# tau_guess <- 1
# est_params <- list(beta=beta_guess,mu=mu_guess,tau=tau_guess)
#
#
# data$mask <- FALSE
# data <- masking(data,params)
# data <- inverse_masking(data,params)
#
#
# print("")
# temp_df <- data.frame(data[-3])
#
# print("Gamma Estimation Unmasked")
#
# estimated_gamma <- expectation_gamma(data,unknown,params)
#
# ggplot(data.frame(estimated_gamma,unknown$gamma),aes(x=estimated_gamma,y=unknown.gamma))+geom_point(alpha=0.5)
# naive_estimate = 0
# print(paste("Correct estimated:",sum(round(estimated_gamma)==unknown$gamma)))
# print(paste("Naive Estimate:",sum(naive_estimate == unknown$gamma)))
#
# print(paste("RMSE with correct params:",round(sqrt(mean((estimated_gamma-unknown$gamma)^2)),3)))
# print(paste("Naive RMSE:",round(sqrt(mean((naive_estimate -unknown$gamma)^2)),3)))
#
#
# temp <- plot_fitting(data,params,unknown,title="Unmasked")
#
#
# data$mask <- TRUE
# data <- masking(data,params)
# data <- inverse_masking(data,params)
#
# print(paste("Percent of data masked:",round(sum(data$mask)/length(data$mask)*100,2)))
#
# temp <- plot_fitting(data,params,unknown,title="Masked")
#
#
#
#
#
# # est_params$beta <- est_params$beta*0
# # data$mask <-TRUE
# # data <- masking(data,params)
# #
# # data <- inverse_masking(data,params)
# # temp_df <- data.frame(data[-3])
# #
# # estimated_gamma <- expectation_gamma(data,unknown,params)
# # features <- data$full_x
# # response <- estimated_gamma
# # df <- data.frame(features,response)
# # testing <- glm(response ~ X0 + X1 + X2 + X3 + X4 + X5 -1, data =df, family = "binomial")
# # print(summary(testing))
# #
# #
# #
# # data$mask <-FALSE
# # data <- masking(data,params)
# #
# # data <- inverse_masking(data,params)
# #
# #
# # estimated_gamma <- expectation_gamma(data,unknown,params)
# # features <- data$full_x
# # response <- estimated_gamma
# # df <- data.frame(features,response)
# # testing <- glm(response ~ X0 + X1 + X2 + X3 + X4 + X5 -1, data =df, family = "binomial")
# # print(summary(testing))
# # print(unknown$beta)
#
#

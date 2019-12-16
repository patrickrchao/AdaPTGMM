# library(dplyr)
# library(ggplot2)
# library(tidyr)
# library(ggthemr)
# library("adaptMT")
# library(stringr)
# options(error =
#           function()
#           {
#             traceback(2);
#             if(!interactive()) quit('no', status = 1, runLast = FALSE)
#           })
#
# data(estrogen)
# p_values <- as.numeric(estrogen$pvals)
# x <-  as.numeric(estrogen$ord_high)
# x <- (x-mean(x))/sd(x)
# #ggplot(data.frame(x,p_values),aes(x,p_values))+geom_point(alpha=0.1)
#
# num_df = 8
# model <- create_model(x,p_values,num_df,iterations=10,alpha_m = 0.1,zeta = 0.2,lambda=0.4)
# data <- model$data
# params <- model$params
#
# data$mask <- TRUE
# data <- masking(data,params)
# data <- inverse_masking(data,params)
# print(paste("Percent of data masked:",round(sum(data$mask)/length(data$mask)*100,2)))
# plot_masking_function(data,params)
#
# beta_guess <- rep(0,num_df)
# mu_guess <- 2
# var_guess <- 1
# est_params <- list(beta=beta_guess,mu=mu_guess,var=var_guess)
#
# output <- AdaPTGMM(data,est_params,params)
# estrogen_gmm_log <- output$fdr_log
# est_params <- output$est_params
# rejections <- output$rejections
#
# plot(data$x,expit(data$full_x%*% est_params$beta),ylab  = "P(Gamma = 1)",main = "Gamma Probability for Covariates",xlab="Covariate")
# ggsave(paste0("Images/","Estrogen_Estimated_Spline",".png"),width=20,height=15,dpi=200,units="cm")
# print(paste("Estimated mu",est_params$mu))
# print(paste("Estimated var",est_params$var))
#
#
# df <- data.frame(data$x,data$p_values)
# for (curr_alpha in seq(0.9, 0.01, -0.01)) {
#   curr_rej <- rejections%>% filter(Alpha<curr_alpha+0.001, Alpha>curr_alpha-0.001)
#   df[as.logical(curr_rej[-1]),] %>% ggplot(aes(data.x,data.p_values)) + geom_point(alpha=0.4)+
#     labs(title = paste0("Rejected p values for alpha = ",curr_alpha))+
#     xlab("Covariate")+ylab("p value")+ylim(0,0.8)+xlim(-1.732,1.732)
#   filename <- str_pad(round((0.9-curr_alpha)*100),width=3,side="left",pad="0")
#   ggsave(paste0("Images/Rejection_Plots/",filename,".png"),width=12,height=8,dpi=200,units="cm")
# }
#
# # Adapt GLM
# pvals <- as.numeric(estrogen$pvals)
# adapt_x <- data.frame(x = as.numeric(estrogen$ord_high))
# # Define the exponential family for AdaPT (Section 4)
# dist <- beta_family()
#
# # Run adapt_glm
# library("splines")
# formulas <- paste0("ns(x, df = ", 5:10, ")")
#
#
# res <- adapt_glm(x = adapt_x, pvals = pvals, pi_formulas = formulas,
#                  mu_formulas = formulas, dist = dist, nfits = 10)
#
#
# # Plot the threshold curve and the level curves of local FDR
# #plot_1d_thresh(res, x,pvals,alpha = 0.1, "P-Value Thresholds")
# estrogen_adapt_log <- data.frame(res$nrejs,res$alphas,"AdaPT")
# colnames(estrogen_adapt_log) <- c("Rejected","FDPHat","Type")
#
#
#
#
#
# full_log <- rbind(estrogen_gmm_log[c("Rejected","FDPHat","Type")],estrogen_adapt_log)
#
# print(full_log %>% filter(FDPHat<0.301)%>%ggplot(aes(x=FDPHat,y=Rejected,fill=Type,color=Type))+geom_line(show.legend=TRUE) +
#         labs(title = "Rejections over FDP Hat for Estrogen")+
#         xlab("FDP Hat")+ylab("Number of Rejections"))
#
#
#
# ggsave(paste0("Images/","Estrogen_Rejections",".png"),width=20,height=15,dpi=800,units="cm")
#
# print(full_log %>% filter(FDPHat<0.101)%>%ggplot(aes(x=FDPHat,y=Rejected,fill=Type,color=Type))+geom_line(show.legend=TRUE) +
#         labs(title = "Rejections over FDP for Estrogen")+
#         xlab("FDP Hat")+ylab("Number of Rejections"))
#
#
#
# ggsave(paste0("Images/","Estrogen_Rejections_Zoom",".png"),width=20,height=15,dpi=800,units="cm")
#

# library(dplyr)
# library(ggplot2)
# library(tidyr)
# library(ggthemr)
#
# # require(devtools)
# # install_github("f1kidd/fmlogit")
# # library(fmlogit)
# #options(error=browser)
# options(error =
#           function()
#           {
#             traceback(2);
#             if(!interactive()) quit('no', status = 1, runLast = FALSE)
#           })
#
# num_df = 3
# num_classes = 3
# data <- generate_data(5000,num_df,num_classes = num_classes)
#
# known <- data$known
# unknown <- data$unknown
# print(paste("True args, mu: [",paste(unknown$mu, collapse = " "),"] var: [",paste(unknown$var, collapse=" "),"]",sep=""))
# print(paste("Percent of Null:",round(sum(unknown$theta<=0)/length(unknown$theta)*100,2)))
# x <- known$x
# #spline_x <- known$spline_x
# p_values <- known$p_values
# z <- known$z
#
# model <- create_model(x,p_values,num_df,iterations=200,alpha_m = 0.5,zeta = 1,lambda=0.5,tent=TRUE,num_classes=num_classes)
# #model <- create_model(x,p_values,num_df,iterations=200,alpha_m = 0.05,zeta = 1,lambda=0.5,tent=FALSE,num_classes=num_classes)
# #model <- create_model_interval(x,z,num_df,iterations=50,alpha_m = 0.05,zeta = 1,lambda=0.4,tent=FALSE,intervals=c(-1,1))
# data <- model$data
# args <- model$args
#
#
# data$mask <- TRUE
# data <- masking(data,args)
# data <- inverse_masking(data,args)
# #plot_x_p_value_masking(data,args)
#
# plot_masking_function(data,args,"AdaPTGMM_Masking_Function")
#
# print(paste("Percent of data masked:",round(sum(data$mask)/length(data$mask)*100,2)))
#
# likelihood(data,unknown,args,optimal_param=TRUE)
# temp <- plot_fitting(data,args,unknown,title="Masked")
#
# params <- initialize_params(num_classes,num_df)
#
# start_time <- Sys.time()
# output <- AdaPTGMM(data,params,args,calc_actual_FDP = TRUE,unknown)
# end_time = Sys.time()
# elapsed = end_time - start_time
# print(paste0("Total Time for GMM Adapt: ",elapsed))
# gmm_log <- output$fdr_log
# params <- output$params
# rejections <- output$rejections
# print(ggplot(gmm_log,aes(x=FDPHat,y=Rejected))+geom_line())
#
#
# # Adapt GLM
# library("adaptMT")
# adapt_x <- data.frame(x = as.numeric(data$x))
# pvals <- data$p_values
# # Define the exponential family for AdaPT (Section 4)
# dist <- beta_family()
#
# # Run adapt_glm
# library("splines")
# formulas <- paste0("ns(x, df = ", 5:10, ")")
#
# start_time <- Sys.time()
# res <- adapt_glm(x = adapt_x, pvals = pvals, pi_formulas = formulas,
#                  mu_formulas = formulas, dist = dist, nfits = 10)
#
# end_time = Sys.time()
# elapsed = end_time - start_time
# print(paste0("Total Time for GLM Adapt: ",elapsed))
#
# # Plot the threshold curve and the level curves of local FDR
# #plot_1d_thresh(res, x,pvals,alpha = 0.1, "P-Value Thresholds")
# adapt_fdr_log <- data.frame(res$nrejs,res$alphas,"AdaPT")
# colnames(adapt_fdr_log) <- c("Rejected","FDPHat","Type")
#
#
#
#
# # adapt_mask_model <- create_model(x,p_values,num_df,iterations=25,alpha_m = 0.5,zeta = 1,lambda=0.5,tent=TRUE)
# # data <- adapt_mask_model$data
# # args <- adapt_mask_model$args
# # data$mask <- TRUE
# # data <- masking(data,args)
# # data <- inverse_masking(data,args)
# # #plot_x_p_value_masking(data,args)
# #
# # plot_masking_function(data,args)
# #
# #
# # params <- initialize_params(num_classes,num_df)
# #
# #
# # output <- AdaPTGMM(data,params,args,calc_actual_FDP = TRUE,unknown)
# # adapt_mask_gmm_log <- output$fdr_log
# # adapt_mask_gmm_log$Type <- "AdaPTGMM AdaPT Mask"
# #plot_fitting(data,args,unknown,title="Masked")
#
# full_log <- rbind(gmm_log[c("Rejected","FDPHat","Type")],adapt_fdr_log)#,adapt_mask_gmm_log[c("Rejected","FDPHat","Type")])
# ggthemr('fresh')
# print(full_log %>% filter(FDPHat<0.301)%>%ggplot(aes(x=FDPHat,y=Rejected,fill=Type,color=Type))+geom_line(show.legend=TRUE) +
#         labs(title = "Rejections over FDP Hat for Gaussian")+
#         xlab("FDP Hat")+ylab("Number of Rejections"))
#
# ggsave(paste0("Images/","Gaussian_Rejections",".png"),width=20,height=15,dpi=200,units="cm")
#
# print(full_log %>% filter(FDPHat<0.101)%>%ggplot(aes(x=FDPHat,y=Rejected,fill=Type,color=Type))+geom_line(show.legend=TRUE) +
#         labs(title = "Rejections over FDP Hat for Gaussian")+
#         xlab("FDP Hat")+ylab("Number of Rejections"))
#
# ggsave(paste0("Images/","Gaussian_Rejections_Zoom",".png"),width=20,height=15,dpi=200,units="cm")
#
#

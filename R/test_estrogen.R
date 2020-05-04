library(dplyr)
library(ggplot2)
library(tidyr)
library(ggthemr)

# require(devtools)
# install_github("f1kidd/fmlogit")
# library(fmlogit)
#options(error=browser)
options(error =function(){traceback(2);if(!interactive()) quit('no', status = 1, runLast = FALSE)})


## Parameters

num_df = 6
num_classes = 2
interval = FALSE
iterations = 5
alpha_m = 0.1 #0.05
zeta = 0.2#0.1
lambda = 0.4 # 0.4
tent = FALSE

#Assert alpha_m/zeta + lambda < 1


data(estrogen)
p_values <- as.numeric(estrogen$pvals)
x <-  as.numeric(estrogen$ord_high)
x <- (x-mean(x))/sd(x)



output <- AdaPTGMM(x,p_values=p_values,num_df=num_df,iterations=iterations,alpha_m = alpha_m,zeta = zeta,lambda=lambda,tent=tent,num_classes=num_classes)
gmm_log <- output$fdr_log
params <- output$params
rejections <- output$rejections
print(ggplot(gmm_log,aes(x=FDPHat,y=Rejected))+geom_line())


# Adapt GLM
library("adaptMT")
adapt_x <- data.frame(x = as.numeric(x))
pvals <- p_values
# Define the exponential family for AdaPT (Section 4)
dist <- beta_family()

# Run adapt_glm
library("splines")
formulas <- paste0("ns(x, df = ", 5:10, ")")

res <- adapt_glm(x = adapt_x, pvals = pvals, pi_formulas = formulas,
                 mu_formulas = formulas, dist = dist, nfits = 10)

# Plot the threshold curve and the level curves of local FDR
#plot_1d_thresh(res, x,pvals,alpha = 0.1, "P-Value Thresholds")
adapt_fdr_log <- data.frame(res$nrejs,res$alphas,"AdaPT")
colnames(adapt_fdr_log) <- c("Rejected","FDPHat","Type")


# AdaPTGMM with AdaPT Masking
#adapt_model <- create_model(x,p_values,num_df,iterations=iterations,alpha_m = 0.49,zeta = 0.99,lambda=0.5,tent=TRUE,num_classes=num_classes)
#plot_masking_function(adapt_model$data,adapt_model$args)
output <- AdaPTGMM(x,p_values=p_values,num_df=num_df,iterations=iterations,alpha_m = 0.5,zeta = 1,lambda=0.5,tent=TRUE,num_classes=num_classes,)

# adapt_mask_model <- create_model(x,p_values,num_df,iterations=25,alpha_m = 0.5,zeta = 1,lambda=0.5,tent=TRUE)
# data <- adapt_mask_model$data
# args <- adapt_mask_model$args
# data$mask <- TRUE
# data <- masking(data,args)
# data <- inverse_masking(data,args)
# #plot_x_p_value_masking(data,args)
#
# plot_masking_function(data,args)
#
#
# params <- initialize_params(num_classes,num_df)
#
#
# output <- AdaPTGMM(data,params,args,calc_actual_FDP = TRUE,unknown)
adapt_mask_gmm_log <- output$fdr_log
adapt_mask_gmm_log$Type <- "AdaPTGMM AdaPT Mask"


full_log <- rbind(gmm_log[c("Rejected","FDPHat","Type")],adapt_fdr_log,adapt_mask_gmm_log[c("Rejected","FDPHat","Type")])

ggthemr('fresh')
print(full_log %>% filter(FDPHat<0.301)%>%ggplot(aes(x=FDPHat,y=Rejected,fill=Type,color=Type))+geom_line(show.legend=TRUE) +
        labs(title = "Rejections over FDP Hat for Estrogen")+
        xlab("FDP Hat")+ylab("Number of Rejections"))

ggsave(paste0("Images/","Estrogen_Rejections",".png"),width=20,height=15,dpi=200,units="cm")

print(full_log %>% filter(FDPHat<0.101)%>%ggplot(aes(x=FDPHat,y=Rejected,fill=Type,color=Type))+geom_line(show.legend=TRUE) +
        labs(title = "Rejections over FDP Hat for Estrogen")+
        xlab("FDP Hat")+ylab("Number of Rejections"))

ggsave(paste0("Images/","Estrogen_Rejections_Zoom",".png"),width=20,height=15,dpi=200,units="cm")






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
# args <- model$args
#
# data$mask <- TRUE
# data <- masking(data,args)
# data <- inverse_masking(data,args)
# print(paste("Percent of data masked:",round(sum(data$mask)/length(data$mask)*100,2)))
# plot_masking_function(data,args)
#
# beta_guess <- rep(0,num_df)
# mu_guess <- 2
# var_guess <- 1
# params <- list(beta=beta_guess,mu=mu_guess,var=var_guess)
#
# output <- AdaPTGMM(data,params,args)
# estrogen_gmm_log <- output$fdr_log
# params <- output$params
# rejections <- output$rejections
#
# plot(data$x,expit(data$full_x%*% params$beta),ylab  = "P(Gamma = 1)",main = "Gamma Probability for Covariates",xlab="Covariate")
# ggsave(paste0("Images/","Estrogen_Estimated_Spline",".png"),width=20,height=15,dpi=200,units="cm")
# print(paste("Estimated mu",params$mu))
# print(paste("Estimated var",params$var))
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

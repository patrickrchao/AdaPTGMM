# library(dplyr)
# library(ggplot2)
# library(tidyr)
# library(ggthemr)
# options("scipen"=100, "digits"=4)
# options(error =
#           function()
#           {
#             traceback(2);
#             if(!interactive()) quit('no', status = 1, runLast = FALSE)
#           })
#
# num_df = 10
#
# # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# # @@@@@@@@@@@@@@@@@ EXPERIMENT 1  @@@@@@@@@@@@@@@@@@
# # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# #Convergence of multiple initializations to single dataset (should converge to MLE, hopefully all similar)
#
# data <- generate_data(num_samples = 5000, num_dim = 1,num_df = num_df)
# known <- data$known
# unknown <- data$unknown
# print(paste("True args, mu:",unknown$mu," var:",unknown$var))
# num_null <- sum(abs(unknown$theta)<1)
# print(paste("Percent of Null:",round(num_null/length(unknown$theta)*100,2)))
# print(paste("Number of Null:",num_null))
# x <- known$x
# #spline_x <- known$spline_x
# z <- known$z
# #z <- known$z
#
# model <- create_model_interval(x,z,num_df,iterations=200,alpha_m = 0.05,zeta = 0.1,lambda=0.4,tent=FALSE,
#                                intervals=c(-1,1))
# data <- model$data
# args <- model$args
# x <- data$x
# z <- data$z
#
# data$mask <- TRUE
# data <- masking(data,args)
# data <- inverse_masking(data,args)
#
# plot_fitting(data,args,unknown,title="Masked")
#
#
# # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# # @@@@@@@@@@@@@@@@@ EXPERIMENT 2  @@@@@@@@@@@@@@@@@@
# # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# #Variance of MLE to multiple draws of data same parameters (should all be similar hopefully)
# # num_dim = 1
# # beta <-  matrix(sample(-2:2,num_dim*(num_df),replace=TRUE),ncol=1)
# # beta[1] <- -2
# # mu <- 3
# # var <- 4
# #
# # title <- "Masked"
# # true_beta <- beta
# # true_mu <- mu
# # true_var <- var
# #
# # beta_seq = data.frame()
# # mu_seq = data.frame()
# # var_seq = data.frame()
# #
# # #all_data$mask = all_data$p_values<0.2 | all_data$p_values>0.8
# # #all_data$masked_p_value = all_data$mask*(1-all_data$p_values)+(1-all_data$mask)*all_data$p_values
# # for(i in seq(8)){
# #   data <- generate_data(num_samples = 10000, num_dim = 1,num_df = num_df, beta=beta, mu = mu, tau = mu)
# #   known <- data$known
# #   unknown <- data$unknown
# #
# #   x <- known$x
# #   #spline_x <- known$spline_x
# #   z <- known$z
# #   #z <- known$z
# #
# #   model <- create_model_interval(x,z,num_df,iterations=200,alpha_m = 0.05,zeta = 0.1,lambda=0.4,tent=FALSE,
# #                                  intervals=c(-1,1))
# #   data <- model$data
# #   args <- model$args
# #   x <- data$x
# #   z <- data$z
# #
# #   data$mask <- TRUE
# #   data <- masking(data,args)
# #   data <- inverse_masking(data,args)
# #
# #   beta_guess = sample(0:0,args$num_df,replace=TRUE)#+true_beta
# #   mu_guess = sample(2,size=1) #true_mu#
# #   var_guess = sample(2,size=1) #true_var
# #   #mu_guess = true_mu#
# #   #var_guess = true_var
# #
# #   params <- list(beta=beta_guess,mu=mu_guess,var=var_guess)
# #   out = fit_parameters(data,params,args,beta_seq,mu_seq,var_seq)
# #   beta_seq = out$beta_seq
# #   mu_seq = out$mu_seq
# #   var_seq = out$var_seq
# # }
# # beta_seq$Trial <- as.factor(beta_seq$Trial)
# # mu_seq$Trial <- as.factor(mu_seq$Trial)
# # var_seq$Trial <- as.factor(var_seq$Trial)
# # x_names <- x_col_names(ncol(data$full_x))
# # if(length(unknown)>0){
# #   true_beta_df <- data.frame(t(true_beta))
# #   colnames(true_beta_df) <- x_names
# #   true_beta_df <- gather(true_beta_df,key="coord",value="Value",x_names)
# #   true_beta_df$coord_new <- factor(true_beta_df$coord,levels=x_names)
# # }
# #
# #
# # # Beta convergence plot
# # plot_beta_seq <- beta_seq %>% gather(key="coord",value="Value",x_names)
# # plot_beta_seq$coord_new <- factor(plot_beta_seq$coord,levels=x_names)
# # curr_title <-  paste0(title," Beta Coefficient Convergence")
# # if(is.numeric(true_beta)){
# #   print(plot_beta_seq %>% ggplot(aes(x=Iteration,y=Value,color=Trial))+geom_line()+
# #           theme(legend.position = "none")+facet_wrap(~coord_new)+
# #           geom_hline(data=true_beta_df,aes(yintercept = Value),linetype="dotted",size=1,color="black")+
# #           ggtitle(curr_title)+
# #           theme(axis.text=element_text(size=15),axis.title=element_text(size=16),plot.title = element_text(size=20)))
# # }else{
# #   print(plot_beta_seq %>% ggplot(aes(x=Iteration,y=Value,color=Trial))+geom_line()+theme(legend.position = "none")+facet_wrap(~coord_new)+ggtitle(curr_title)+
# #           theme(axis.text=element_text(size=15),axis.title=element_text(size=16),plot.title = element_text(size=20)))
# #
# # }
# # ggsave(paste0("Images/",gsub(" ","_",curr_title),".png"),width=20,height=15,dpi=200,units="cm")
# # beta_seq$Trial <- as.numeric(beta_seq$Trial)
# #
# # # Mu convergence plot
# #
# # plot_mu_seq <- mu_seq %>% gather(key="coord",value="Value","Mu")
# # curr_title <-  paste0(title," Mu Convergence")
# # if(length(unknown)>0){
# #   print(plot_mu_seq %>% ggplot(aes(x=Iteration,y=Value,color=Trial))+geom_line()+theme(legend.position = "none")+geom_hline(yintercept = true_mu,linetype="dotted",size=1,color="black")+ggtitle(curr_title)+
# #           theme(axis.text=element_text(size=15),axis.title=element_text(size=16),plot.title = element_text(size=20)))
# # }else{
# #   print(plot_beta_seq %>% ggplot(aes(x=Iteration,y=Value,color=Trial))+geom_line()+theme(legend.position = "none")+facet_wrap(~coord_new)+ggtitle(curr_title)+
# #           theme(axis.text=element_text(size=15),axis.title=element_text(size=16),plot.title = element_text(size=20)))
# # }
# # ggsave(paste0("Images/",gsub(" ","_",curr_title),".png"),width=20,height=15,dpi=200,units="cm")
# # # var convergence plot
# # plot_var_seq <- var_seq %>% gather(key="coord",value="Value","Var")
# # curr_title <-  paste0(title," Variance Convergence")
# # if(is.numeric(true_var)){
# #   print(plot_var_seq %>% ggplot(aes(x=Iteration,y=Value,color=Trial))+geom_line()+theme(legend.position = "none")+geom_hline(yintercept = true_var,linetype="dotted",size=1,color="black")+ggtitle(curr_title)+
# #           theme(axis.text=element_text(size=15),axis.title=element_text(size=16),plot.title = element_text(size=20)))
# # }
# # ggsave(paste0("Images/",gsub(" ","_",curr_title),".png"),width=20,height=15,dpi=200,units="cm")
# # #ggsave(paste0("Images/",curr_title,".png"),width=20,height=15,dpi=800,units="cm")
# #
# # actual <- expit(data$full_x%*% true_beta)
# # final_betas <-filter(beta_seq,Iteration ==args$iterations)[,c(-(ncol(beta_seq)-1),-ncol(beta_seq))]
# # predictions <- expit(data$full_x %*% t(data.matrix(final_betas)))
# # plot_splines <-data.frame(data$x,predictions,actual) %>% gather(key="Trial",value="gamma_mean",-data.x,-actual)
# #
# # curr_title <-  paste0(title," Spline Estimation")
# # print(plot_splines%>% ggplot(aes(x=data.x,y=gamma_mean,color=Trial),size=1)+labs(x="X")+
# #         geom_line()+geom_line(aes(x=data.x,y=actual),linetype="dotted",color="black",size=1)+
# #         ggtitle(curr_title)+
# #         theme(axis.text=element_text(size=15),axis.title=element_text(size=16),plot.title = element_text(size=20))
# # )
# # ggsave(paste0("Images/",gsub(" ","_",curr_title),".png"),width=20,height=15,dpi=200,units="cm")

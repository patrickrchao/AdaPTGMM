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

num_df = 3
num_classes = 4
num_data_points = 5000
interval = TRUE
iterations = 5
alpha_m = 0.1 #0.05
zeta = 0.2#0.1
lambda = 0.1 # 0.4
tent = FALSE

run_fitting_experiments = FALSE
#Assert alpha_m/zeta + lambda < 1


all_data <- generate_data(num_data_points,num_df,num_classes = num_classes,interval=interval)
known <- all_data$known
unknown <- all_data$unknown
x <- known$x
p_values <- known$p_values
z <- known$z

print(paste("True args, mu: [",paste(unknown$mu, collapse = " "),"] var: [",paste(unknown$var, collapse=" "),"]",sep=""))
if(interval){
  print(paste("Percent of Null:",round(sum(abs(unknown$theta)<=1)/length(unknown$theta)*100,2)))
}else{
  print(paste("Percent of Null:",round(sum(unknown$theta<=0)/length(unknown$theta)*100,2)))
}


# Print initial likelihood

# likelihood(true_model,optimal_param=TRUE)
# #likelihood(revealed_data,unknown,args,optimal_param=TRUE)
#
# if(run_fitting_experiments){
#   plot_fitting(data,args,true_model,title="Masked")
# }


output <- AdaPTGMM(x,p_values=p_values,z=z,interval=interval,num_df=num_df,iterations=iterations,alpha_m = alpha_m,zeta = zeta,lambda=lambda,tent=tent,num_classes=num_classes,calc_actual_FDP = TRUE,unknown=unknown)
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
output <- AdaPTGMM(x,p_values=p_values,num_df=num_df,iterations=iterations,alpha_m = 0.5,zeta = 1,lambda=0.5,tent=TRUE,num_classes=num_classes,calc_actual_FDP = TRUE,unknown=unknown)

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
        labs(title = "Rejections over FDP Hat for Gaussian")+
        xlab("FDP Hat")+ylab("Number of Rejections"))

ggsave(paste0("Images/","Gaussian_Rejections",".png"),width=20,height=15,dpi=200,units="cm")

print(full_log %>% filter(FDPHat<0.101)%>%ggplot(aes(x=FDPHat,y=Rejected,fill=Type,color=Type))+geom_line(show.legend=TRUE) +
        labs(title = "Rejections over FDP Hat for Gaussian")+
        xlab("FDP Hat")+ylab("Number of Rejections"))

ggsave(paste0("Images/","Gaussian_Rejections_Zoom",".png"),width=20,height=15,dpi=200,units="cm")



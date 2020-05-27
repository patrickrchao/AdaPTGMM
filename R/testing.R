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
num_classes = 2
num_data_points = 200
interval = FALSE
iterations = 5
alpha_m = 0.1 #0.05
zeta = 0.2
lambda = 0.4 # 0.4
tent = FALSE

run_initialization_experiments = FALSE
#Assert alpha_m/zeta + lambda < 1
#plot_masking_function(alpha_m,lambda,zeta,tent,title="Masking Function")
all_data <- generate_data(num_data_points,num_df= 3,num_classes = 2,interval=interval)
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
if(run_initialization_experiments){
  initialization_experiment(x,p_values=p_values,z=z,interval=interval,num_df=num_df,iterations=50,alpha_m = alpha_m,zeta = zeta,lambda=lambda,tent=tent,num_classes=num_classes,unknown=unknown)
}


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
adapt_fdr_log <- data.frame(res$nrejs,res$alphas,"AdaPT")
colnames(adapt_fdr_log) <- c("Rejected","FDPHat","Type")


# AdaPTGMM with AdaPT Masking
output <- AdaPTGMM(x,p_values=p_values,num_df=num_df,iterations=iterations,alpha_m = 0.5,zeta = 1,lambda=0.5,tent=TRUE,num_classes=num_classes,calc_actual_FDP = TRUE,unknown=unknown)

adapt_mask_gmm_log <- output$fdr_log

adapt_mask_gmm_log$Type <- "GMM Sym."

file_name <- "Gaussian_Rejections"
plot_title <- "Rejections over FDP Hat for"
if(interval){
  file_name <- paste(file_name,"Interval",sep="_")
  plot_title <- paste(plot_title,"Interval")
}else{
  plot_title <- paste(plot_title,"One Sided")
}
full_log <- rbind(gmm_log[c("Rejected","FDPHat","Type")],adapt_fdr_log,adapt_mask_gmm_log[c("Rejected","FDPHat","Type")])

ggthemr('fresh')
print(full_log %>% filter(FDPHat<0.301)%>%ggplot(aes(x=FDPHat,y=Rejected,fill=Type,color=Type))+geom_line(show.legend=TRUE) +
        labs(title = plot_title)+theme(text = element_text(size=20))+
        xlab("FDP Hat")+ylab("Number of Rejections"))

ggsave(paste0("Images/",file_name,".png"),width=20,height=15,dpi=200,units="cm")

print(full_log %>% filter(FDPHat<0.101)%>%ggplot(aes(x=FDPHat,y=Rejected,fill=Type,color=Type))+geom_line(show.legend=TRUE) +
        labs(title = plot_title)+theme(text = element_text(size=20))+
        xlab("FDP Hat")+ylab("Number of Rejections"))

ggsave(paste0("Images/",file_name,"_Zoom",".png"),width=20,height=15,dpi=200,units="cm")


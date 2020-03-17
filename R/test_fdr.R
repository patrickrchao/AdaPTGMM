library(dplyr)
library(ggplot2)
library(tidyr)
library(ggthemr)
options(error =
          function()
          {
            traceback(2);
            if(!interactive()) quit('no', status = 1, runLast = FALSE)
          })

num_df = 10

data <- generate_data(5000,num_df)
known <- data$known
unknown <- data$unknown
print(paste("True params, mu:",unknown$mu," var:",unknown$var))
print(paste("Percent of Null:",round(sum(unknown$theta<=0)/length(unknown$theta)*100,2)))
x <- known$x
#spline_x <- known$spline_x
p_values <- known$p_values
#z <- known$z

model <- create_model(x,p_values,num_df,iterations=50,alpha_m = 0.05,zeta = 0.1,lambda=0.4,tent=FALSE)
#model <- create_model(x,p_values,num_df,iterations=25,alpha_m = 0.5,zeta = 1,lambda=0.5,tent=TRUE)
data <- model$data
params <- model$params


data$mask <- TRUE
data <- masking(data,params)
data <- inverse_masking(data,params)
#plot_x_p_value_masking(data,params)

plot_masking_function(data,params,"AdaPTGMM_Masking_Function")

print(paste("Percent of data masked:",round(sum(data$mask)/length(data$mask)*100,2)))

likelihood(data,unknown,params,optimal_param=TRUE)
temp <- plot_fitting(data,params,unknown,title="Masked")



beta_guess <- rep(0,num_df)
mu_guess <- c(0,2)
var_guess <- c(1,1)
est_params <- list(beta=beta_guess,mu=mu_guess,var=var_guess)

start_time <- Sys.time()
output <- AdaPTGMM(data,est_params,params,calc_actual_FDP = TRUE,unknown)
end_time = Sys.time()
elapsed = end_time - start_time
print(paste0("Total Time for GMM Adapt: ",elapsed))
gmm_log <- output$fdr_log
est_params <- output$est_params
rejections <- output$rejections
print(ggplot(gmm_log,aes(x=FDPHat,y=Rejected))+geom_line())


# Adapt GLM
library("adaptMT")
adapt_x <- data.frame(x = as.numeric(data$x))
pvals <- data$p_values
# Define the exponential family for AdaPT (Section 4)
dist <- beta_family()

# Run adapt_glm
library("splines")
formulas <- paste0("ns(x, df = ", 5:10, ")")

start_time <- Sys.time()
res <- adapt_glm(x = adapt_x, pvals = pvals, pi_formulas = formulas,
                 mu_formulas = formulas, dist = dist, nfits = 10)

end_time = Sys.time()
elapsed = end_time - start_time
print(paste0("Total Time for GMM Adapt: ",elapsed))

# Plot the threshold curve and the level curves of local FDR
#plot_1d_thresh(res, x,pvals,alpha = 0.1, "P-Value Thresholds")
adapt_fdr_log <- data.frame(res$nrejs,res$alphas,"AdaPT")
colnames(adapt_fdr_log) <- c("Rejected","FDPHat","Type")




adapt_mask_model <- create_model(x,p_values,num_df,iterations=25,alpha_m = 0.5,zeta = 1,lambda=0.5,tent=TRUE)
data <- adapt_mask_model$data
params <- adapt_mask_model$params
data$mask <- TRUE
data <- masking(data,params)
data <- inverse_masking(data,params)
#plot_x_p_value_masking(data,params)

plot_masking_function(data,params)


beta_guess <- rep(0,num_df)
mu_guess <- c(0,2)
var_guess <- c(1,1)
est_params <- list(beta=beta_guess,mu=mu_guess,var=var_guess)


output <- AdaPTGMM(data,est_params,params,calc_actual_FDP = TRUE,unknown)
adapt_mask_gmm_log <- output$fdr_log
adapt_mask_gmm_log$Type <- "AdaPTGMM AdaPT Mask"
#plot_fitting(data,params,unknown,title="Masked")

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




fit_parameters <- function(data,est_params,params,beta_seq=data.frame(),mu_seq=data.frame(),var_seq=data.frame()){

  if("Trial" %in% colnames(beta_seq)){
    trial = max(beta_seq$Trial)+1
  }else{
    trial = 0
  }
  x_names <- paste("X",0:(ncol(data$full_x)-1),sep="")

  n_dim = length(est_params$beta)
  beta_seq <- rbind(beta_seq,c(est_params$beta,0,trial))
  mu_seq <- rbind(mu_seq,c(est_params$mu,0,trial))
  var_seq <- rbind(var_seq,c(est_params$var,0,trial))

  for(i in seq(params$iterations)){

    est_params <- update_parameters(data,est_params,params)

    beta_seq <- rbind(beta_seq,c(est_params$beta,i,trial))
    mu_seq <- rbind(mu_seq,c(est_params$mu,i,trial))
    var_seq <- rbind(var_seq,c(est_params$var,i,trial))

  }

  colnames(beta_seq) <- c(x_names,"Iteration","Trial")
  colnames(mu_seq) <- c("Mu","Iteration","Trial")
  colnames(var_seq) <- c("Var","Iteration","Trial")
  out <- list()
  out$best_beta = est_params$beta
  out$beta_seq = beta_seq

  out$best_mu = est_params$mu
  out$mu_seq = mu_seq

  out$best_var = est_params$var
  out$var_seq = var_seq
  return(out)
}


plot_fitting<- function(data,params,unknown=TRUE,num_trials=8,title){

  if(length(unknown)>0){
    true_beta <- unknown$beta
    true_mu <- unknown$mu
    true_var <- unknown$var
  }
  beta_seq = data.frame()
  mu_seq = data.frame()
  var_seq = data.frame()

  #all_data$mask = all_data$p_values<0.2 | all_data$p_values>0.8
  #all_data$masked_p_value = all_data$mask*(1-all_data$p_values)+(1-all_data$mask)*all_data$p_values
  for(i in seq(num_trials)){
    beta_guess = sample(-1:1,params$num_df,replace=TRUE)#+true_beta
    beta_guess[1] = sample(-4:-1,1)
    #mu_guess = sample(2:5,size=1) #true_mu#
    #var_guess = sample(1:3,size=1) #true_var
    mu_guess = true_mu#
    var_guess = true_var
    est_params <- list(beta=beta_guess,mu=mu_guess,var=var_guess)
    out = fit_parameters(data,est_params,params,beta_seq,mu_seq,var_seq)
    beta_seq = out$beta_seq
    mu_seq = out$mu_seq
    var_seq = out$var_seq
  }

  beta_seq$Trial <- as.factor(beta_seq$Trial)
  mu_seq$Trial <- as.factor(mu_seq$Trial)
  var_seq$Trial <- as.factor(var_seq$Trial)
  x_names <- x_col_names(ncol(data$full_x))

  if(length(unknown)>0){
    true_beta_df <- data.frame(t(true_beta))
    colnames(true_beta_df) <- x_names
    true_beta_df <- gather(true_beta_df,key="coord",value="Value",x_names)
    true_beta_df$coord_new <- factor(true_beta_df$coord,levels=x_names)
  }


  # Beta convergence plot
  plot_beta_seq <- beta_seq %>% gather(key="coord",value="Value",x_names)
  plot_beta_seq$coord_new <- factor(plot_beta_seq$coord,levels=x_names)
  curr_title <-  paste0(title," Beta Coefficient Convergence")
  if(is.numeric(true_beta)){
    print(plot_beta_seq %>% ggplot(aes(x=Iteration,y=Value,color=Trial))+geom_line()+
            theme(legend.position = "none")+facet_wrap(~coord_new)+
            geom_hline(data=true_beta_df,aes(yintercept = Value),linetype="dotted",size=1,color="black")+
            ggtitle(curr_title)+
            theme(axis.text=element_text(size=15),axis.title=element_text(size=16),plot.title = element_text(size=20)))
  }else{
    print(plot_beta_seq %>% ggplot(aes(x=Iteration,y=Value,color=Trial))+geom_line()+theme(legend.position = "none")+facet_wrap(~coord_new)+ggtitle(curr_title)+
            theme(axis.text=element_text(size=15),axis.title=element_text(size=16),plot.title = element_text(size=20)))

  }
  ggsave(paste0("Images/",gsub(" ","_",curr_title),".png"),width=20,height=15,dpi=200,units="cm")
  beta_seq$Trial <- as.numeric(beta_seq$Trial)

  # Mu convergence plot

  plot_mu_seq <- mu_seq %>% gather(key="coord",value="Value","Mu")
  curr_title <-  paste0(title," Mu Convergence")
  if(length(unknown)>0){
    print(plot_mu_seq %>% ggplot(aes(x=Iteration,y=Value,color=Trial))+geom_line()+theme(legend.position = "none")+geom_hline(yintercept = true_mu,linetype="dotted",size=1,color="black")+ggtitle(curr_title)+
            theme(axis.text=element_text(size=15),axis.title=element_text(size=16),plot.title = element_text(size=20)))
  }else{
    print(plot_beta_seq %>% ggplot(aes(x=Iteration,y=Value,color=Trial))+geom_line()+theme(legend.position = "none")+facet_wrap(~coord_new)+ggtitle(curr_title)+
            theme(axis.text=element_text(size=15),axis.title=element_text(size=16),plot.title = element_text(size=20)))
  }
  ggsave(paste0("Images/",gsub(" ","_",curr_title),".png"),width=20,height=15,dpi=200,units="cm")
  # var convergence plot
  plot_var_seq <- var_seq %>% gather(key="coord",value="Value","Var")
  curr_title <-  paste0(title," Variance Convergence")
  if(is.numeric(true_var)){
    print(plot_var_seq %>% ggplot(aes(x=Iteration,y=Value,color=Trial))+geom_line()+theme(legend.position = "none")+geom_hline(yintercept = true_var,linetype="dotted",size=1,color="black")+ggtitle(curr_title)+
            theme(axis.text=element_text(size=15),axis.title=element_text(size=16),plot.title = element_text(size=20)))
  }
  ggsave(paste0("Images/",gsub(" ","_",curr_title),".png"),width=20,height=15,dpi=200,units="cm")
  #ggsave(paste0("Images/",curr_title,".png"),width=20,height=15,dpi=800,units="cm")

  actual <- expit(data$full_x%*% true_beta)
  final_betas <-filter(beta_seq,Iteration ==params$iterations)[,c(-(ncol(beta_seq)-1),-ncol(beta_seq))]
  predictions <- expit(data$full_x %*% t(data.matrix(final_betas)))
  plot_splines <-data.frame(data$x,predictions,actual) %>% gather(key="Trial",value="gamma_mean",-data.x,-actual)

  curr_title <-  paste0(title," Spline Estimation")
  print(plot_splines%>% ggplot(aes(x=data.x,y=gamma_mean,color=Trial),size=1)+labs(x="X")+
          geom_line()+geom_line(aes(x=data.x,y=actual),linetype="dotted",color="black",size=1)+
          ggtitle(curr_title)+
          theme(axis.text=element_text(size=15),axis.title=element_text(size=16),plot.title = element_text(size=20))
  )
  ggsave(paste0("Images/",gsub(" ","_",curr_title),".png"),width=20,height=15,dpi=200,units="cm")
  #browser()
  output <- list()
  output$beta_seq <- beta_seq
  output$mu_seq <- mu_seq
  output$var_seq <- var_seq
  return(beta_seq)
}


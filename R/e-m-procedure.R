
fit_parameters <- function(data,est_params,params,beta_seq=data.frame(),mu_seq=data.frame(),var_seq=data.frame()){

  if("Trial" %in% colnames(beta_seq)){
    trial = max(beta_seq$Trial)+1
  }else{
    trial = 0
  }
  beta_names <- coefficient_names("Beta",ncol(data$full_x))
  mu_names <- coefficient_names("Mu",params$num_classes)
  var_names <- coefficient_names("Var",params$num_classes)

  beta_seq <- rbind(beta_seq,c(est_params$beta,0,trial))
  mu_seq <-   rbind(mu_seq,  c(est_params$mu,  0,trial))
  var_seq <-  rbind(var_seq, c(est_params$var, 0,trial))
  # can speed this up
  for(i in seq(params$iterations)){

    est_params <- update_parameters(data,est_params,params)

    beta_seq <- rbind(beta_seq,c(est_params$beta,i,trial))
    mu_seq <-   rbind(mu_seq,  c(est_params$mu,  i,trial))
    var_seq <-  rbind(var_seq, c(est_params$var, i,trial))

  }

  colnames(beta_seq) <- c(beta_names,"Iteration","Trial")
  colnames(mu_seq) <- c(mu_names,"Iteration","Trial")
  colnames(var_seq) <- c(var_names,"Iteration","Trial")
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

  for(i in seq(num_trials)){
    beta_guess = sample(-1:1,params$num_df,replace=TRUE)#+true_beta
    beta_guess[1] = sample(-4:-1,1)
    #mu_guess = sample(2:5,size=1)
    #var_guess = sample(1:3,size=1)
    mu_guess <-  c(0,sample(3:5,size=params$num_classes-1,replace=FALSE))
    var_guess <-  c(1,sample(3:8,size=params$num_classes-1,replace=TRUE))
    #mu_guess = true_mu#
    #var_guess = true_var

    est_params <- list(beta=beta_guess,mu=mu_guess,var=var_guess)
    out = fit_parameters(data,est_params,params,beta_seq,mu_seq,var_seq)
    beta_seq = out$beta_seq
    mu_seq = out$mu_seq
    var_seq = out$var_seq
  }

  beta_seq$Trial <- as.factor(beta_seq$Trial)
  mu_seq$Trial <- as.factor(mu_seq$Trial)
  var_seq$Trial <- as.factor(var_seq$Trial)

  plot_variable_convergence(beta_seq,paste(title,"Beta Coefficient Convergence"),true_beta)
  plot_variable_convergence(mu_seq,paste(title,"Mu Coefficient Convergence"),true_mu,ignore_first_coef = TRUE)
  plot_variable_convergence(var_seq,paste(title,"Var Coefficient Convergence"),true_var,ignore_first_coef = TRUE)

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

  output <- list()
  output$beta_seq <- beta_seq
  output$mu_seq <- mu_seq
  output$var_seq <- var_seq
  return(beta_seq)
}

plot_variable_convergence <- function(df,title,true_value,ignore_first_coef = FALSE){

  if (ignore_first_coef){
    df <- df[-1]
    if(is.numeric(true_value) || is.data.frame(true_value)){
      true_value <- true_value[-1]
    }
  }
  col_names <- colnames(df)[1:(length(colnames(df))-2)]
  plot_df <- df %>% gather(key="coord",value="Value",col_names)
  if(length(col_names)>1){
    plot_df$coord <- factor(plot_df$coord,levels=col_names)
  }
  #curr_title <-  paste0(title," Beta Coefficient Convergence")
  if(is.numeric(true_value) || is.data.frame(true_value)){
    true_value_df<- data.frame(t(true_value))
    colnames(true_value_df) <- col_names
    true_value_df <- gather(true_value_df,key="coord",value="Value",col_names)
    true_value_df$coord <- factor(true_value_df$coord,levels=col_names)

    plot_obj <- plot_df %>% ggplot(aes(x=Iteration,y=Value,color=Trial))+geom_line()+
            theme(legend.position = "none")+facet_wrap(~coord)+ggtitle(title)+
            geom_hline(data=true_value_df,aes(yintercept = Value),linetype="dotted",size=1,color="black")+
            theme(axis.text=element_text(size=15),axis.title=element_text(size=16),plot.title = element_text(size=20))
  }else{
    plot_obj <- plot_df %>% ggplot(aes(x=Iteration,y=Value,color=Trial))+geom_line()+theme(legend.position = "none")+facet_wrap(~coord)+ggtitle(title)+
            theme(axis.text=element_text(size=15),axis.title=element_text(size=16),plot.title = element_text(size=20))

  }
  if (length(col_names)==1){
    plot_obj <- plot_obj +
      theme(
        strip.text.y = element_blank(),
        strip.text.x = element_blank()
      )
  }
  print(plot_obj)
  ggsave(paste0("Images/",gsub(" ","_",title),".png"),width=20,height=15,dpi=200,units="cm")
}



initialization_experiment <- function(x,
                                      interval=FALSE,
                                      p_values=FALSE,
                                      z=FALSE,
                                      num_df=10,
                                      alpha_m=0.05,
                                      zeta=0.1,
                                      lambda=0.4,
                                      spline=TRUE,
                                      iterations=20,
                                      tent=FALSE,
                                      num_classes = 2,
                                      unknown=FALSE,num_trials=7,title){


  if(interval){
    model <- create_model_interval(x,z,num_df,iterations=iterations,alpha_m = alpha_m,zeta = zeta,lambda=lambda,tent=tent,num_classes=num_classes,intervals=c(-1,1))
  }else{
    model <- create_model(x,p_values,num_df,iterations=iterations,alpha_m = alpha_m,zeta = zeta,lambda=lambda,tent=tent,num_classes=num_classes)
    true_args <- model$args
    true_args$num_df <- ncol(unknown$beta)
    true_args$num_classes <- nrow(unknown$beta)
    true_data <- model$data
    true_data$full_x <- unknown$full_x
    true_model <-  create_model_from_args(true_data,true_args,unknown)
  }
  models <- vector(mode="list", length=num_trials)
  names(models) <- 1:num_trials
  for(trial in 1:num_trials){
    if(interval){
      model <- create_model_interval(x,z,num_df,iterations=iterations,alpha_m = alpha_m,zeta = zeta,lambda=lambda,tent=tent,num_classes=num_classes,intervals=c(-1,1))
    }else{
      model <- create_model(x,p_values,num_df,iterations=iterations,alpha_m = alpha_m,zeta = zeta,lambda=lambda,tent=tent,num_classes=num_classes)

    }
    model <- fit_parameters(model)
    models[[trial]] <- model

  }

  plot_theta_density(models,true_model)
  plot_theta_density(models,true_model)
  plot_theta_density(models,true_model)
  plot_theta_density(models,true_model)
  plot_theta_density(models,true_model)
}

plot_theta_density_from_params <- function(model,unknown){
  true_args <- model$args
  true_args$num_df <- ncol(unknown$beta)
  true_args$num_classes <- nrow(unknown$beta)
  true_data <- model$data
  true_data$full_x <- unknown$full_x
  true_model <-  create_model_from_args(true_data,true_args,unknown)
  models <- list()
  models$"1" <- model
  plot_theta_density(models,true_model)
}
old_initialization_experiment <- function(data,args,true_model=FALSE,num_trials=8,title){

  if(!is.logical(true_model)){
    true_params <- true_model$params
    true_beta <- true_args$beta
    true_mu <- true_args$mu
    true_var <- true_args$var
  }

  beta_seq = data.frame()
  mu_seq = data.frame()
  var_seq = data.frame()

  for(i in seq(num_trials)){
    print("New Iteration")

    model <- create_model_from_args(data,args)
    print("Initial Likelihood")
    likelihood(model)
    out = fit_parameters_logging(model,beta_seq,mu_seq,var_seq)
    model <- out$model
    beta_seq = out$beta_seq
    mu_seq = out$mu_seq
    var_seq = out$var_seq
    params <- model$params
    likelihood(model)

    plot_density(model,true_model)
    plot_density(model,true_model)

  }

  beta_seq$Trial <- as.factor(beta_seq$Trial)
  mu_seq$Trial <- as.factor(mu_seq$Trial)
  var_seq$Trial <- as.factor(var_seq$Trial)

  plot_variable_convergence(beta_seq,paste(title,"Beta Coefficient Convergence"),true_beta)
  plot_variable_convergence(mu_seq,paste(title,"Mu Coefficient Convergence"),true_mu,ignore_first_coef = TRUE)
  plot_variable_convergence(var_seq,paste(title,"Var Coefficient Convergence"),true_var,ignore_first_coef = TRUE)

  final_betas <-filter(beta_seq,Iteration ==args$iterations)
  actual <- probability_from_spline(data$full_x,true_beta,args$num_classes)
  colnames(actual) <- paste0("Class_",0:(args$num_classes-1))

  trial_predictions <- list()
  for(curr_trial in 0:(num_trials-1)){
    prob <- probability_from_spline(data$full_x,t(subset(final_betas,Trial==curr_trial,select=-c(Class,Iteration,Trial))),args$num_classes)
    trial_predictions[toString(curr_trial)] <- data.frame(rowSums(prob))
  }

  trial_predictions <- data.frame(trial_predictions)
  colnames(trial_predictions) <- 0:(num_trials-1)
  for (curr_class in 0:(args$num_classes-1)){
    curr_beta <- subset(final_betas,Class==curr_class,select=-c(Class,Iteration,Trial))
    curr_actual <- actual[,paste0("Class_",curr_class)]

    predictions <- probability_from_spline(data$full_x, t(curr_beta),args$num_classes)/trial_predictions
    colnames(predictions) <- 0:(num_trials-1)
    plot_splines <- data.frame(data$x,predictions,curr_actual)
    colnames(plot_splines) <-  c("X",paste0("Trial_",0:(num_trials-1)),"Actual")
    plot_splines <- plot_splines %>% gather(key="Trial",value="Class_Probability",-X,-Actual)
    curr_title <-  paste(title,"Spline Estimation Class",curr_class)
    print(plot_splines%>% ggplot(aes(x=X,y=Class_Probability,color=Trial),size=1)+labs(x="X")+
            geom_line()+geom_line(aes(x=X,y=Actual),linetype='solid',color="black",size=1.5)+
            ggtitle(curr_title)+
            theme(axis.text=element_text(size=15),axis.title=element_text(size=16),plot.title = element_text(size=20))
    )
    ggsave(paste0("Images/",gsub(" ","_",curr_title),".png"),width=20,height=15,dpi=200,units="cm")
  }
  #temp <- split(beta_seq,beta_seq$Trial)
  # #beta_seq$Trial <- factor(beta_seq$Trial)
  # final_betas <- filter(beta_seq,Iteration ==args$iterations)
  # predictions <- data.frame()
  # for (curr_trial in 0:(num_trials-1)){
  #   curr_beta <- subset(final_betas,Trial==curr_trial,select=-c(Class,Trial,Iteration))
  #   curr_predictions <- data.frame(exp(data$full_x %*% t(curr_beta)))
  #   curr_predictions$Trial <- curr_trial
  #   predictions <- rbind(predictions,curr_predictions)
  # }
  # colnames(predictions) <- paste0("Class_",0:(args$num_classes-1))
  # # temp <- filter(beta_seq,Iteration ==args$iterations) %>% group_by("Trial",add = FALSE)
  # # temp <- temp[,1:5]
  # # final <- temp %>% do(probability_from_spline(.,data$full_x,args$num_classes))
  # actual <- expit(data$full_x%*% true_beta)
  # colnames(actual) <- paste0("True_Class_",0:(args$num_classes-1))
  # #final_betas <-filter(beta_seq,Iteration ==args$iterations)[,c(-(ncol(beta_seq)-1),-ncol(beta_seq))]
  # #predictions <- expit(data$full_x %*% t(data.matrix(final_betas)))
  # # plot_splines <-data.frame(data$x,predictions,actual) %>% gather(key="Trial",value="gamma_mean",-data.x,-actual)
  # predicted_splines <- data.frame(data$x,predictions)%>% gather(key="Class",value="Class_Prob",-data.x)
  # true_splines <- data.frame(actual) %>% gather(key="Class",value="True_Class_Prob")
  # plot_splines <- cbind(predicted_splines,true_splines)
  # curr_title <-  paste0(title," Spline Estimation")
  # print(plot_splines%>% ggplot(aes(x=data.x,y=gamma_mean,color=Trial),size=1)+labs(x="X")+facet_wrap(~class)+
  #         geom_line()+geom_line(aes(x=data.x,y=actual),linetype="dotted",color="black",size=1)+
  #         ggtitle(curr_title)+
  #         theme(axis.text=element_text(size=15),axis.title=element_text(size=16),plot.title = element_text(size=20))
  # )
  # ggsave(paste0("Images/",gsub(" ","_",curr_title),".png"),width=20,height=15,dpi=200,units="cm")

  output <- list()
  output$beta_seq <- beta_seq
  output$mu_seq <- mu_seq
  output$var_seq <- var_seq
  return(beta_seq)
}

plot_variable_convergence <- function(df,title,true_value,ignore_first_coef = FALSE){

  if("Class" %in% colnames(df)){
    max_class = max(df$Class)
    # Can be parallelized
    for (curr_class in 0:max_class){

      subset_df <- subset(df,Class==curr_class,select=-c(Class))
      plot_variable_convergence(subset_df,paste(title,"Class", curr_class),true_value[,curr_class+1],ignore_first_coef)
    }

  }else{
    if (ignore_first_coef){
      df <- df[-1]
      if(is.numeric(true_value) || is.data.frame(true_value)){
        true_value <- true_value[-1]
      }
    }

    col_names <- colnames(df)[1:(length(colnames(df))-2)]
    plot_df <- df %>% gather(key="coord", value="Value", col_names)
    if(length(col_names) > 1){
      plot_df$coord <- factor(plot_df$coord, levels=col_names)
    }

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

}

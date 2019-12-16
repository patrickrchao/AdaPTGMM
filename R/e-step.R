library(glmnet)
expectation_gamma<- function(data,est_params,params){

  prob <- calculate_probabilities(data,est_params,params)
  numerator <- prob$expit_prob * (prob$small_prob_alt + prob$big_prob_alt)
  denominator <- numerator + (1-prob$expit_prob) * (prob$small_prob_null + prob$big_prob_null)
  expected_gamma <- numerator/denominator
  return(expected_gamma)
}

fit_beta <- function(x,gammas){

  logistic_glm_data <- data.frame(x,gammas)
  x_names <- paste("X",0:(ncol(x)-1),sep="")
  colnames(logistic_glm_data) <- c(x_names,"gammas")
  # -1 to ignore intercept
  formula <- paste("gammas ~ ",paste0(colnames(logistic_glm_data)[c(1:ncol(x))],collapse=" + ")," -1",sep="")

  glm_beta <- suppressWarnings(glm(formula,data=logistic_glm_data,family=binomial()))
  beta <- glm_beta$coefficients
  return(beta)
}

# fit_beta_ridge <- function(x,gammas){
#
#   lambdas <- 10^seq(-8, -9, by = -.7)
#
#   glm_beta <- cv.glmnet(x[,-1], cbind(gammas,1), alpha = 0, family="binomial",lambda = lambdas,intercept=TRUE)
#
#   #glm_beta <- suppressWarnings(glm(formula,data=logistic_glm_data,family=binomial()))
#   beta <- as.matrix(coef(glm_beta))
#   row.names(beta)[1] <- "X0"
#   return(beta)
# }

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

    gammas <- expectation_gamma(data,est_params,params)
    curr_beta <- fit_beta(data$full_x,gammas)


    w_ia <- calculate_w(data,est_params,params)

    if(params$testing_interval){
      for(iter in seq(5)){


        gradients <- calculate_gradients(data, est_params,params,w_ia)
        # We use gradient ascent since we want to maximize the log likelihood
        est_params$mu <- est_params$mu + 0.5 * gradients$mu
        est_params$var <- max(est_params$var + 0.5 * gradients$var,0.1)
        #print(paste(est_params$mu,est_params$var))
      }
    }else{
      est_params$mu <- weighted_mean(w_ia$w_ia,w_ia$z)
      est_params$var <- weighted_mean(w_ia$w_ia,(w_ia$z-est_params$mu)^2)
    }
    est_params$beta <- curr_beta
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
    mu_guess = sample(2:5,size=1) #true_mu#
    var_guess = sample(3:6,size=1) #true_var
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
            ggtitle(curr_title))
  }else{
    print(plot_beta_seq %>% ggplot(aes(x=Iteration,y=Value,color=Trial))+geom_line()+theme(legend.position = "none")+facet_wrap(~coord_new)+ggtitle(curr_title))

  }
  ggsave(paste0("Images/",gsub(" ","_",curr_title),".png"),width=20,height=15,dpi=200,units="cm")
  beta_seq$Trial <- as.numeric(beta_seq$Trial)

  # Mu convergence plot

  plot_mu_seq <- mu_seq %>% gather(key="coord",value="Value","Mu")
  curr_title <-  paste0(title," Mu Convergence")
  if(length(unknown)>0){
    print(plot_mu_seq %>% ggplot(aes(x=Iteration,y=Value,color=Trial))+geom_line()+theme(legend.position = "none")+geom_hline(yintercept = true_mu,linetype="dotted",size=1,color="black")+ggtitle(curr_title))
  }else{
    print(plot_beta_seq %>% ggplot(aes(x=Iteration,y=Value,color=Trial))+geom_line()+theme(legend.position = "none")+facet_wrap(~coord_new)+ggtitle(curr_title))
  }
  ggsave(paste0("Images/",gsub(" ","_",curr_title),".png"),width=20,height=15,dpi=200,units="cm")
  # var convergence plot
  plot_var_seq <- var_seq %>% gather(key="coord",value="Value","Var")
  curr_title <-  paste0(title," Variance Convergence")
  if(is.numeric(true_var)){
    print(plot_var_seq %>% ggplot(aes(x=Iteration,y=Value,color=Trial))+geom_line()+theme(legend.position = "none")+geom_hline(yintercept = true_var,linetype="dotted",size=1,color="black")+ggtitle(curr_title))
  }
  ggsave(paste0("Images/",gsub(" ","_",curr_title),".png"),width=20,height=15,dpi=200,units="cm")
  #ggsave(paste0("Images/",curr_title,".png"),width=20,height=15,dpi=800,units="cm")

  actual <- expit(data$full_x%*% true_beta)
  final_betas <-filter(beta_seq,Iteration ==params$iterations)[,c(-(ncol(beta_seq)-1),-ncol(beta_seq))]
  predictions <- expit(data$full_x %*% t(data.matrix(final_betas)))
  plot_splines <-data.frame(data$x,predictions,actual) %>% gather(key="Trial",value="gamma_mean",-data.x,-actual)

  curr_title <-  paste0(title," Spline Estimation")
  print(plot_splines%>% ggplot(aes(x=data.x,y=gamma_mean,color=Trial),size=1)+geom_line()+geom_line(aes(x=data.x,y=actual),linetype="dotted",color="black",size=1)+ggtitle(curr_title))
  ggsave(paste0("Images/",gsub(" ","_",curr_title),".png"),width=20,height=15,dpi=200,units="cm")
  #browser()
  return(beta_seq)
}

# this function returns dataframe
# one dimensional case only contains
# w_{i,big}, w_{i,small} (k assumed to be one)
calculate_w <- function(data,est_params,params){

  prob <- calculate_probabilities(data,est_params,params)


  mu <- est_params$mu
  var <- est_params$var
  beta <- est_params$beta
  small_z = data$small_z
  big_z = data$big_z

  # P[p_small|gamma=1]
  #prob_small_1 <- gaussian_pdf(small_z,mu,var)/gaussian_pdf(small_z,0,1)
  # P[p_big|gamma=1]
  #prob_big_1 <- gaussian_pdf(big_z,mu,var)/gaussian_pdf(big_z,0,1)

  # P[z_small|gamma = 0]
  #prob_small_0 <- 1#gaussian_pdf(small_z,0,1)
  # P[z_big|gamma = 0]
  #prob_big_0 <- 1#gaussian_pdf(big_z,0,1)

  prob_class_1 <- prob$expit_prob
  prob_class_0 <- 1-prob_class_1
  output <- data.frame()

  names <- c("w_ia","z","a")

  # w_{i,small}
  #numerator <- prob_class_1*prob_small_1

  numerator <- prob_class_1 * prob$small_prob_alt

  #denominator <- prob_class_1*(prob_small_1+prob_big_1/params$zeta)+prob_class_0*(prob_small_0+prob_big_0/params$zeta)
  denominator <- prob_class_1 * (prob$small_prob_alt + prob$big_prob_alt) + prob_class_0 * (prob$small_prob_null + prob$big_prob_null)
  w_small <- data.frame(numerator/denominator,small_z,"small")#rep("small",nrow(data$full_x)))
  colnames(w_small) <- names

  # w_{i,big}
  numerator <- prob_class_1*prob$big_prob_alt#/params$zeta
  w_big <- data.frame(numerator/denominator,big_z,"big")#,nrow(data$full_x)))

  colnames(w_big) <- names
  output <- rbind(w_small,w_big)
  return(output)
}

calculate_gradients <- function(data,est_params,params,w_ia){
  gradients <- list()
  mu <- est_params$mu
  var <- est_params$var
  z <- w_ia$z
  #temp_term <- exp(4*est_params$mu * w_ia$z)

  coth_term <- coth(mu*z/var)
  gradients$mu <- mean(w_ia$w_ia / var* (z*coth_term-mu))
  #mean(2*w_ia$w_ia*(-est_params$mu+w_ia$z*(temp_term+1)/(temp_term-1)))

  #pdf_1 <- gaussian_pdf(w_ia$z,est_params$mu,est_params$var)
  #pdf_2 <- gaussian_pdf(-w_ia$z,est_params$mu,est_params$var)
  #gradients$var <- mean(w_ia$w_ia/(2*est_params$var)*
  #  ( (w_ia$z-est_params$mu)^2*pdf_1-(w_ia$z+est_params$mu)^2*pdf_2/(pdf_1-pdf_2)/(est_params$var) - 1)
  #  )
  gradients$var <- mean(w_ia$w_ia/(2*var^2)*
                          (mu^2-2*mu*w_ia$w_ia*coth_term-var+z^2))
  return(gradients)
}


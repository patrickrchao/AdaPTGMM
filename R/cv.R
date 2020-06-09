cv_params <- function(x,interval,p_values=FALSE,z=FALSE,alpha_m,zeta,lambda,spline,tent,selection="BIC"){
  if(selection=="CV"){
    train_percent <- 0.6
    iterations <- 15
    num_training <- round(length(x)*train_percent)
    training_indices <- sample(seq(1,length(x)), num_training)
    train_x <- x[training_indices]
    valid_x <- x[-training_indices]

    if(interval){
      train_z <- z[training_indices]
      valid_z <- z[-training_indices]
    }else{
      train_p <- p_values[training_indices]
      valid_p <- p_values[-training_indices]
    }
  }else{
    train_x <- x
    if(interval){
      train_z <- z
    }else{
      train_p <- p_values
    }
  }

  df_range <- c(2,3,4,5)
  class_range <- c(2,3,4)
  grid <-  data.frame(matrix(0,nrow=length(df_range),ncol=length(class_range)))#expand.grid(df=df_range,class=class_range)
  colnames(grid) <- paste("class",class_range,sep="_")
  rownames(grid) <- paste("df",df_range,sep="_")

  for(df in df_range){
    for(class in class_range){
      if(interval){
        pretrain_model <- create_model_interval(train_x,train_z,df,iterations=iterations,alpha_m = alpha_m,zeta = zeta,lambda=lambda,tent=tent,num_classes=class,intervals=c(-1,1))
      }else{
        pretrain_model <- create_model(train_x,train_p,df,iterations=iterations,alpha_m = alpha_m,zeta = zeta,lambda=lambda,tent=tent,num_classes=class)
      }

      trained_model <- fit_parameters(pretrain_model)

      if(selection == "CV"){

        if(interval){
          valid_model <- create_model_interval(valid_x,valid_z,df,iterations=iterations,alpha_m = alpha_m,zeta = zeta,lambda=lambda,tent=tent,num_classes=class,intervals=c(-1,1))
        }else{
          valid_model <- create_model(valid_x,valid_p,df,iterations=iterations,alpha_m = alpha_m,zeta = zeta,lambda=lambda,tent=tent,num_classes=class)
        }

        valid_model$params <- trained_model$params
        log_lik <- log_likelihood(valid_model,verbose=FALSE)

        log_lik <- log_likelihood(valid_model,verbose=FALSE)
        grid[paste0("df_",df),paste0("class_",class)] <- grid[paste0("df_",df),paste0("class_",class)] + log_lik
      }else if(selection == "BIC"){
        log_lik <- log_likelihood(trained_model,verbose=FALSE)
        N <- length(train_x)
        # df degrees of freedom in beta
        # (class-1)*2 degrees of freedom in mu and tau^2+1
        d <- df+(class-1)*2

        grid[paste0("df_",df),paste0("class_",class)] <- log_lik - log(N)*d/2
      }else{
        stop("Invalid cross validation metric. Terminating now")
      }

    }
  }

  best_params <- which(grid == max(grid), arr.ind = TRUE)
  cv_params <- list(num_df=df_range[best_params[1]],num_classes=class_range[best_params[2]] )
  return(cv_params)
}

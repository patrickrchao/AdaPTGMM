plot_nn_masking <- function(res,x,z=NULL,pvals,alpha,title,grid_density=100){
  x <- data.frame(x=x)
  orig_df <- data.frame(x=x,pvals=pvals)
  if(!("model" %in% names(res))){
    stop("Missing model outputs. Must specify `return_all_models=TRUE` in AdaPTGMM.")
  }
  model <- .select_model(res$model,alpha,res$alphas)
  model$data$mask <- (model$data$a != "NONE")
  data <- model$data
  x_range <- seq(min(data$x),max(data$x),length.out=grid_density)
  y_range <- seq(1e-8,1-1e-8,length.out=grid_density)

  grid <- expand.grid(x=x_range, pvals=y_range)
  full_x <- data.frame(x=c(x$x,grid$x))
  full_pvals <- c(pvals,grid$pvals)

  model$args$n <- nrow(full_x)
  data <- construct_data(x=full_x,p=full_pvals,args=model$args,z=model$args$p_to_z(full_pvals))

  data$class_prob <- class_prob(model$params$beta,model$args$nclasses,x=data$x)
  model$data <- data
  w_ika <- e_step_w_ika(model,include_z=FALSE)
  odds <- big_over_small_prob(model,w_ika)
  odds_per_alpha <- res$odds_per_alpha
  # TODO: Need check that alpha exists
  odds_threshold <- odds_per_alpha[abs(res$alphas-alpha)<1e-5]

  new_mask <- data$mask
  new_mask[odds>odds_threshold] <- FALSE

  a <- data$a
  a[!new_mask] <- "NONE"
  df <- data.frame(x=data$x,pvals=data$pvals,a=a)
  col <- list(s="#EA3323",b="#1401F5",NONE=NA)
  colnames(df) <- c("x","pvals","a")
  min_x <- min(df$x)*0.99
  max_x <- max(df$x)*1.01

  orig_df$a <- a[1:nrow(orig_df)]

  outline <- data.frame(x=c(min_x,min_x,max_x,max_x),y=c(0,1,1,0))
  filepath = paste0("Images/Voronoi/Intermediate",alpha,".pdf")
  out <- ggplot(df,aes(x,pvals,fill=a))+geom_voronoi(color=NA,outline=outline)+
    geom_point(data=orig_df)+
    scale_x_continuous(limits = c(min_x,max_x), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-0.01,1.01), expand = c(0, 0))+theme(legend.position = "none")+
   labs(x=expression("predictor x"[i]),y=expression("p-value p"[i]))+ggtitle("AdaPTGMM (Intermediate Stage)")+
    scale_fill_manual(values =  col)

  ggsave(filename = filepath, plot=out,width=7.5, height=6)
}

.select_model <- function(all_models,alpha,all_alphas){
  model <- all_models
  #Check within 1e-5 since doubles are not stored exactly
  model$params <- model$params[abs(all_alphas-alpha)<1e-5][[1]]
  return(model)
}
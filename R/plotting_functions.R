plot_1d_thresh <- function(res,x,pvals,alpha,title=NULL,grid_density=200,xlab="x"){
  # Include code to load grid and lattice
  if(!("model" %in% names(res))){
    stop("Missing model outputs. Must specify `return_all_models=TRUE` in AdaPTGMM.")
  }


  # Preset Colors
  colors <-  list(s= "#FADEDE",
                  s_border = "#EA3323",
                  b = "#B5D7E4",
                  b_border = "#1401F5",
                  n = "#FFFFFF",
                  n_points ="#000000"
  )

  x_name <- colnames(x)[1]
  orig_df <- data.frame(x=x,pvals=pvals)

  colnames(orig_df)[1] <- x_name

  #Extract model and odds threshold at alpha level
  out <- .select_model_and_odds(res$model,res$odds_per_alpha,alpha,res$alphas,res$nrejs)
  model <- out$model
  odds_threshold <- out$odds_threshold
  nrejs <- out$nrejs
  if(is.null(title)){
    title <- paste("alpha =",alpha,", nrejs =",nrejs)
  }



  data <- model$data
  args <- model$args
  points_of_interest <- c(args$alpha_m, args$lambda, args$alpha_m * args$zeta + args$lambda)
  x_grid <- seq(min(data$x),max(data$x),length.out=grid_density)
  x_grid <- sort(unique(c(x_grid,points_of_interest)))
  y_grid <- seq(1e-8,1-1e-8,length.out=grid_density)
  grid <- expand.grid(x=x_grid, pvals=y_grid)

  #Combine grid and original dataset to compute which values would be masked
  full_x <- rbind(x,grid["x"])
  full_pvals <- c(pvals,grid$pvals)
  model$args$n <-  nrow(full_x)

  # Determine masking
  full_a <- .determine_masking(model,full_x,full_pvals,odds_threshold)
  orig_df$a <- full_a[1:nrow(orig_df)]
  grid$a <- full_a[(nrow(orig_df) + 1) : length(full_a)]

  # Add levels and colors to the data frames
  grid <- .append_columns(grid,colors)
  orig_df <- .append_columns(orig_df,colors)
  filepath = paste0("Images/Voronoi/Intermediate",alpha,".pdf")

  pdf("Images/Intermediate_Thresholds.pdf",width=7,height=5)
  fig <- levelplot(a_level~x*pvals,grid, at=0:3,col.regions = c(colors[["s"]],colors[["n"]],colors[["n"]],colors[["b"]]),colorkey=F,contour=T,
                   xlab = list(xlab,cex=1.5), ylab = list("p-values",cex=1.5),main=list(title,cex=1.5), scales=list(y=list(at=c(0,0.5,1))),
                   panel=function(...){
                     panel.levelplot(...)
                     grid.points(orig_df$x,orig_df$pvals, pch=16,size=unit(0.2,"char"),gp=gpar(col=orig_df$a_col))

                   },interpolate=T)
  print(fig)
  dev.off()
  return(fig)
}
.select_model_and_odds <- function(all_models,odds_per_alpha,alpha,all_alphas,nrejs){
  model <- all_models
  #Check within 1e-5 since doubles are not stored exactly
  index <-abs(all_alphas-alpha) < 1e-5
  if(!any(index)){
    stop("Invalid alpha inputted.")
  }
  model$params <- model$params[index][[1]]

  odds_threshold <- odds_per_alpha[index]

  nrejs <- nrejs[index]
  out <- list(model=model,odds_threshold = odds_threshold, nrejs = nrejs)
  return(out)
}

#' Given the model and full x/p value data frames, determine which points would be masked or unmasked
.determine_masking <- function(model,full_x,full_pvals,odds_threshold){
  n <- model$args$n

  data <- construct_data(x=full_x,p=full_pvals,args=model$args,z=model$args$p_to_z(full_pvals))
  nclasses <- model$args$nclasses
  model_type <-  model$args$model_type
  multinom_data <- head(full_x,nclasses)
  multinom_data$class <- 0:(nclasses-1)
  # Recreate the beta model
  if(model_type == "glm"){
    beta <- nnet::multinom(model$args$beta_formula, data=multinom_data,Wts=model$params$beta,maxit=0,trace=F)
  }else{
    stop(paste("plotting not implemented for model type",model$args$model_type))
  }
  data$class_prob <- class_prob(beta,nclasses,n,model_type,x=full_x)

  model$data <- data

  w_ika <- e_step_w_ika(model,include_z=FALSE)
  odds <- big_over_small_prob(model,w_ika)

  new_mask <- data$mask
  new_mask[odds>odds_threshold] <- FALSE

  a <- data$a
  a[!new_mask] <- "NONE"
  return(a)
}

.append_columns <- function(df,colors){
  df$a_level <- 0
  df$a_level[df$a == "s"] <- 1
  df$a_level[df$a == "NONE"] <- 2
  df$a_level[df$a == "b"] <- 3

  df$a_col <- colors[["n"]]
  df$a_col[df$a == "s"] <- colors[["s_border"]]
  df$a_col[df$a == "NONE"] <- colors[["n_points"]]
  df$a_col[df$a == "b"] <- colors[["b_border"]]
  return(df)
}



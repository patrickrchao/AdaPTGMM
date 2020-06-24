plot_nn_masking <- function(res,x,p,alpha,title){
  if(alpha)
  data <- model$data
  a <- data$a
  a[!data$mask] <- "NONE"
  print(sum(data$mask))
  df <- data.frame(x=data$x,pvals=data$pvals,a=a)
  # df$col <- "#000000"
  # df$col[df$a == "s"] <- "#FF0000"
  # df$col[df$a == "b"] <- "#00FF00"
  col <- list(s="#EA3323",b="#1401F5",NONE=NA)
  colnames(df) <- c("x","pvals","a")
  min_x <- min(df$x)*0.99
  max_x <- max(df$x)*1.01

  outline <- data.frame(x=c(min_x,min_x,max_x,max_x),y=c(0,1,1,0))
  filepath = paste0("Images/Voronoi/Intermediate",alpha,".pdf")
  out <- ggplot(df,aes(x,pvals,fill=a))+geom_voronoi(color=NA,outline=outline)+geom_point()+
    scale_x_continuous(limits = c(min_x,max_x), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-0.01,1.01), expand = c(0, 0))+theme(legend.position = "none")+
   labs(x=expression("predictor x"[i]),y=expression("p-value p"[i]))+ggtitle("AdaPTGMM (Intermediate Stage)")+
    scale_fill_manual(values =  col)


  ggsave(filename = filepath, plot=out,width=7.5, height=6)
}


plot_analyst_view <- function(data,args,params,reveal_threshold,t,big_odds){

  alpha_m <- args$alpha_m
  lambda <- args$lambda
  zeta <- args$zeta

  num_x <- 500
  num_y <- 200
  x_poly <- seq(min(data$x),max(data$x),length.out=num_x)


  red_y_range <- seq(1e-5,alpha_m-1e-5,length.out=num_y)
  blue_y_range <- seq(lambda+1e-5,lambda + alpha_m/zeta-1e-5,length.out=num_y)
  grid <- expand.grid(x=x_poly, p_values=c(red_y_range,blue_y_range))

  colnames(grid) <- c("x","p_values")
  #spline_x <- generate_spline(grid$x,args$num_df)


  spline_x <- args$spline_func(grid$x)#predict(args$spline_func,grid$x)
  #spline_x <- cbind(rep(1,length(grid$x)),spline_x)

  grid_data <- list()
  grid_data$x <- grid$x
  grid_data$p_values <- grid$p_values
  grid_data$full_x <- spline_x
  grid_data <- p_value_preprocessing(grid_data,args)
  grid$prob <- decision(grid_data,params,args)
  grid$color <- "BLUE"
  grid$color[grid$p_values<alpha_m] <- "RED"
  x_plot <- c(min(x_poly),x_poly,max(x_poly))

  red_poly <- grid[grid$color=="RED"& grid$prob<reveal_threshold,] %>% group_by(x) %>% summarise(p_values=max(p_values))
  red_y <- c(0,red_poly$p_values,0)
  red_x <- red_poly$x
  red_x <- c(min(red_x),red_x,max(red_x))

  if(args$tent){
    blue_poly <- grid[grid$color=="BLUE"& grid$prob<reveal_threshold,] %>% group_by(x) %>% summarise(p_values=min(p_values))
    blue_y <- c(lambda+alpha_m/zeta,blue_poly$p_values,lambda+alpha_m/zeta)

  }else{
    blue_poly <- grid[grid$color=="BLUE"& grid$prob<reveal_threshold,] %>% group_by(x) %>% summarise(p_values=max(p_values))
    blue_y <- c(lambda,blue_poly$p_values,lambda)
  }
  blue_x <- blue_poly$x
  blue_x <- c(min(blue_x),blue_x,max(blue_x))



  a <- as.factor(data$a)
  color <- rep("black",length(a))
  color[a=="s"] <- "red"
  color[a=="b"] <- "blue"
  color[!data$mask] <- "black"
  pdf(file=paste0("Images/Intermediate_",t,".pdf"))
  plot(data$x,data$p_values,col=color,
       main = "AdaPTGMM (Intermediate Stage)",ylab = expression("p-value p"[i]),xlab=expression("predictor x"[i]),
       xaxt = "n",yaxt="n",ylim=c(0,1),yaxs="i",pch=19,cex.main=2, cex.lab=1.45,cex.axis=2,cex=1.5,type="n")


  if(length(red_x) != length(red_y)){
    browser()
  }
  polygon(red_x,red_y,col = "#FADEDE",
          border = "#EA3323",
          lwd = 1)

  polygon(blue_x,blue_y,col = "#B5D7E4",
          border = "#1401F5",
          lwd = 1)
  #par(new = TRUE)
  points(data$x,data$p_values,col=color,
         main = "AdaPTGMM (Intermediate Stage)",ylab = expression("p-value p"[i]),xlab=expression("predictor x"[i]),type="p",
         xaxt = "n",yaxt="n",ylim=c(0,1),yaxs="i",pch=19,cex.main=2, cex.lab=1.45,cex.axis=2,cex=1.5,add=TRUE)

  ytick<-c(0,alpha_m,lambda,lambda+alpha_m/zeta,1)
  axis(side=2, at=ytick, labels = TRUE)


  dev.off()



  appended_x <- c(data$x,data$x[data$mask&data$a=="s"],data$x[data$mask&data$a=="b"])
  appended_p_values <- c(data$p_values,data$big_p_values[data$mask & data$a=="s"],data$small_p_values[data$mask & data$a=="b"])

  color <- rep("black",length(data$x))
  color[data$mask] <- "purple"
  color <- c(color,rep("purple",sum(data$mask)))
  pch_values  <- rep(19,length(data$x))
  pch_values[data$mask] <- 1
  pch_values <- c(pch_values,rep(1,sum(data$mask)))

  pdf(file=paste0("Images/Analyst_",t,".pdf"))
  plot(appended_x,appended_p_values,col=color,
       main = "AdaPTGMM (Analyst View)",ylab = expression("p-value p"[i]),xlab=expression("predictor x"[i]),type="n",
       xaxt = "n",yaxt="n",ylim=c(0,1),yaxs="i",pch=pch_values,cex.main=2, cex.lab=1.45,cex.axis=2,cex=1.5)

  polygon(red_x,red_y,col = "#FADEDE",
          border = "#EA3323",
          lwd = 1)

  polygon(blue_x,blue_y,col = "#B5D7E4",
          border = "#1401F5",
          lwd = 1)
  points(appended_x,appended_p_values,col=color,
         main = "AdaPTGMM (Analyst View)",ylab = expression("p-value p"[i]),xlab=expression("predictor x"[i]),type="p",
         xaxt = "n",yaxt="n",ylim=c(0,1),yaxs="i",pch=pch_values,cex.main=2, cex.lab=1.45,cex.axis=2,cex=1.5)

  ytick<-c(0,alpha_m,lambda,lambda+alpha_m/zeta,1)
  axis(side=2, at=ytick, labels = TRUE)
  dev.off()
}

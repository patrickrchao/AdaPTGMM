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

plot_masking_shape <- function(alpha_m,lambda,zeta,title=""){
  set.seed(1)

  cov <- c(runif(50),runif(7,0.5,0.8),runif(3,0,0.2),runif(10,0,1))
  pvals <- c(runif(50),rep(1e-2,10),runif(7,0,0.1),runif(3,0,0.3))

  df <- data.frame(x=cov,pvals=pvals,color="#000000")
  df$color <- as.character(df$color)
  df <- df[order(cov),]

  red <- "#FADEDE"
  red_border <- "#EA3323"
  blue <- "#B5D7E4"
  blue_border <- "#1401F5"
  pcts <- c(95,75,20)
  out <- lapply(pcts,function(x)list(1))
  n <- length(cov) # at least 2
  #func varies from 0 to 1.5

  func <- function(x){(5/(x+2)+exp(-(x-0.6)^2/(2*0.03))-1.5)*0.7-0.002}
  values <- func(df$x)
  lend <- func(0)
  rend <- func(1)
  for(ind in 1:length(pcts)){

    mask_pct <- pcts[ind]
    x_values <- c(1,0,0,df$x,1)
    y_red <- alpha_m*mask_pct/100*c(0,0,lend,values,rend)
    y_blue <- lambda+alpha_m*zeta - y_red*zeta

    red_poly <- data.frame(x=x_values,y=y_red,color=red,border =red_border)
    blue_poly <- data.frame(x=x_values,y=y_blue,color=blue,border = blue_border)

    n_poly <- length(y_red)
    df$color <- "#000000"
    df$color[df$pvals<y_red[4:(n_poly-1)]] <- red_border
    df$color[df$pvals>y_blue[3:(n+2)] & df$pvals < (alpha_m*zeta + lambda)] <- blue_border

    out[[ind]] <- ggplot(data=df,aes(x=x,y=pvals,color=color)) +
      geom_polygon(data=red_poly,aes(x=x,y=y,fill=color,color=border))+
      geom_polygon(data=blue_poly,aes(x=x,y=y,fill=color,color=border))+
      geom_point() + scale_fill_identity()+scale_color_identity()+
      theme_classic()+
      scale_x_continuous(limits=c(0,1),expand=c(0,0),breaks=NULL) +
      scale_y_continuous(limits=c(0,1),expand=c(0,0),breaks=c(0,1)) +
      labs(x = expression("predictor x"[i]), y=expression("p-value p"[i]))+
      theme(
        plot.title = element_text(size=18,color="black"),
        axis.title.x = element_blank(),#element_text(size=18, face="bold",color="black"),
        axis.title.y = element_text( size=18, face="bold",color="black"),
        axis.text.x = element_text(size=16,color="black"),
        axis.text.y = element_text(size=16,color="black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)
      )+
      coord_cartesian(xlim = c(0,1),ylim=c(0,1), # This focuses the x-axis on the range of interest
                      clip = 'off')

    if(ind >1){
      out[[ind]] <- out[[ind]] +
        theme(

        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        plot.margin = unit(c(0.5,0.5,0.3,0.5), "cm"))
    }else{
      if(alpha_m!=0.5 & lambda !=0.5){
      out[[ind]] <- out[[ind]] +
        annotate(geom="text", y=lambda, x=0, label=expression(lambda),hjust=1.3,size=6)+
        annotate(geom="text", y=alpha_m, x=0, label=expression(alpha[m]),hjust=1.3,size=6)+
        annotate(geom="text", y=lambda+alpha_m*zeta, x=0, label=bquote(lambda+ zeta*alpha[m]), hjust=1.15,size=6)
      }
      out[[ind]] <- out[[ind]] +
        theme(plot.margin = unit(c(0.5,0.5,0.3,1.5), "cm"))
    }

    if(ind ==2){
      out[[ind]] <- out[[ind]] +
        theme(
          axis.title.x = element_text(size=18, face="bold",color="black"),
           plot.title = element_text(hjust = 0.5)
      )+labs(title = title)
    }

    if(alpha_m != 0.5 & lambda != 0.5){
      out[[ind]]<- out[[ind]] +
        geom_segment(aes(y = alpha_m,yend =alpha_m, x=0,xend=1),linetype="dashed",color=red_border,size=0.5) +
        geom_segment(aes(y = lambda,yend =lambda, x=0,xend=1),linetype="dashed",color=blue_border,size=0.5)
    }else{
      out[[ind]] <- out[[ind]] + geom_segment(aes(y = 0.5-0.001,yend =0.5-0.001, x=0,xend=1),linetype="dashed",color=red_border,size=1) +
        geom_segment(aes(y =0.5+0.001,yend =0.5+0.001, x=0,xend=1),linetype="dashed",color=blue_border,size=0.7)
    }

   # out[[ind]] <- out[[ind]]+ geom_segment(aes(y = lambda+ alpha_m*zeta,yend =lambda+ alpha_m*zeta, x=0,xend=1),linetype="dashed",color=blue_border,size=0.5)
  }


  gA <- ggplotGrob(out[[1]])
  gB <- ggplotGrob(out[[2]])
  gC <- ggplotGrob(out[[3]])
  grid::grid.newpage()

  ## Initiate writing to PDF file
  filepath = paste0("Images/Masking_figures/",title,"_shape.pdf")
  #filepath = paste0("Images/Masking_figures/adapt_",mask_pct,".pdf")

  pdf(filepath, height = 4, width = 12)

  ## Create a graphical object g here
  full <- grid::grid.draw(cbind(gA, gB, gC))

  ## Stop writing to the PDF file
  dev.off()
  # library(cowplot)
  # full <- plot_grid(out[[1]], out[[2]],out[[3]], ncol=3, align="v")

  #full <- grid.arrange(gA,, ncol=length(out))




    #
    #
    # red_poly <- data.frame(x=x_values,y=y_red_gmm,color="#FADEDE",border ="#EA3323")
    # blue_poly <- data.frame(x=x_values,y=y_blue_gmm,color="#B5D7E4",border="#1401F5")
    # df <- data.frame(x=x,pvals=pvals)
    # out <- ggplot(data=df,aes(x=x,y=pvals)) +
    #   geom_polygon(data=red_poly,aes(x=x,y=y,fill=color,color=border))+
    #   geom_polygon(data=blue_poly,aes(x=x,y=y,fill=color,color=border))+
    #   geom_point() + scale_fill_identity()+scale_color_identity()+
    #   theme_classic()+
    #   scale_x_continuous(limits=c(0,1),expand=c(0,0),breaks=NULL) +
    #   scale_y_continuous(limits=c(0,1),expand=c(0,0),breaks=c(0,alpha_m,lambda,1)) +
    #   labs(x = expression("predictor x"[i]), y=expression("p-value p"[i]),title = "AdaPTGMM Masking")+
    #   theme(
    #     plot.title = element_text(size=14),
    #     axis.title.x = element_text(size=14, face="bold"),
    #     axis.title.y = element_text( size=14, face="bold"),
    #     panel.border = element_rect(colour = "black", fill=NA, size=1)
    #   )
    #
    #
    # filepath = paste0("Images/Masking_figures/gmm_",mask_pct,".pdf")
    # ggsave(filepath,plot=out,width=7.5,height=6)
  #}
}


plot_masking_function <- function(alpha_m,lambda,zeta,tent=TRUE,title="Masking_Function"){
  breakpoints <- c(0,alpha_m,lambda,lambda + alpha_m*zeta,1)

  breaks.minor <- c(alpha_m,lambda,lambda+alpha_m*zeta) #defines the edges of the categories (second label set I need)
  labels.minor <- c("alpha[m]","lambda","lambda20")
   # c(expression(alpha[m]),expression(lambda),expression(lambda),expression(lambda))#
    #
  red <- "#FADEDE"
  red_border <- "#EA3323"
  blue <- "#B5D7E4"
  blue_border <- "#1401F5"

  output <- ggplot() +
    theme_classic()+
    scale_x_continuous(limits=c(0,1),expand=c(0,0),breaks=c(0,alpha_m,lambda,lambda+alpha_m*zeta,1)) +
    scale_y_continuous(limits=c(0,1),expand=c(0,0),breaks=c(0,alpha_m,lambda,lambda+alpha_m*zeta,1)) +
    labs(x = expression("p-value p"[i]), y=bquote("Masked " ~ tilde(p)[i]),title = title)+
  theme(
    plot.title = element_text(size=18),
    axis.title.x = element_text(size=20, face="bold",vjust=1),
    axis.title.y = element_text( size=20, face="bold",vjust=3),
    #axis.text = element_text( size=16, face="bold"),
    panel.border = element_rect(colour = "black", fill=NA, size=1),
    plot.margin = unit(c(0.5,0.5,1,1.5), "cm"),
    axis.text.y = element_blank(),
    axis.text.x = element_blank()
  )+
    annotate("rect", xmin=breakpoints[1], xmax=breakpoints[2], ymin=0, ymax=1,alpha=1,fill=red)+
    annotate("rect", xmin=breakpoints[3], xmax=breakpoints[4], ymin=0, ymax=1,alpha=1,fill=blue)+
    scale_fill_identity()+scale_color_identity()+
    geom_point(alpha=0.7)+
    geom_segment(aes(x = breakpoints[4],xend =breakpoints[4], y=0,yend=1),linetype="solid",color=blue_border,size=1) +


    geom_segment(aes(x = breakpoints[1],xend =breakpoints[3], y=breakpoints[1],yend=breakpoints[3]),size=1.5)+

    geom_segment(aes(x = breakpoints[4],xend =breakpoints[5], y=breakpoints[4],yend=breakpoints[5]),size=1.5)+
    theme(text = element_text(size=25))+
    coord_cartesian(xlim = c(0,1),ylim=c(0,1), # This focuses the x-axis on the range of interest
                   clip = 'off')+
    annotate(geom="text", y=0, x=0, label="0.0",hjust=1.7, size=6)+
    annotate(geom="text", y=0, x=0, label="0.0",vjust=1.7, size=6)+
    annotate(geom="text", y=1, x=0, label="1.0",hjust=1.7, size=6)+
    annotate(geom="text", y=0, x=1, label="1.0",vjust=1.7, size=6)
    #geom_text(aes(x = breaks.minor, label = as.expression(labels.minor), y = 0),
    #          size = 4, col = 'black',vjust=2)
  #data = data.frame(br = breaks.minor,lab=labels.minor)
  if(tent){
    output <- output + geom_segment(aes(x = breakpoints[3],xend =breakpoints[4], y=breakpoints[2],yend=breakpoints[1]),size=1.5)
  }else{
    output <- output + geom_segment(aes(x = breakpoints[3],xend =breakpoints[4], y=breakpoints[1],yend=breakpoints[2]),size=1.5)
  }

    if(alpha_m != 0.5 & lambda != 0.5){
      output <- output +
        annotate(geom="text", x=lambda, y=0, label=expression(lambda),vjust=1.7,size=6)+
        annotate(geom="text", x=alpha_m, y=0, label=expression(alpha[m]),vjust=1.7,size=6)+
        annotate(geom="text", y=lambda, x=0, label=expression(lambda),hjust=1.7,size=6)+
        annotate(geom="text", y=alpha_m, x=0, label=expression(alpha[m]),hjust=1.4,size=6)+
        geom_segment(aes(x = breakpoints[2],xend =breakpoints[2], y=0,yend=1),linetype="solid",color=red_border,size=1) +
        geom_segment(aes(x = breakpoints[3],xend =breakpoints[3], y=0,yend=1),linetype="solid",color=blue_border,size=1)
    }else{
        output <- output + geom_segment(aes(x = breakpoints[2]-0.001,xend =breakpoints[2]-0.001, y=0,yend=1),linetype="solid",color=red_border,size=1) +
        geom_segment(aes(x = breakpoints[3]+0.001,xend =breakpoints[3]+0.001, y=0,yend=1),linetype="solid",color=blue_border,size=0.7)
    }
  if(zeta*alpha_m+lambda < 1){
    output <- output +
      annotate(geom="text", x=lambda+alpha_m*zeta, y=0, label=bquote(lambda+ zeta*alpha[m]), vjust=1.7,size=6)+
      annotate(geom="text", y=lambda+alpha_m*zeta, x=0, label=bquote(lambda+ zeta*alpha[m]), hjust=1.15,size=6)
  }






  filepath = paste0("Images/Masking_figures/function_",zeta,".pdf")
  ggsave(filepath,plot=output,width=7.5,height=6)
}

# plot_masking_function(0.2,0.3,3,title="AdaPTGMM Masking")
# plot_masking_function(0.5,0.5,1,title="AdaPT Masking")
#
# plot_masking_shape(0.2,0.3,3,title="AdaPTGMM")
# plot_masking_shape(0.5,0.5,1,title="AdaPT")


library(reshape2)
plot_masking_function <- function(alpha_m,lambda,zeta,tent,title="Masking_Function"){
  # tent <- args$tent
  # alpha_m <- args$alpha_m
  # lambda <- args$lambda
  # zeta <- args$zeta
  ggthemr('fresh')
  breakpoints <- c(0,alpha_m,lambda,lambda + alpha_m/zeta,1)
  output <- ggplot() +
          labs(title = title)+
          xlab("p-value")+ylab("Masked p-value") + xlim(0,1)+ylim(0,1)+
          annotate("rect", xmin=breakpoints[1], xmax=breakpoints[2], ymin=0, ymax=1,alpha=0.1)+
          annotate("rect", xmin=breakpoints[3], xmax=breakpoints[4], ymin=0, ymax=1,alpha=0.1)+
          geom_point(alpha=0.7)+
          geom_segment(aes(x = breakpoints[2],xend =breakpoints[2], y=0,yend=1),linetype="dotted",color="tomato",size=1.5) +
          geom_segment(aes(x = breakpoints[3],xend =breakpoints[3], y=0,yend=1),linetype="dotted",color="tomato",size=1.5) +
          geom_segment(aes(x = breakpoints[4],xend =breakpoints[4], y=0,yend=1),linetype="dotted",color="tomato",size=1.5) +


          geom_segment(aes(x = breakpoints[1],xend =breakpoints[3], y=breakpoints[1],yend=breakpoints[3]),size=1.5)+

          geom_segment(aes(x = breakpoints[4],xend =breakpoints[5], y=breakpoints[4],yend=breakpoints[5]),size=1.5)+
          theme(text = element_text(size=25))

  if(tent){
    output <- output + geom_segment(aes(x = breakpoints[3],xend =breakpoints[4], y=breakpoints[2],yend=breakpoints[1]),size=1.5)
  }else{
    output <- output + geom_segment(aes(x = breakpoints[3],xend =breakpoints[4], y=breakpoints[1],yend=breakpoints[2]),size=1.5)
  }
  print(output)
  ggsave(paste0("Images/",gsub(" ","_",title),".png"),width=20,height=15,dpi=200,units="cm")
}

plot_x_p_value_masking <- function(data,args){
  print(ggplot(data.frame(Covariates= data$x,p_values = data$p_values),aes(Covariates,p_values)) + geom_point()+
    stat_function(fun=helper_fdr_curve_bottom,args=list(args=args,shift=min(data$x),range=max(data$x)-min(data$x)),color="blue")+
    stat_function(fun=helper_fdr_curve_top,args=list(args=args,shift=min(data$x),range=max(data$x)-min(data$x))),color="red")+
    geom_ribbon(aes(ymin=0, ymax=1))

}

helper_sqrt <- function(x){
  return(sqrt(2)/sqrt(x*10+2))
}
helper_fdr_curve_bottom <- function(x,args,shift,range){
  return(helper_sqrt(x/range-shift)*args$alpha_m)
}


helper_fdr_curve_top <- function(x,args,shift,range){
 bottom = helper_fdr_curve_bottom(x,args,shift,range)
 scaled = bottom/args$zeta
 #reflected = -(scaled-args$alpha_m/args$zeta)+args$alpha/args$zeta
 shifted = scaled + args$lambda
 return(shifted)
}


plot_theta_density <- function(models,true_model){
  random_index <- sample(length(x), 1)
  theta_range <- seq(-4,10,0.05)
  true_params <- true_model$params
  true_args <- true_model$args
  true_num_classes <- true_args$num_classes
  true_spline_x <- true_model$data$full_x[random_index, ]

  true_class_prob <- probability_from_spline(true_spline_x, true_params$beta, true_args$num_classes)
  theta_density <- data.frame("Theta"=theta_range)
  theta_density$True_Density <- rep(0,length(theta_range))
  print(paste("True class 0 prob:",true_class_prob[1]))

  for(curr_class in 1:(true_num_classes)){

    theta_density$True_Density <- theta_density$True_Density + gaussian_pdf(theta_range,true_params$mu[curr_class],true_params$var[curr_class])*true_class_prob[curr_class]

  }

  for(trial in 1:length(models)){
    model <- models[[trial]]

    params <- model$params

    num_classes <- model$args$num_classes

    x_value <- model$data$x[random_index]
    title <- paste("Mixture density at",round(x_value,3))
    spline_x <- model$data$full_x[random_index, ]


    class_prob <- probability_from_spline(spline_x, params$beta, model$args$num_classes)
    #class_prob <- class_prob/sum(class_prob)

    col_name <- paste0("Init. ",trial)

    theta_density[,col_name] <- rep(0,length(theta_range))
    for(curr_class in 1:(num_classes)){

      theta_density[,col_name] <- theta_density[,col_name]  + gaussian_pdf(theta_range,params$mu[curr_class],params$var[curr_class])*class_prob[curr_class]
      #theta_density$True_Density <- theta_density$True_Density + gaussian_pdf(theta_range,true_args$mu[curr_class],true_args$var[curr_class])*true_class_prob[curr_class]
    }



  }
  print(paste("Est class 0 prob:",class_prob[1]))
  theta_density <- theta_density %>% gather(key="Initializations",value="Density",-True_Density,-Theta)
  output <- theta_density %>% ggplot( aes(x = Theta, y = Density,col=Initializations)) +
    geom_line() + ggtitle(title) + geom_line(aes(x=Theta,y=True_Density),linetype='dotted',color="black",size=1.5)+
    theme(axis.text=element_text(size=15),axis.title=element_text(size=16),plot.title = element_text(size=20),legend.position = "none")

  # plot_splines%>% ggplot(aes(x=X,y=Class_Probability,color=Trial),size=1)+labs(x="X")+
  #   geom_line()+geom_line(aes(x=X,y=Actual),linetype='solid',color="black",size=1.5)+
  #   ggtitle(curr_title)+
  #   theme(axis.text=element_text(size=15),axis.title=element_text(size=16),plot.title = element_text(size=20)
  # output <- theta_density %>% ggplot()+geom_line( aes(x = Theta, y = Estimated_Density), color = "blue") +
  #   geom_line( aes(x = Theta, y = True_Density), color = "red")

  print(output)
  ggsave(paste0("Images/",gsub(" ","_",title),".png"),width=20,height=15,dpi=200,units="cm")

  # print("True")
  # print(true_args$mu)
  # print(true_args$var)
  # print("Est")
  # print(params$mu)
  # print(params$var)
  # print("True class prob")
  # print(true_class_prob)

  #print("Est class prob")
 # print(class_prob)
    # theme(legend.position = "none")+facet_wrap(~coord)+ggtitle(title)+
    # geom_hline(data=true_value_df,aes(yintercept = Value),linetype="dotted",size=1,color="black")+
    # theme(axis.text=element_text(size=15),axis.title=element_text(size=16),plot.title = element_text(size=20))
}


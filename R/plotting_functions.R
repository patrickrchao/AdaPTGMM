plot_masking_function <- function(data,params){
  plot(data$p_values,data$masked_p_i)
  ggthemr('fresh')
  breakpoints <- c(0,params$alpha_m,params$lambda,params$lambda + params$alpha_m/params$zeta,1)
  output <- ggplot() +
          labs(title = "Masking Function")+
          xlab("p-value")+ylab("Masked p-value") + xlim(0,1)+ylim(0,1)+
          annotate("rect", xmin=breakpoints[1], xmax=breakpoints[2], ymin=0, ymax=1,alpha=0.1)+
          annotate("rect", xmin=breakpoints[3], xmax=breakpoints[4], ymin=0, ymax=1,alpha=0.1)+
          geom_point(alpha=0.7)+
          geom_segment(aes(x = breakpoints[2],xend =breakpoints[2], y=0,yend=1),linetype="dotted",color="tomato",size=1.5) +
          geom_segment(aes(x = breakpoints[3],xend =breakpoints[3], y=0,yend=1),linetype="dotted",color="tomato",size=1.5) +
          geom_segment(aes(x = breakpoints[4],xend =breakpoints[4], y=0,yend=1),linetype="dotted",color="tomato",size=1.5) +


          geom_segment(aes(x = breakpoints[1],xend =breakpoints[3], y=breakpoints[1],yend=breakpoints[3]),size=1.5)+

          geom_segment(aes(x = breakpoints[4],xend =breakpoints[5], y=breakpoints[4],yend=breakpoints[5]),size=1.5)

  if(params$tent){
    output <- output + geom_segment(aes(x = breakpoints[3],xend =breakpoints[4], y=breakpoints[2],yend=breakpoints[1]),size=1.5)
  }else{
    output <- output + geom_segment(aes(x = breakpoints[3],xend =breakpoints[4], y=breakpoints[1],yend=breakpoints[2]),size=1.5)
  }
  print(output)
}

plot_x_p_value_masking <- function(data,params){
  print(ggplot(data.frame(Covariates= data$x,p_values = data$p_values),aes(Covariates,p_values)) + geom_point()+
    stat_function(fun=helper_fdr_curve_bottom,args=list(params=params,shift=min(data$x),range=max(data$x)-min(data$x)),color="blue")+
    stat_function(fun=helper_fdr_curve_top,args=list(params=params,shift=min(data$x),range=max(data$x)-min(data$x))),color="red")+
    geom_ribbon(aes(ymin=0, ymax=))

}

helper_sqrt <- function(x){
  return(sqrt(2)/sqrt(x*10+2))
}
helper_fdr_curve_bottom <- function(x,params,shift,range){
  return(helper_sqrt(x/range-shift)*params$alpha_m)
}


helper_fdr_curve_top <- function(x,params,shift,range){
 bottom = helper_fdr_curve_bottom(x,params,shift,range)
 scaled = bottom/params$zeta
 #reflected = -(scaled-params$alpha_m/params$zeta)+params$alpha/params$zeta
 shifted = scaled + params$lambda
 return(shifted)
}

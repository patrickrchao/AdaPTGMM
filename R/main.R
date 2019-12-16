# alpha_m=0.05,zeta=0.1,lambda=0.4,
create_model <- function(x,p_values,num_df=10,alpha_m=0.05,zeta=0.1,lambda=0.4,spline=TRUE,iterations=20,tent=FALSE){
  check_equal_length(x,p_values)

  z_to_p <- function(z) 1-pnorm(z)
  p_to_z <- function(p) -qnorm(p)

  params <- list(alpha_m=alpha_m,zeta = zeta,lambda=lambda,spline=spline,iterations=iterations,num_df=num_df,tent=tent,
                 testing_interval=FALSE,p_to_z =p_to_z,z_to_p=z_to_p)
  if(spline){
    full_x <- generate_spline(x,num_df)
  }
  data <- list(x = x, z=-qnorm(p_values), full_x=full_x,p_values=p_values,mask=TRUE)
  model <- structure(list(params=params,data=data),class="model")
  return(model)
}


create_model_interval <- function(x,z,num_df=10,alpha_m=0.05,zeta=0.1,lambda=0.4,spline=TRUE,iterations=20,tent=FALSE,
                         intervals = c(-1,1)){

  check_equal_length(x,z)


  #intervals is either a list of length 2 or a data frame
  # data frame rows corresponds to hypothesis index
  # left_end and right_end correspond to endpoints of interval null

  if(is.vector(intervals)){
    if(length(intervals)!=2){
      print(paste0("Interval does not have length two:",intervals,". Terminating now."))
      stop()
    }
    if(intervals[1]>=intervals[2]){
      print(paste0("Interval hypothesis is invalid:",intervals[1],">=",intervals[2],". Terminating now."))
      stop()
    }
    #intervals = data.frame(z=z,left_end=intervals[1],right_end=intervals[2])
    interval_center = (intervals[1]+intervals[2])/2
    radius = intervals[2]-interval_center
    centered_z = abs(z-interval_center)
    z_to_p <- function(centered_z,radius) pnorm(centered_z+radius,lower.tail=FALSE)+pnorm(-centered_z+radius)
    z_to_p_rad <- function(centered_z) z_to_p(centered_z,radius)
    p_to_z_inv <- inverse(z_to_p_rad,lower = 0)
    p_to_z <- function(centered_z) unlist(mapply(p_to_z_inv,centered_z))
    #row.names(intervals) <- seq(1,length(x))
  #}else if(nrows(intervals)==length(z) && ncols(intervals) == 2 && all(intervals[,1]<intervals[,2])){
  #  intervals["z"] = z
  }else{
    print(paste0("Intervals is an invalid vector. Terminating now."))
    stop()
  }
  params <- list(alpha_m=alpha_m,zeta = zeta,lambda=lambda,spline=spline,iterations=iterations,num_df=num_df,tent=tent,
                 interval_radius = radius,testing_interval=TRUE,z_to_p=z_to_p_rad,p_to_z=p_to_z)
  if(spline){
    full_x <- generate_spline(x,num_df)
  }
  data <- list(x = x, z=centered_z, full_x=full_x,p_values=params$z_to_p(centered_z),mask=TRUE)
  model <- structure(list(params=params,data=data),class="model")
  return(model)
}

check_equal_length <- function(x,y){
  if(length(x)!= length(y)){
    print(paste0("Inputs do not have the same length:",length(x)," and ",length(y),". Terminating now."))
    stop()
  }
}

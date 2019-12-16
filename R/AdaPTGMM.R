AdaPTGMM <- function(data,est_params,params,calc_actual_FDP=FALSE,unknown=FALSE) {
  data$mask <- TRUE
  data <- masking(data,params)
  data <- inverse_masking(data,params)

  x_names <- colnames(data$full_x)
  fdr_log <- data.frame()

  A_t = data$mask & data$p_values>params$lambda
  R_t = data$mask & data$p_values<params$alpha_m
  size_A_t = sum(A_t)
  size_R_t = sum(R_t)


  fdphat <- (size_A_t + 1) / max(size_R_t, 1)*params$zeta
  print(paste0("alpha N/A"," FDPhat: ",round(fdphat, 3)," A_t: ",size_A_t," size_R_t: ",size_R_t))

  fdphat <- 1
  count = 0
  min_fdp <- 1
  rejections <- data.frame()
  for (t in seq(0.9, 0.01, -0.01)) {

    while (min_fdp > t) {
      if (count %% 10 == 0){
        out <- fit_parameters(data,est_params,params)
        est_params$beta <- out$best_beta
        est_params$var <- out$best_var
        est_params$mu <- out$best_mu
      }
      count = count + 1
      if (size_R_t == 0) {
        break
      }

      big_odds = decision(data,est_params,params)
      big_odds = big_odds * data$mask
      max_odds = max(big_odds)

      data$mask = as.logical(data$mask * (big_odds < quantile(big_odds,.99)))
      #data$mask = as.logical(data$mask * (big_odds < max_odds- 1e-2))
      data <- masking(data,params)
      data <- inverse_masking(data,params)


      A_t = data$mask & data$p_values>params$lambda
      R_t = data$mask & data$p_values<params$alpha_m
      size_A_t = sum(A_t)
      size_R_t = sum(R_t)
      fdphat <- (size_A_t + 1) / max(size_R_t, 1)*params$zeta
      if (fdphat < min_fdp) {
        min_fdp = min(min_fdp, fdphat)
      }
    }

    print(paste0("alpha: ",t," FDPhat: ",round(fdphat, 3)," A_t: ",size_A_t,
                 " size_R_t: ",size_R_t," minfdp ",min_fdp))
    if(calc_actual_FDP){
      actual_fdp <- sum(R_t*(unknown$theta == 0)) / max(size_R_t, 1)
      curr_row <- c(size_A_t,size_R_t,t,actual_fdp)
    }else{
      curr_row <-c(size_A_t,size_R_t,t)

    }

    fdr_log <- rbind(fdr_log, curr_row)

    rejections <- rbind(rejections,c(t, R_t))
  }
  colnames(rejections) <- c("Alpha",paste("Hypo. ",1:(length(data$x)),sep=""))
  if(calc_actual_FDP){
    colnames(fdr_log) <-c("Accepted","Rejected","FDPHat","Actual_FDP")
  }else{
    colnames(fdr_log) <-c("Accepted","Rejected","FDPHat")
  }
  fdr_log$Type <- "AdaPTGMM"
  output <- list(fdr_log=fdr_log,est_params=est_params,rejections = rejections)
  return(output)
}

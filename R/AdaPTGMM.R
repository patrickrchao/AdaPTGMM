AdaPTGMM <- function(x,
                     interval=FALSE,
                     p_values=FALSE,
                     z=FALSE,
                     num_df=10,
                     alpha_m=0.05,
                     zeta=0.1,
                     lambda=0.4,
                     spline=TRUE,
                     iterations=20,
                     tent=FALSE,
                     num_classes = 2,
                     calc_actual_FDP=FALSE,
                     unknown=FALSE,
                     likelihood=FALSE,
                     alphas = seq(0.9, 0.01, -0.01), selection="BIC") {
  cv_fit_params <- cv_params(x,interval,p_values,z,alpha_m,zeta,lambda,spline,tent,selection)
  num_df <- cv_fit_params$num_df
  num_classes <- cv_fit_params$num_classes
  print(paste("Cross Validation Parameters, df:",num_df,"classes:",num_classes))
  if(interval){
    model <- create_model_interval(x,z,num_df,iterations=iterations,alpha_m = alpha_m,zeta = zeta,lambda=lambda,tent=tent,num_classes=num_classes,intervals=c(-1,1))
  }else{
    model <- create_model(x,p_values,num_df,iterations=iterations,alpha_m = alpha_m,zeta = zeta,lambda=lambda,tent=tent,num_classes=num_classes)
  }
  data <- model$data
  args <- model$args
  params <- model$params
  # Increment number of columns by 1 if log actual FDP
  fdr_log <- data.frame(matrix(ncol = 3 + calc_actual_FDP, nrow = length(alphas)))

  rejections <- data.frame(matrix(ncol = 1 + length(x), nrow = length(alphas)))
  #hist(data$p_values,breaks=30,xlab="p values",main = "Histogram of p values")
  A_t = data$mask & data$p_values>args$lambda
  R_t = data$mask & data$p_values<args$alpha_m
  size_A_t = sum(A_t)
  size_R_t = sum(R_t)

  fdphat <- (size_A_t + 1) / max(size_R_t, 1)*args$zeta
  print(paste0("alpha N/A"," FDPhat: ",round(fdphat, 3)," A_t: ",size_A_t," size_R_t: ",size_R_t))

  fdphat <- 1
  count <- 1
  min_fdp <- 1

  model$args$iterations = 20
  model <- fit_parameters(model)
  model$args$iterations <- iterations
  for (t in alphas) {
    while (min_fdp > t) {
      if (size_R_t == 0) {
        min_fdp <- (size_A_t + 1) / max(size_R_t, 1)*args$zeta
        if(calc_actual_FDP){

          actual_fdp <- 0

          curr_row <- c(size_A_t,size_R_t,min_fdp,actual_fdp)
        }else{
          curr_row <-c(size_A_t,size_R_t,min_fdp)

        }
        fdr_log[count,] <- curr_row

        rejections[count,] <-  c(t,R_t)
        break
      }
      model <- fit_parameters(model)
      #likelihood(model)
      #browser()

      data <- model$data
      params <- model$params




      big_odds = decision(data,params,args)
      if(sum(is.na(big_odds))){
        browser()
      }
      big_odds = big_odds * data$mask
      reveal_threshold <- quantile(big_odds[data$mask],.97)
      revealed_p_values <- data$p_values[data$mask & (big_odds>=reveal_threshold)]
      #print(revealed_p_values)
      if(abs(revealed_p_values-0.405)<0.005 | revealed_p_values < 0.01){
        #browser()
      }
      data$mask <- as.logical(data$mask & (big_odds < reveal_threshold))
      if((count+1) %% 1==0){

        #plot_analyst_view(data,args,params,reveal_threshold,t,big_odds)
      }



      data <- reveal(data,args)

      A_t = data$mask & data$p_values > args$lambda
      R_t = data$mask & data$p_values < args$alpha_m
      size_A_t = sum(A_t)
      size_R_t = sum(R_t)
      fdphat <- (size_A_t + 1) / max(size_R_t, 1)*args$zeta
      #print(paste(sum(data$mask),round(fdphat,3),size_A_t,size_R_t))
      #print(paste(count,min_fdp,fdphat))
      if (fdphat < min_fdp) {
        min_fdp = min(min_fdp, fdphat)
      }

      model$data <- data

    }

    #print(paste("alpha:",round(t,2),"FDPhat:",round(fdphat, 3),"A_t:",size_A_t,
     #           "Num Rejections:",size_R_t,"minfdp",round(min_fdp,4)))

    if(calc_actual_FDP){
      if(interval){
        actual_fdp <- sum(R_t*(abs(unknown$theta)<=1.5)) / max(size_R_t, 1)
      }else{
        actual_fdp <- sum(R_t*(unknown$theta <= 0)) / max(size_R_t, 1)
      }
      curr_row <- c(size_A_t,size_R_t,t,actual_fdp)
      #curr_row <- c(size_A_t,size_R_t,min_fdp,actual_fdp)
    }else{
      curr_row <- c(size_A_t,size_R_t,t)
      #curr_row <-c(size_A_t,size_R_t,min_fdp)

    }

    fdr_log[count,] <- curr_row

    rejections[count,] <-  c(t,R_t)
    count <- count + 1
  }


  colnames(rejections) <- c("Alpha",paste("Hypo. ",1:(length(data$x)),sep=""))
  if(calc_actual_FDP){
    colnames(fdr_log) <-c("Accepted","Rejected","FDPHat","Actual_FDP")
    plot(fdr_log$FDPHat[fdr_log$FDPHat<0.4],fdr_log$Actual_FDP[fdr_log$FDPHat<0.4],xlab  = "FDPHat",ylab="True FDP")
    abline(a=0,b=1)
    #hist(unknown$theta[as.logical(rejections[81,-1])])
    #print(sum(abs(unknown$theta[as.logical(rejections[81,-1])])<=1)/length(unknown$theta[as.logical(rejections[81,-1])]))
  }else{
    colnames(fdr_log) <-c("Accepted","Rejected","FDPHat")
  }
  fdr_log$Type <- "GMM"
  output <- list(fdr_log=fdr_log,params=params,rejections = rejections)
 # likelihood(model)

  return(output)
}



# AdaPTGMM_temp <- function(x,interval=FALSE,p_values=FALSE,z=FALSE,num_df=10,alpha_m=0.05,zeta=0.1,lambda=0.4,spline=TRUE,iterations=20,tent=FALSE,num_classes = 2,calc_actual_FDP=FALSE,unknown=FALSE,likelihood=FALSE) {
# # adapt <- function(x, pvals, models,
# #                   dist = beta_family(),
# #                   s0 = rep(0.45, length(pvals)),
# #                   alphas = seq(0.01, 1, 0.01),
# #                   params0 = list(pix = NULL, mux = NULL),
# #                   nfits = 20, nms = 1,
# #                   niter_fit = 10, tol = 1e-4,
# #                   niter_ms = 20, cr = "BIC",
# #                   fs = TRUE,
# #                   verbose = list(print = TRUE, fit = FALSE, ms = TRUE)
# # ){
#
#   if(interval){
#     model <- create_model_interval(x,z,num_df,iterations=iterations,alpha_m = alpha_m,zeta = zeta,lambda=lambda,tent=tent,num_classes=num_classes)
#   }else{
#     model <- create_model(x,p_values,num_df,iterations=iterations,alpha_m = alpha_m,zeta = zeta,lambda=lambda,tent=tent,num_classes=num_classes)
#   }
#   alphas <- seq(0.9, 0.01, -0.01)
#   data <- model$data
#   args <- model$args
#   params <- model$params
#   ## Check if 'pvals' is a vector of values in [0, 1]
#   if (!is.numeric(p_values) || min(p_values) < 0 || max(p_values) > 1){
#     stop("Invalid p-values")
#   }
#
#   ## Check if the size of 'x' matches that of 'pvals'
#   if (nrow(x) != length(pvals)){
#     stop("'x' must have the same rows as the length of 'pvals'")
#   }
#
#   ## Check if necessary packages are installed.
#   #check_pkgs(models)
#
#   ## When a single model is provided, set 'nms' to be NULL
#   # if (class(models) == "adapt_model"){
#   #   model <- models
#   #   nms <- NULL
#   # } else if (is.list(models)){
#   #   types <- sapply(models, class)
#   #   if (any(types != "adapt_model")){
#   #     stop("All models should be of class \'adapt_model\'.")
#   #   }
#   #   if (!is.integer(nms) || nms <= 0){
#   #     nms <- 1
#   #   } else if (nms > nfits){
#   #     warning("Model selection cannot be more frequent than model fitting. Set \'nms\' to \'nfits\'")
#   #     nms <- nfits
#   #   }
#   #   if (length(models) == 1){
#   #     model <- models[[1]]
#   #     nms <- NULL
#   #   }
#   # }
#
#   ## Create time stamps when model is fitted or model selection is performed
#   s0 <- rep(args$alpha_m,length(p_values))
#   if(args$tent){
#     nmasks <- sum(p_values <= s0) + sum(p_values >= args$lambda+(args$alpha_m-s0)/args$zeta)#sum(p_values <= s0) + sum(p_values >= 1 - s0)
#     A <- sum(p_values >= args$lambda+(args$alpha_m-s0)/args$zeta)
#   }else{
#     nmasks <- sum(p_values <= s0) + sum(p_values >= args$lambda+s0/args$zeta)
#     A <- sum(p_values >= args$lambda+s0/args$zeta)
#   }
#   R <- sum(p_values <= s0)
#   stamps <- create_stamps(nmasks, iterations, 1)
#
#   # ## Create root arguments to simplify fitting and model selection
#   # fit_args_root <- list(
#   #   x = x, pvals = pvals, dist = dist,
#   #   niter = niter_fit, tol = tol,
#   #   verbose = verbose$fit, type = Mstep_type
#   # )
#   # if (any(stamps[, 2] == "ms")){
#   #   ms_args_root <- list(
#   #     x = x, pvals = pvals, dist = dist, models = models,
#   #     cr = cr, niter = niter_ms, tol = tol,
#   #     verbose = verbose$ms, type = Mstep_type
#   #   )
#   # }
#
#   ## Initialization
#   n <- length(p_values)
#   s <- s0
#
#   minfdp <- fdp_hat(A, R, fs) # initial FDPhat
#
#   ## Remove the alphas greater than the initial FDPhat, except the smallest one among them
#   alphas <- sort(alphas)
#   if (min(alphas) >= minfdp){
#     warning("Task completed! Initial \'s0\' has guaranteed FDR control for all alpha's in \'alphas\'.")
#     alphaind <- 0
#   } else if (max(alphas) < minfdp){
#     alphaind <- length(alphas)
#   } else {
#     alphaind <- max(which(alphas <= minfdp))
#   }
#
#   m <- length(alphas)
#   nrejs_return <- rep(0, m) # number of rejections
#   s_return <- matrix(0, n, m) # threshold
#   params_return <- list() # parameters (including pix and mux)
#   model_list <- list() # all selected models
#   info_list <- list() # other information (df, vi, etc.)
#   if(args$tent){
#     reveal_order <- which((p_values > s) & (p_values <= args$lambda+(args$alpha_m-s0)/args$zeta))
#   }else{
#     reveal_order <- which((p_values > s) & (pvals <  args$lambda+s0/args$zeta))
#   }
#
#    # the order to be revealed
#   if (length(reveal_order) > 0){
#     init_pvals <- p_values[reveal_order]
#     reveal_order <- reveal_order[order(init_pvals, decreasing = TRUE)]
#   }
#   fdp_return <- c(rep(Inf, length(reveal_order)), minfdp) # fdphat along the whole path
#
#   if (m > alphaind){
#     nrejs_return[(alphaind + 1):m] <- R
#     s_return[, (alphaind + 1):m] <- s0
#   }
#   ## alphas <- alphas[1:m]
#
#   for (i in 1:(nrow(stamps) - 1)){
#     if (alphaind == 0){
#
#       print("Task completed!")
#
#       break
#     }
#
#     alpha <- alphas[alphaind]
#     # mask <- (pvals <= s) | (pvals >= 1 - s)
#     mask <- rep(TRUE, n)
#     mask[reveal_order] <- FALSE
#     nmasks <- sum(mask)
#
#     if(args$tent){
#       nmasks <- sum(p_values <= s) + sum(p_values >= args$lambda+(args$alpha_m-s)/args$zeta)#sum(p_values <= s0) + sum(p_values >= 1 - s0)
#       A <- sum(p_values >= args$lambda+(args$alpha_m-s)/args$zeta)
#     }else{
#       nmasks <- sum(p_values <= s) + sum(p_values >= args$lambda+s/args$zeta)
#       A <- sum(p_values >= args$lambda+s/args$zeta)
#     }
#     R <- sum(p_values <= s)
#
#     start <- stamps[i, 1]
#     end <- stamps[i + 1, 1]
#     nreveals <- min(nmasks, end - start) # number of hypotheses to be revealed
#     type <- stamps[i, 2] # type of model update
#
#     ## Model selection or model fitting
#     if (type == "ms"){
#       ms_args <- c(
#         list(s = s, params0 = params),
#         ms_args_root
#       )
#       ## Use "EM_mix_ms" from "EM-mix-ms.R"
#       ms_res <- do.call(EM_mix_ms, ms_args)
#       params <- ms_res$params
#       model <- ms_res$model
#       modinfo <- ms_res$info
#     } else if (type == "fit"){
#       fit_args <- c(
#         list(s = s, params0 = params, model = model),
#         fit_args_root
#       )
#       ## Use "EM_mix" from "EM-mix.R"
#       fit_res <- do.call(EM_mix, fit_args)
#       params <- fit_res$params
#       modinfo <- fit_res$info
#     }
#
#
#
#     # ## Estimate local FDR
#     # lfdr <- compute_lfdr_mix(
#     #   pmin(pvals, 1 - pvals),
#     #   dist, params, lfdr_type)
#     # ## Find the top "nreveals" hypotheses with highest lfdr
#     # lfdr[!mask] <- -Inf
#     # inds <- order(lfdr, decreasing = TRUE)[1:nreveals]
#     # reveal_order <- c(reveal_order, inds)
#     # ## Shortcut to calculate FDPhat after revealing the hypotheses one by one
#     # Adecre <- cumsum(pvals[inds] >= 1 - s[inds])
#     # Rdecre <- cumsum(pvals[inds] <= s[inds])
#     # fdp <- fdp_hat(A - Adecre, R - Rdecre, fs)
#     # fdp_return <- c(fdp_return, fdp)
#     # fdp <- pmin(fdp, minfdp)
#     # ## Calculate the current minimum FDPhat
#     # minfdp <- min(fdp)
#     A_t = data$mask & p_values > args$lambda
#     R_t = data$mask & p_values < args$alpha_m
#     size_A_t = sum(A_t)
#     size_R_t = sum(R_t)
#     fdp <- (size_A_t + 1) / max(size_R_t, 1)*args$zeta
#     #print(paste(count,min_fdp,fdphat))
#     minfdp = min(minfdp, fdp)
#
#     while (alphaind > 0){# check if lower FDR level is achieved
#       alpha <- alphas[alphaind]
#       if (any(fdp <= alpha)){
#         breakpoint <- which(fdp <= alpha)[1]
#         lfdr_lev <- lfdr[inds[breakpoint]]
#         snew <- compute_threshold_mix(dist, params, lfdr_lev, lfdr_type)
#         snew <- pmin(s, snew)
#
#         ## Sanity check to avoid rounding errors
#         tmp_pvals <- pvals[inds[1:breakpoint]]
#         tmp_pvals <- pmin(tmp_pvals, 1 - tmp_pvals)
#         tmp_inds <- which(tmp_pvals <= snew[inds[1:breakpoint]])
#         if (length(tmp_inds) > 0){
#           snew[inds[tmp_inds]] <- pmin(snew[inds[tmp_inds]], tmp_pvals[tmp_inds] - 1e-15)
#         }
#
#         s_return[, alphaind] <- snew
#
#         fdpnew <- fdp[breakpoint]
#         Rnew <- sum(pvals <= snew)
#         nrejs_return[alphaind] <- Rnew
#         if (verbose$print){
#           cat(paste0(
#             "alpha = ", alpha, ": FDPhat ",
#             round(fdpnew, 4), ", Number of Rej. ",
#             Rnew, "\n"))
#         }
#
#         alphaind <- alphaind - 1
#       } else {
#         break
#       }
#     }
#
#     if (alphaind == 0){ # check again to save computation
#       if (verbose$print){
#         cat("Task completed!\n")
#       }
#       break
#     }
#
#     ## Update s(x)
#     final_lfdr_lev <- lfdr[tail(inds, 1)]
#     snew <- compute_threshold_mix(dist, params, final_lfdr_lev, lfdr_type)
#     s <- pmin(s, snew)
#
#     ## Sanity check to avoid rounding errors
#     tmp_pvals <- pvals[inds]
#     tmp_pvals <- pmin(tmp_pvals, 1 - tmp_pvals)
#     tmp_inds <- which(tmp_pvals <= snew[inds])
#     if (length(tmp_inds) > 0){
#       s[inds[tmp_inds]] <- pmin(s[inds[tmp_inds]], tmp_pvals[tmp_inds] - 1e-15)
#     }
#   }
#
#   remain_inds <- (1:n)[-reveal_order]
#   if (length(remain_inds) > 0){
#     tmp_pvals <- pvals[remain_inds]
#     tmp_pvals <- pmin(tmp_pvals, 1 - tmp_pvals)
#     remain_reveal_order <- remain_inds[order(tmp_pvals, decreasing = TRUE)]
#     reveal_order <- c(reveal_order, remain_reveal_order)
#     fdp_return <- c(fdp_return, rep(minfdp, length(remain_inds)))
#   }
#
#   rejs_return <- apply(s_return, 2, function(s){which(pvals <= s)})
#
#   qvals <- rep(1, n)
#   qvals[reveal_order] <- ifelse(pvals[reveal_order] <= s0[reveal_order],
#                                 cummin(fdp_return[1:n]), rep(Inf, n))
#
#   args <- list(nfits = nfits, nms = nms,
#                niter_fit = niter_fit, niter_ms = niter_ms,
#                tol = tol, cr = cr)
#
#   res <- structure(
#     list(nrejs = nrejs_return,
#          rejs = rejs_return,
#          s = s_return,
#          params = params_return,
#          qvals = qvals,
#          order = reveal_order,
#          alphas = alphas,
#          dist = dist,
#          models = model_list,
#          info = info_list,
#          args = args),
#     class = "adapt")
#   return(res)
# }
#
#
# create_stamps <- function(nmasks, nfits, nms){
#   fit_stamps <- c(seq(0, nmasks, floor(nmasks / nfits)))[1:nfits]
#   stamps <- data.frame(stamp = fit_stamps, type = "fit",
#                        stringsAsFactors = FALSE)
#   if (!is.null(nms)){
#     ms_inds <- seq(1, nfits, floor(nfits / nms))[1:nms]
#     stamps[ms_inds, 2] <- "ms"
#   }
#   stamps <- rbind(stamps, data.frame(stamp = nmasks, type = "end"))
#   return(stamps)
# }
#
# fdp_hat <- function(A, R, fs = TRUE){
#   (as.numeric(fs) + A) / pmax(1, R)
# }

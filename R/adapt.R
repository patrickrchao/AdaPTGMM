#' Adaptive P-Value Thresholding for Multiple Hypothesis Testing using Gaussian Mixture Models
#'
#' @description Fits a Gaussian Mixture model to the distribution of test statistics and
#' returns rejections and fitted parameters.
#'
#' @param x Dataframe of covariates TODO: (specify type)
#' @param pvals Vector of p-values (supply either pvals or test statistics)
#' @param z Vector of test statistics, required if \code{testing}='\code{interval}'.
#' @param se Vector of standard errors, if left blank when given test statistics, the standard errors are assumed to be 1.
#' @param testing The form of testing procedure, either "\code{one_sided}" or "\code{interval}". Default is "\code{one_sided}".
#' @param rendpoint Corresponds to right endpoint of null hypothesis interval. Required if \code{testing}=\code{'interval'}.
#' @param lendpoint Corresponds to left endpoint of null hypothesis interval. If interval testing and \code{lendpoint} is blank,
#' \code{lendpoint} will be assumed to be \code{-rendpoint}.
#' @param beta_formulas TODO: Fill this in
#' @param model_type Type of model used for modeling beta, options include \code{gam} and \code{glm}. Default is \code{glm}.
#' @param nclasses Vector of number of classes in Gaussian Mixture model. The vector corresponds to the possible
#' number of classes to select in the model selection procedure. Minimum number of classes is 2.
#' Note: recommended to use <5 classes. Default is c(2,3,4). The greater the number of degrees of freedom the longer it takes the EM procedure to fit, and the
#' longer the list of possible values, the longer the model selection procedure takes.
#' @param niter_fit Number of iterations of EM per model update.
#' @param niter_ms Number of iterations of EM in model selection.
#' @param nfits Number of model fitting steps.
#' @param alpha_m The maximum possible rejected p-value. We recommend \eqn{0.01\le \alpha_m \le 0.1}, default is 0.1.
#' @param zeta Controls minimum possible number of rejections. We recommend large values of zeta, \code{\zeta>1}, for situations with low numbers of rejections.
#' If \code{zeta}\ge 1/\alpha}, the desired FDR level, any number of rejections is possible.
#' @param lambda Controls where p-values are mirrored, boundary of blue region. TODO: Fix wording. We recommend \eqn{0.3\le\lambda\le 0.5}, default is 0.4
#' This is the most expensive part of the procedure, we recommend smaller number (<5) of iterations for larger problems. Default is 10.
#' @param masking_shape Controls the shape of the masking function, either "\code{tent}" or "\code{comb}" masking functions. Default is "\code{tent}".
#' @param alphas Vector of FDR levels of interest. Default is [0.01,0.02,...,0.89,0.9].
#' @param target_alpha_level Desired FDR level to optimize the procedure over, i.e.
#' @param cr Type of selection criterion in model_selection. Options include "\code{BIC}", "\code{AIC}", "\code{AICC}", \code{HIC}. Default is "\code{AIC}".
#' @param tol Positive scalar for early stopping if mu and tau do not update by more than \code{tol}.
#' @param initialization Type of initialization procedure to use. Options include "\code{kmeans}" or "\code{random}", where
#' "\code{random}" corresponds to drawing uniformly from predefined intervals. Default is "\code{kmeans}".
#' @param return_all_models Boolean, whether to return all models used at various alpha levels. Default \code{FALSE}.
#' Required \code{TRUE} for plot_nn_masking. Warning, can be expensive to store all models for large problems.
#' @param intercept_model Boolean. Include intercept only model in the model selection, default is \code{TRUE}.
#' @param verbose Boolean. Include print statements at each stage of the procedure.
#' @details
#'  The constraint on these masking function parameters is
#' \deqn{0< \alpha_m \le \lambda <\lambda+ \alpha_m\zeta\le 1.}
#' Setting \code{alpha_m} to 0.5, \code{lambda} to 0.5, \code{zeta} to 1, and \code{masking_shape} to "\code{tent}" results in the AdaPT masking function.
#'
#' @export

adapt_gmm <- function(x = NULL,
                      pvals = NULL,
                      z = NULL,
                      se = NULL,
                      testing = "one_sided",
                      rendpoint = NULL,
                      lendpoint = NULL,
                      beta_formulas = NULL,
                      model_type = "nnet",
                      nclasses = c(2,3,4),
                      niter_fit = 2,
                      niter_ms = 10,
                      nfits = 20,
                      alpha_m = NULL,
                      zeta = NULL,
                      lambda = NULL,
                      masking_shape = "tent",
                      alphas = seq(0.01, 1, 0.01),
                      target_alpha_level = NULL,
                      cr = "AIC",
                      tol = 1e-4,
                      initialization = "kmeans",
                      intercept_model = TRUE,
                      return_all_models = FALSE,
                      ){
  options(error =function(){traceback(2);if(!interactive()) quit('no', status = 1, runLast = FALSE)})

  n = nrow(x)

  masking_params <- select_masking_params(n,alpha_m,zeta,lambda,set_default_target(target_alpha_level,alphas))
  initialization <- select_initialization(masking_params, initialization, testing)
  .input_checks(x, pvals, z, se, testing, model_type,rendpoint, lendpoint, nclasses, niter_fit, niter_ms, nfits, masking_params, masking_shape, alphas,cr)

  args <- construct_args(testing,model_type,rendpoint,lendpoint,masking_params,masking_shape,niter_fit,niter_ms,nfits,n,initialization)

  data <- construct_data(x,pvals,z,se,args)

  beta_formulas <- clean_beta_formulas(beta_formulas,intercept_model)

  model <- model_selection(data,args,beta_formulas,nclasses,cr,initialization)


  data <- model$data
  args <- model$args


  n <- args$n
  n_alphas <- length(alphas)

  fdr_log <- data.frame(matrix(ncol = 3 , nrow = n_alphas))
  rejections <- data.frame(matrix(ncol = 1 + n, nrow = n_alphas))

  colnames(fdr_log) <-c("Alpha","Accepted","Rejected")
  colnames(rejections) <- c("Alpha",paste("Hypo. ",1:n,sep=""))

  pvals <- data$pvals
  values <- compute_fdphat(data,args)
  fdphat <- values$fdphat
  min_fdp <- fdphat

  sorted <- sort(alphas,decreasing=TRUE,index.return=TRUE)
  sorted_alphas <- sorted$x
  sorted_indices <- sorted$ix

  refitting_constant <- max(floor(sum(data$mask)/(nfits+1)),1)
  nrevealed <- 0
  big_odds <-  big_over_small_prob(model)
  to_reveal_order <- order(big_odds,decreasing=TRUE)
  reveal_order_index <- 1

  qvals <- rep(Inf, n)
  rejs <- rep(list(integer(0)), n_alphas)
  nrejs <- rep(0, n_alphas)
  all_params <- rep(list(), n_alphas)
  odds_per_alpha <- rep(Inf,n_alphas)

  reveal_hypo <- NULL
  rejection_set <- which(data$mask & (data$pvals <= args$alpha_m))
  A_t <-  sum(data$mask & data$pvals > args$lambda)
  R_t <-  length(rejection_set)
  fdphat <- ((A_t + 1)/args$zeta) / max(R_t, 1)
  min_fdp <- fdphat




  # while (min_fdp > 0 & R_t > 0 & sum(data$mask) > 0) {
  #
  #   #if((nrevealed %% refitting_constant) == 0 & nrevealed > 0){
  #
  #   model$data <- data
  #  # model <- model_selection(data,args,beta_formulas,nclasses,cr,initialization)
  #   model <- EM(model)
  #   big_odds <-  big_over_small_prob(model)
  #   to_reveal_order <- order(big_odds, decreasing=TRUE)[1:sum(data$mask)]
  #   #browser()
  #   future_A_t <- A_t - cumsum(data$pvals[to_reveal_order] > args$lambda)
  #   future_R_t <- R_t - cumsum(data$pvals[to_reveal_order] < args$alpha_m)
  #   future_fdphat <- (future_A_t + 1)/args$zeta / pmax(future_R_t,1)
  #     #reveal_order_index <- 1
  #   #}
  #   n_reveal <- min(refitting_constant, length(to_reveal_order))
  #   reveal_hypo <- to_reveal_order[1:n_reveal]
  #   A_t <- future_A_t[n_reveal]
  #   R_t <- future_R_t[n_reveal]
  #   data$mask[reveal_hypo] <- FALSE
  #   #fdphat <- min(future_fdphat[1:length(reveal_hypo)])
  #   for(reveal_index in seq(1:n_reveal)){
  #
  #     hypo_index <- to_reveal_order[reveal_index]
  #     if(data$a[hypo_index] == "s" | data$a[hypo_index] == "s_neg"){
  #       qvals[reveal_hypo] <- min_fdp
  #     }
  #     min_fdp <- min(min_fdp,future_fdphat[reveal_index])
  #   }
  # }
  #
  #
  #    # nrevealed <- nrevealed + 1
  #   #  reveal_order_index <- reveal_order_index + 1
  #   #  fdphat <- (A_t + 1) / args$zeta / max(R_t, 1)
  #
  #    # min_fdp <- min(min_fdp,fdphat)
  #  # }
  # #model$data <- data
  # for (index in seq(1:n_alphas)) {
  #   alpha <- sorted_alphas[index]
  #   rejection_set <- which(data$pvals < args$alpha_m & qvals < alpha)
  #   nrejs[sorted_indices[index]] <- length(rejection_set)
  #   rejs[[sorted_indices[index]]] <- rejection_set
  #   #qvals[rejection_set] <- min(min_fdp,qvals[rejection_set])
  # }
  # # if(is.null(reveal_hypo)){
  # #   odds_per_alpha[sorted_indices[index]] <- Inf
  # # }else{
  # #   odds_per_alpha[sorted_indices[index]] <- big_odds[reveal_hypo]
  # # }
  #
  # if(return_all_models){
  #   all_params[[sorted_indices[index]]] <- model$params
  # }
  # cat(paste0("alpha = " , alpha,
  #            ": FDPhat ",round(min_fdp, 4),
  #             ", Number of Rej. ",R_t,"\n"))






  for (index in seq(1:n_alphas)) {
    alpha <- sorted_alphas[index]
    while (min_fdp > alpha & R_t > 0) {

      if((nrevealed %% refitting_constant) == 0 & nrevealed > 0){

        model$data <- data
        # model <- model_selection(data,args,beta_formulas,nclasses,cr,initialization)
        em_out <- EM(model,return_w_ika = TRUE)
        model <- em_out$model
        w_ika <- em_out$w_ika
        big_odds <-  big_over_small_prob(model,w_ika)


        to_reveal_order <- order(big_odds, decreasing=TRUE)
        reveal_order_index <- 1

      }
      reveal_hypo <- to_reveal_order[reveal_order_index]
      data$mask[reveal_hypo] <- FALSE

      if(data$a[reveal_hypo] == "s" | data$a[reveal_hypo] == "s_neg"){
        qvals[reveal_hypo] <- min_fdp
        R_t <- R_t - 1

      }else{
        A_t <- A_t - 1
      }

      nrevealed <- nrevealed + 1
      reveal_order_index <- reveal_order_index + 1
      fdphat <- (A_t + 1) / args$zeta / max(R_t, 1)

      min_fdp <- min(min_fdp,fdphat)
    }
    model$data <- data
    rejection_set <- which(data$mask & (data$pvals < args$alpha_m))
    nrejs[sorted_indices[index]] <- R_t
    rejs[[sorted_indices[index]]] <- rejection_set
    qvals[rejection_set] <- min(min_fdp,qvals[rejection_set])
    if(is.null(reveal_hypo)){
      odds_per_alpha[sorted_indices[index]] <- Inf
    }else{
      odds_per_alpha[sorted_indices[index]] <- big_odds[reveal_hypo]
    }

    if(return_all_models){
      all_params[[sorted_indices[index]]] <- model$params
    }
    cat(paste0("alpha = " , alpha,
               ": FDPhat ",round(min_fdp, 4),
               ", Number of Rej. ",R_t,"\n"))

  }

  cat("Complete.\n")

  output <- list(nrejs=nrejs, rejs=rejs, params=model$params, qvals=qvals,alphas=alphas,
                 odds_per_alpha=odds_per_alpha, args = model$args)
  if(return_all_models){
    output$model <- model
    output$model$params <- all_params
  }
  return(output)
}

#' Compute FDPHat
#'
#' @param data data for model
#' @param args args for model
#'
#' @details Computes FDPHat with (A_t+1)/(max(R_t,1))
#' @return List with A_t, R_t, and estimated FDPHat
#' @noRd
compute_fdphat <- function(data,args){
  pvals <- data$pvals
  mask <- data$mask

  rejs <- mask & pvals < args$alpha_m
  A_t <-  sum(mask & pvals > args$lambda)
  R_t <-  sum(rejs)
  fdphat <- (A_t + 1) / max(R_t, 1) / args$zeta
  output <- list(A_t=A_t, R_t=R_t, rejs = which(as.logical(rejs)), fdphat=fdphat)
  return(output)
}



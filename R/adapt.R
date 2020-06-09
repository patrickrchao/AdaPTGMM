#' Adaptive P-Value Thresholding for Multiple Hypothesis Testing using Gaussian Mixture Models
#'
#' @description Fits a Gaussian Mixture model to the distribution of test statistics and
#' returns rejections and fitted parameters.
#'
#' @param x Vector of covariates TODO: (specify type)
#' @param p_values Vector of p-values (supply either p_values or test statistics)
#' @param z Vector of test statistics, TODO: not sure about how to deal with this yet
#' @param testing What form of testing procedure, either "\code{one_sided}" or "\code{interval}". Default is "\code{one_sided}".
#' @param ndf Vector of degrees of freedom in the spline basis for estimating \eqn{\beta}. Minimum degrees of freedom is 1,
#' representing no relationship between covariates and hypotheses. Note: Recommended to use <7 degrees of freedom.
#' Default is c(1,3,4). TODO: cannot use ndf=2 for now.
#' @param nclasses Vector of number of classes in Gaussian Mixture model. Minimum number of classes is 2.
#' Note: recommended to use <5 classes. Default is c(2,3,4).
#' @param niter Number of iterations per fitting of the expectation maximization algorithm.
#' @param alpha_m The maximum possible rejected p-value. We recommend \eqn{0.01\le \alpha_m \le 0.1}, default is 0.1.
#' @param zeta Controls minimum possible number of rejections, we recommend small values of zeta with low total number of samples.
#' If \code{zeta}=\eqn{\alpha} the desired FDR level, any number of rejections is possible.
#' @param lambda Controls where p-values are mirrored, boundary of blue region. TODO: Fix wording. We recommend \eqn{0.3\le\lambda\le 0.5}, default is 0.4
#' This is the most expensive part of the procedure, we recommend smaller number (<5) of iterations for larger problems. Default is 10.
#' @param masking_shape Controls the shape of the masking function, either "\code{tent}" or "\code{comb}" masking functions. Default is "\code{tent}".
#' @param alphas Vector of FDR levels of interest. Default is [0.9,0.89,...,0.02,0.01].
#' @param selection Type of selection procedure in model_selection. Options include "\code{BIC}", "\code{AIC}", "\code{cross_validation}". Default is "\code{BIC}".
#' @details
#'  The constraint on these masking function parameters is
#' \deqn{0< \alpha_m \le \lambda <\lambda+ \alpha_m/\zeta\le 1.}
#' Setting \code{alpha_m} to 0.5, \code{lambda} to 0.5, \code{zeta} to 1, and \code{masking_shape} to "\code{tent}" results in the AdaPT default masking function.
#'
#' @export

adapt_gmm <- function(x=NULL,
                      p_values=NULL,
                      z=NULL,
                      testing="one_sided",
                      ndf=c(1,3,4),
                      nclasses=c(2,3,4),
                      niter=5,
                      alpha_m=0.05,
                      zeta=0.1,
                      lambda=0.4,
                      masking_shape="tent",
                      alphas = seq(0.01, 1, 0.01), selection="BIC"){

  options(error =function(){traceback(2);if(!interactive()) quit('no', status = 1, runLast = FALSE)})
  input_checks(x,p_values,z,ndf,nclasses,niter,alpha_m,zeta,lambda,masking_shape,alphas)
  args <- construct_args(testing,alpha_m,zeta,lambda,masking_shape,niter,n=length(x))
  #hist(p_values,breaks=30,xlab="p values",main = "Histogram of p values")
  data <- construct_data(x,p_values,z,args)
  model <- model_selection(data,args,ndf,nclasses,selection)

  data <- model$data
  args <- model$args

  n <- args$n
  n_alphas <- length(alphas)

  fdr_log <- data.frame(matrix(ncol = 3 , nrow = n_alphas))
  rejections <- data.frame(matrix(ncol = 1 + n, nrow = n_alphas))

  colnames(fdr_log) <-c("Alpha","Accepted","Rejected")
  colnames(rejections) <- c("Alpha",paste("Hypo. ",1:n,sep=""))


  p_values <- data$p_values
  values <- compute_fdphat(data,args)
  fdphat <- values$fdphat
  min_fdp <- values$fdphat
  sorted <- sort(alphas,decreasing=TRUE,index.return=TRUE)

  sorted_alphas <- sorted$x
  sorted_indices <- sorted$ix

  for (index in seq(1:n_alphas)) {

    alpha <- sorted_alphas[index]
    while (min_fdp > alpha & values$R_t > 0) {
      model <- EM(model)
      big_odds = big_over_small_prob(model)
      reveal_threshold <- quantile(big_odds[data$mask],.97)

      reveal_indices <- big_odds >= reveal_threshold
      data <- reveal(data,reveal_indices)

      values <- compute_fdphat(data,args)
      fdphat <- values$fdphat
      min_fdp <- min(min_fdp,fdphat)
      model$data <- data

    }

    A_t <- values$A_t
    R_t <- values$R_t
    fdr_log[sorted_indices[index],] <- c(alpha,A_t,R_t)
    rejections[sorted_indices[index],] <- c(alpha,values$rejs)

    print(paste("alpha:",round(alpha,2),"FDPhat:",round(fdphat, 3),"A_t:",A_t,
                "Num Rejections:",R_t,"minfdp",round(min_fdp,4)))
  }
  output <- list(fdr_log = fdr_log, params=model$params, rejections = rejections)
  #plot(fdr_log$Alpha,fdr_log$Rejected)
  return(output)
}

#' Compute FDPHat
#'
#' @param data data for model
#' @param args args for model
#'
#' @details Computes FDPHat with (A_t+1)/(max(R_t,1))
#' @return List with A_t, R_t, and estimated FDPHat
compute_fdphat <- function(data,args){
  p_values <- data$p_values
  mask <- data$mask

  rejs <- mask & p_values < args$alpha_m
  A_t <-  sum(mask & p_values > args$lambda)
  R_t <-  sum(rejs)
  fdphat <- (A_t + 1) / max(R_t, 1)*args$zeta
  output <- list(A_t=A_t, R_t=R_t,rejs = as.logical(rejs), fdphat=fdphat)
  return(output)
}



#' Perform Expectation step for gamma
#'
#' We need to sum over a in the w_ika computation
#' @return gamma probability
#' P[\gamma_i = k | x_i, \tilde p_i]= \sum_a w_{ika}
#' \sum_{a} P[\gamma_i = k | x_i] P[a_i=a_ia,\tilde p_i | \gamma_i=k]/\sum_{a',k'} P[\gamma_i = k' | x_i] P[a_i=a_ia',\tilde p_i | \gamma_i=k']
#' @noRd
e_step_gamma <- function(model,w_ika){
  gammas <- subset(marginalize(w_ika,"a"),select=-c(i))
  gammas <- gammas[order(gammas$class),]

  return(gammas)
}


#' Perform Expectation step for w_ika
#'
#' @param model model class
#' @param prev_w_ika Previous w_ika table, providing this speeds up w_ika computation because we only need to update
#' the value of each w_ika
#' @param include_z boolean corresponding to whether to include the 'z' column in the output,
#' FALSE for likelihood computation
#' @param agg_over_hypotheses boolean corresponding to whether to sum over a and class,
#' TRUE for likelihood computation
#'
#' @return DataTable, columns include a, i, k, value
#' a corresponds to big/small
#' i corresponds to hypothesis number
#' k corresponds to which gaussian class
#' value is the computed probability
#'
#' Rows correspond to hypotheses
#'
#' P[gamma_i=k,a_i=a_ia | x_i, \tilde p_i]=
#' {P[a_i=a,\tilde p_i | \gamma=k]P[\gamma=k | x_i]\zeta^1{a_ia=b}}/ {\sum_{a',k'} P[a_i=a',\tilde p_i | \gamma=k']P[\gamma=k' | x_i]  \zeta^1{a_ia'=b}}
#' @noRd
e_step_w_ika <- function(model, prev_w_ika = NULL, include_z = TRUE, agg_over_hypotheses = FALSE){

  # need to iterate over a, k
  # use helper in each case
  data <- model$data
  args <- model$args
  all_a <- args$all_a
  params <- model$params
  nclasses <- args$nclasses
  num_a <- length(all_a)
  n <- args$n
  testing <- args$testing



  # if previous w_ika exists, use previous w_ika
  if(is.null(prev_w_ika)){
    # 5 columns corresponding to a, class, i, value, z
    w_ika <- data.table::data.table(matrix(0,nrow=num_a*nclasses*n,ncol=5))
    w_ika <- setNames(w_ika,c("a","class","i","value","z"))

    # Fill in w_ika table
    w_ika$i <- rep(1:n, num_a * nclasses)
    w_ika$a <- rep(all_a,each = n * nclasses)
    w_ika$class <- rep(rep(1:nclasses,each=n), num_a)

    # Add corresponding z to each row
    # Unmasked hypotheses use true z
    # Masked hypotheses use z corresponding to z_small or z_big
    masked_i <- data$mask[w_ika$i]
    unmasked_i <- !masked_i
    if(include_z){
      w_ika$z <- 0
      w_ika$z[unmasked_i]                 <- data$z      [w_ika$i[unmasked_i]]
      w_ika$z[masked_i & w_ika$a == "s"]  <- data$small_z[w_ika$i[masked_i & w_ika$a == "s"]]
      w_ika$z[masked_i & w_ika$a == "b"]  <- data$big_z  [w_ika$i[masked_i & w_ika$a == "b"]]
      w_ika$z[masked_i & (w_ika$a == "s_neg" | w_ika$a == "b_neg")]  <- -1 * w_ika$z[masked_i & (w_ika$a == "s" | w_ika$a == "b")]

    }

  }else{
    w_ika <- prev_w_ika
  }

  # Fill in w_ika values
  count <- 0

  for(a in args$all_a){
    for(k in 1:nclasses){
      start_i <- 1 + count * n
      end_i <- (1 + count) * n

      w_ika[start_i:end_i,"value"] <- w_ika_helper(a, k, data, params$mu, params$var, args$zeta, args$jacobian)
      count <- count + 1
    }
  }
  if(length(args$all_a) == 4){

    masked_i <- data$mask[w_ika$i]
    unmasked_i <- !masked_i

    # Pairwise Reveal
    subset <- (data$z > 0 & data$a == "s") | (data$z < 0 & data$a == "b")
    w_ika$value[subset & (w_ika$a == "b" | w_ika$a=="s_neg") & masked_i] <- 0

    subset <- (data$z < 0 & data$a == "s") | (data$z > 0 & data$a == "b")
    w_ika$value[subset & (w_ika$a == "b_neg" | w_ika$a=="s") & masked_i] <- 0



    # Sign Reveal
 #  w_ika$value[((w_ika$z<0 ) & masked_i)& data$z>0 ] <-  0
  # w_ika$value[((w_ika$z>0 ) & masked_i)& data$z<0 ] <-  0
  }

  # Normalize by total sum, or P[\tilde p_i | x_i]
  #sum over a and gamma and divide by the total
  w_ika <- w_ika[, value:= value/sum(value), by=i]
  #w_ika$value <- #ave(x=w_ika$value,c(w_ika$i),FUN=function(x) x/sum(x))
  return(w_ika)
}


#' Computes P[a_i=a,\tilde p_i | \gamma=k]P[\gamma=k | x_i] for each a,k
#'
#' @param a which a_i the row corresponds to
#' @param class which class the row corresponds to
#' @param data data class
#' @param mu vector of means of the Gaussian models
#' @param var vector of variances of the Gaussian models
#' @param zeta args$zeta, ratio of size of small/big region
#' @param jacobian jacobian probability function, varies for one sided vs interval null
#'
#' For a=b, the probability is multiplied by zeta
#'
#' For unmasked p-values, we compute the probability P[p_i | \gamma=k]P[\gamma=k | x_i]
#' and store it in P[a_i=s,\tilde p_i | \gamma=k]P[\gamma=k | x_i]
#' In this way, we do not include the normalization from zeta and set
#' the probability with a_i=b to zero.
#' @noRd
w_ika_helper <- function(a,class,data,mu,var,zeta,jacobian){
  mask <- data$mask
  se <- data$se

  prob <- rep(0,length(data$pvals))

  sign <- ifelse(a == "s_neg" | a == "b_neg", -1, 1)

  if(a == "b" | a == "b_neg"){
    prob[mask]  <- jacobian(sign * data$big_z[mask],   mu[class], var[class], se=se[mask]) * zeta
  }else if(a == "s_neg" | a == "s"){
    prob[mask]  <- jacobian(sign * data$small_z[mask], mu[class], var[class], se=se[mask])
    if(a == "s"){
      prob[!mask] <- jacobian(data$z[!mask], mu[class], var[class], se=se[!mask])
    }
  }

  # Scale by class probability
  prob <- prob * data$class_prob[, class]
  return(prob)
}

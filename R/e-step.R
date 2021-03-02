#' Perform Expectation step for gamma
#'
#' We need to sum over a in the w_ika computation
#' @return gamma probability
#' P[\gamma_i = k | x_i, \tilde p_i]= \sum_a w_{ika}
#' \sum_{a} P[\gamma_i = k | x_i] P[a_i=a_ia,\tilde p_i | \gamma_i=k]/\sum_{a',k'} P[\gamma_i = k' | x_i] P[a_i=a_ia',\tilde p_i | \gamma_i=k']
#' @noRd
e_step_gamma <- function(w_ika){
  gammas <- subset(marginalize(w_ika,"a"),select=-c(i))
  gammas <- gammas[order(gammas$class),]

  return(gammas)
}


#' Perform Expectation step for w_ika
#'
#' @param model model class
#' @param prev_w_ika Previous w_ika table, providing this speeds up w_ika computation because we only need to update
#' the value of each w_ika
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
e_step_w_ika <- function(model, prev_w_ika = NULL, normalize = TRUE){

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
    w_ika <- data.table::data.table(matrix(0,nrow=num_a*nclasses*n,ncol=7))
    w_ika <- setNames(w_ika,c("a","class","i","value","z","z_base_prob","class_density"))

    # Fill in w_ika table
    w_ika$i <- rep(1:n, num_a * nclasses)
    w_ika$a <- rep(all_a,each = n * nclasses)
    w_ika$class <- rep(rep(1:nclasses,each=n), num_a)

    # Add corresponding z to each row
    # Unmasked hypotheses use true z
    # Masked hypotheses use z corresponding to z_small or z_big
    masked_i <- data$mask[w_ika$i]
    unmasked_i <- !masked_i

    w_ika$z <- 0
    w_ika$z[unmasked_i]                 <- data$z      [w_ika$i[unmasked_i]]
    w_ika$z[masked_i & w_ika$a == "s"]  <- data$small_z[w_ika$i[masked_i & w_ika$a == "s"]]
    w_ika$z[masked_i & w_ika$a == "b"]  <- data$big_z  [w_ika$i[masked_i & w_ika$a == "b"]]
    w_ika$z[masked_i & (w_ika$a == "s_neg" | w_ika$a == "b_neg")]  <- -1 * w_ika$z[masked_i & (w_ika$a == "s" | w_ika$a == "b")]

    w_ika$z_base_prob <- args$base_prob(w_ika$z,data$se)

  }else{
    w_ika <- prev_w_ika
  }

  # Fill in w_ika values
  w_ika <- w_ika_helper(w_ika,args,data,params)
  # Normalize by base probability to account for jacobian
  w_ika$value <- w_ika$class_density / w_ika$z_base_prob

  if(length(args$all_a) == 4){

    masked_i <- data$mask[w_ika$i]
    unmasked_i <- !masked_i

    # Pairwise Reveal
    subset <- (data$z[w_ika$i] > 0 & data$a[w_ika$i] == "s") |
                (data$z[w_ika$i] < 0 & data$a[w_ika$i] == "b")
    w_ika$value[subset & (w_ika$a == "b" | w_ika$a=="s_neg") & masked_i] <- 0

    subset <- (data$z[w_ika$i] < 0 & data$a[w_ika$i] == "s") |
                (data$z[w_ika$i] > 0 & data$a[w_ika$i] == "b")
    w_ika$value[subset & (w_ika$a == "b_neg" | w_ika$a=="s") & masked_i] <- 0

    # Sign Reveal
 #  w_ika$value[((w_ika$z<0 ) & masked_i)& data$z>0 ] <-  0
  # w_ika$value[((w_ika$z>0 ) & masked_i)& data$z<0 ] <-  0
  }

  # Normalize by total sum, or P[\tilde p_i | x_i]
  #sum over a and gamma and divide by the total
  # Does not normalize for log likelihood computation
  if(normalize){
    w_ika <- w_ika[, value:= value/sum(value), by=i]
  }
  if(any(is.na(w_ika))){
      stop("NA value in w_ika table. Stopping.")
  }
  #w_ika$value <- #ave(x=w_ika$value,c(w_ika$i),FUN=function(x) x/sum(x))
  return(w_ika)
}





#' Computes P[a_i=a,\tilde p_i | \gamma=k]\ for each a,k without the denominator normalization
#' TODO: Fix this documentation
#' @w_ika partially filled out w_ika table
#' @param args args class
#' @param data data class
#' @param param params class
#'
#' For a=b/neg_z, the probability is multiplied by zeta
#'
#' For unmasked p-values, we compute the probability P[p_i | \gamma=k]
#' and store it in P[a_i=s,\tilde p_i | \gamma=k]
#' In this way, we do not include the normalization from zeta and set
#' the probability with a_i=b to zero.
#' @noRd
w_ika_helper <- function(w_ika,args,data,params){


  subset <- data$mask[w_ika$i] | w_ika$a == "s"
  w_ika <- w_ika[subset,]

  sign <- rep(1,sum(subset))
  zeta_jacobian <- (w_ika$a == "b" | w_ika$a == "neg_b") * (args$zeta - 1) + 1
  if(length(args$all_a) == 4){
    sign <- (w_ika$a == "s" | w_ika$a == "b") * 2 - 1
  }

  w_ika$class_density <- dnorm(w_ika$z * sign,
                                       mean = params$mu[w_ika$class],
                                       sd = sqrt(params$var[w_ika$class] + data$se[w_ika$i])) *
    zeta_jacobian * data$class_prob[matrix(c(w_ika$i,w_ika$class),ncol=2)]

  return(w_ika)
}

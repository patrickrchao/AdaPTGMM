#' Perform Expectation step for gamma
#'
#' We just need to sum over a in the w_ika computation
#' @return gamma probability
#' P[\gamma_i = k | x_i, \tilde p_i]=
#' \sum_{a} P[\gamma_i = k | x_i] P[a_i=a_ia,\tilde p_i | \gamma_i=k]/\sum_{a',k'} P[\gamma_i = k' | x_i] P[a_i=a_ia',\tilde p_i | \gamma_i=k']
#' @importFrom magrittr "%>%"
#' @noRd
e_step_gamma <- function(model,w_ika){

  gammas  <- w_ika %>%
                dplyr::select(-c(a)) %>% # remove 'a' column
                dplyr::group_by(class,i) %>% #groupby class and hypothesis
                dplyr::summarise(weight=sum(value)) %>%# sum across a_i
                #tidyr::spread(class,value) %>% #reshape into number of hypo by number of classes
                dplyr::select(-c(i)) #remove hypothesis numbering column

  return(gammas)
}


#' Perform Expectation step for w_ika
#'
#' @return Dataframe, columns include a, i, k, value
#' a corresponds to big/small
#' i corresponds to hypothesis number
#' k corresponds to which gaussian class
#' value is the computed probability
#'
#' Rows correspond to hypotheses
#'
#'
#' P[a_i=a,\tilde p_i | \gamma=k]P[\gamma=k | x_i] \zeta^1{a_ia=b}/ {\sum_{a',k'} P[a_i=a',\tilde p_i | \gamma=k']P[\gamma=k' | x_i] \zeta^1{a_ia'=b}}
#' This probability is equivalent to w_ika, or P[gamma_i=k,a_i=a_ia | x_i, \tilde p_i]
e_step_w_ika <- function(model){

  # need to iterate over a, k
  # use helper in each case
  data <- model$data
  args <- model$args
  all_a <- args$all_a
  params <- model$params
  nclasses <- args$nclasses
  num_a <- length(all_a)
  n <- args$n

  all_prob <- expand.grid(a=all_a,class=0:(nclasses-1))

  # prob <- cbind(all_prob,
  #               t(apply(all_prob,MARGIN = 1,FUN=w_ika_helper,
  #                       data = data,
  #                       mu = params$mu,
  #                       var = params$var,
  #                       zeta = args$zeta)))

  prob <- data.frame(matrix(0,nrow=num_a*nclasses*n,ncol=5))
  colnames(prob) <- c("a","class","i","value","z")
  count <- 0
  prob$i <- rep(1:n,num_a*nclasses)
  prob$a <- rep(all_a,each = n*nclasses)
  prob$class <- rep(rep(1:nclasses,each=n),num_a)
  for(a in args$all_a){
    for(k in 0:(nclasses-1)){
      start_ind <- 1+count*n
      end_ind <- (count+1)*n
      #prob[start_ind:end_ind,"a"] <- a
      #prob[start_ind:end_ind,"class"] <- k
      #prob[start_ind:end_ind,"i"] <- 1:n
      prob[start_ind:end_ind,"value"] <- w_ika_helper(a,k,data,params$mu,params$var,args$zeta)
      count <- count + 1
    }
  }

  #prob <- tidyr::gather(prob,"i","value",-a,-class)

  # Normalize by total sum, or P[\tilde p_i | x_i]
  #browser()
  #hypo_sum <- aggregate(prob$value, by=list(i=prob$i), FUN=sum)
  #prob$value <- prob$value / hypo_sum$[prob$i]

  prob <- prob %>% dplyr::group_by(i) %>%
    dplyr::mutate(value = value / sum(value))%>% # sum over a and gamma
    dplyr::ungroup()


  # Add corresponding z to each row
  # Unmasked hypotheses use true z
  # Masked hypotheses use z corresponding to z_small or z_big
  masked_indices <- data$mask[prob$i]
  unmasked_indices <- !masked_indices

  prob$z <- 0
  prob$z[unmasked_indices]                <- data$z      [prob$i[unmasked_indices]]
  prob$z[masked_indices & prob$a == "s"]  <- data$small_z[prob$i[masked_indices & prob$a == "s"]]
  prob$z[masked_indices & prob$a == "b"]  <- data$big_z  [prob$i[masked_indices & prob$a == "b"]]

  return(prob)
}


#' Computes P[a_i=a,\tilde p_i | \gamma=k]P[\gamma=k | x_i] for each a,k
#'
#' For a=b, the probability is divided by zeta
#'
#' For unmasked p-values, we compute the probability P[p_i | \gamma=k]P[\gamma=k | x_i]
#' and store it in P[a_i=s,\tilde p_i | \gamma=k]P[\gamma=k | x_i]
#' In this way, we do not include the normalization from zeta and set
#' the probability with a_i=b to zero.
#'
#'
w_ika_helper <- function(a,class,data,mu,var,zeta){
  mask <- data$mask
  # Add one since class k has parameters at index k+1 (classes begin at 0)
  class_ind <- class + 1

  prob <- rep(0,length(data$x))

  if(a=="s"){
    prob[!mask] <- prob_jacobian(data$z[!mask]    ,mu[class_ind],var[class_ind])
    prob[mask] <- prob_jacobian(data$small_z[mask],mu[class_ind],var[class_ind])
  }else{
    prob[mask] <- prob_jacobian(data$big_z[mask],  mu[class_ind],var[class_ind])/zeta
  }
  prob <- prob * data$class_prob[,class_ind]
  return(prob)
}

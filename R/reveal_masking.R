#' Compute Probability of big vs small p-value from masked data
#'
#' @param model Model class
#' @param w_ika w_ika data frame
#' @return Vector of odds ratio of big/small
#' @details Specifically, we compute P[a_i=b | \tilde p_i,x_i]/P[a_i=s | \tilde p_i,x_i]
#'
#' @noRd
big_over_small_prob <- function(model,w_ika = NULL){
  if(is.null(w_ika)){
    w_ika <- e_step_w_ika(model)
  }
  # marginalize over class and remove hypothesis numbering column
  big_small <- tidyr::spread(marginalize(w_ika,"class"),a,value)
  browser()
  if(model$args$testing == "interval"){
    odds <- (big_small$b + big_small$b_neg)/ (big_small$s + big_small$s_neg)
  }else if(model$args$testing == "one_sided"){
    odds <- (big_small$b) / (big_small$s)
  }


  # Fast Version
  # big_small <- marginalize(w_ika,"class")
  # big_small <- big_small[order(big_small$a),]
  # big <- seq(1,nrow(big_small),2)
  # small <- seq(2,nrow(big_small),2)
  # odds <- big_small$value[big]/big_small$value[small]

  return(odds)

}

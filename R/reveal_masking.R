#' Compute Probability of big vs small p-value from masked data
#'
#' @param model
#'
#' @return Vector of odds ratio of big/small
#' @details Specifically, we compute P[a_i=b | \tilde p_i,x_i]/P[a_i=s | \tilde p_i,x_i]
#'
#' @noRd
big_over_small_prob <- function(model){

  w_ika <- e_step_w_ika(model)

  # marginalize over class and remove hypothesis numbering column
  big_small <- tidyr::spread(marginalize(w_ika,"class"),a,value)
  odds <- big_small$b / big_small$s

  # Fast Version
  # big_small <- marginalize(w_ika,"class")
  # big_small <- big_small[order(big_small$a),]
  # big <- seq(1,nrow(big_small),2)
  # small <- seq(2,nrow(big_small),2)
  # odds <- big_small$value[big]/big_small$value[small]

  return(odds)

}


#' Reveal Hypotheses
#'
#' @param data data class
#' @param indices indices to reveal
#'
#' @return data class
#' @noRd
reveal <- function(data,indices){
  data$mask[indices] <- FALSE
  return(data)
}

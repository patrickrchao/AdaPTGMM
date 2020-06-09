#' Compute Probability of big vs small p-value from masked data
#'
#' @param model
#'
#' @return Vector of odds ratio of big/small
#' @details Specifically, we compute P[a_i=b | \tilde p_i,x_i]/P[a_i=s | \tilde p_i,x_i]
#'
#' @noRd
big_over_small_prob <- function(model){

  prob <- e_step_w_ika(model)

  big_small <- prob %>% dplyr::group_by(a,i) %>% #groupby a and hypothesis
    dplyr::summarise(value=sum(value))  %>% # sum across gamma
    tidyr::spread(a,value) %>%  #create rows for a=s,b
    dplyr::select(-c(i)) #remove hypothesis numbering column

  odds <- big_small$b / big_small$s

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

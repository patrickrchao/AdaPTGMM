# Returns ratio (odds) of the masked p_value belonging to the big value and the small value
# For example, the probability of the p_value being 0.9 to 0.1 (if using symmetric masking)
decision <- function(data, est_params, params) {
  prob <- calculate_conditional_probabilities(data,est_params,params)
  small_prob <- prob$expit_prob * (prob$s_1) + (1 - prob$expit_prob) * prob$s_0
  big_prob <-   prob$expit_prob * (prob$b_1) + (1 - prob$expit_prob) * prob$b_0
  return(big_prob / small_prob)
}



# Returns ratio (odds) of the masked p_value belonging to the big value and the small value
# For example, the probability of the p_value being 0.9 to 0.1 (if using symmetric masking)
decision <- function(data, est_params, params) {
  prob <- calculate_probabilities(data,est_params,params)
  small_prob <- prob$expit_prob * (prob$small_prob_alt) + (1 - prob$expit_prob) * prob$small_prob_null
  big_prob <- prob$expit_prob * (prob$big_prob_alt) + (1 - prob$expit_prob) * prob$big_prob_null
  return(big_prob / small_prob)
}



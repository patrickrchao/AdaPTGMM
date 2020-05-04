# Returns ratio (odds) of the masked p_value belonging to the big value and the small value
# For example, the probability of the p_value being 0.9 to 0.1 (if using symmetric masking)
decision <- function(data, params, args) {

  prob <- calculate_conditional_probabilities(data,params,args)
  small_prob <- 0
  big_prob <- 0

  for (class in 0:(args$num_classes - 1)) {

    for (a in args$all_a) {
      class_prob_str = paste0("class_prob_",class)
      weight_str <- paste0(a,"_",class)
      curr_prob <- prob[[class_prob_str]]*prob[[weight_str]]
      if(str_detect(a,"s")){
        small_prob <- small_prob + curr_prob
      }else{
        big_prob <- big_prob + curr_prob
      }
    }
  }

  return(big_prob / small_prob)
}



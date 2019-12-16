# Returns ratio (odds) of the masked p_value belin
decision <- function(data,est_params,params) {

  mu <- est_params$mu
  var <- est_params$var
  beta <- est_params$beta
  small_z = data$small_z
  big_z = data$big_z
  expit_prob <- expit(data$full_x %*% beta)
  small_prob <- expit_prob*(gaussian_pdf(small_z,mu,var)/
                             gaussian_pdf(small_z,0,1))+ (1-expit_prob)*1

  big_prob <- expit_prob*(gaussian_pdf(big_z,mu,var)/
                              gaussian_pdf(big_z,0,1)/params$zeta)+ (1-expit_prob)/params$zeta
  return(big_prob/small_prob)
}
